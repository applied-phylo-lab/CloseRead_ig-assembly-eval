import os
import sys
import subprocess
import logging
from .steps.data_prep import data_prep
from .steps.convert_primary_bam import convert_primary_bam
from .steps.loci_location import loci_location
from .steps.final_ig_loci import final_ig_loci
from .steps.cigar_processing import cigar_processing
from .steps.coverage_analysis import coverage_analysis
import argparse
from concurrent.futures import ThreadPoolExecutor

def resolve_path(path):
    return os.path.abspath(os.path.expanduser(path))

def run_pipeline_cli():

    parser = argparse.ArgumentParser( description="Run the CloseRead pipeline.", add_help=False)
    parser.add_argument('--help', action='help', help='Show this help message and exit')
    parser.add_argument("-s", "--species", required=True, help="Comma-separated list of species (e.g., species1,species2).")
    parser.add_argument("-d", "--home", required=True, type=resolve_path, help="Path to the working directory.")
    parser.add_argument("-h", "--haploid", required=True, help="Haploid status (True or False).")
    parser.add_argument("-fq", "--fastqdir", required=True, type=str, help="Path to the FASTQ directory.")
    parser.add_argument("-t", "--threads", required=False, default=32, type=int, help="# of threads to use (default: 32).",dest="t")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        '-ig', '--igdetective_home', 
        type=resolve_path, 
        help="Path to the IGDetective directory."
    )
    group.add_argument(
        '-cig', '--customIG', 
        type=str, 
        help="Path to directory containing ${species_name}.customIG.txt."
    )
    args = parser.parse_args()
    run_pipeline(args)


def configure_logging(home):
    """Configure logging for the pipeline."""
    log_dir = os.path.join(home, "logs")
    os.makedirs(log_dir, exist_ok=True)  # Create the logs directory if it doesn't exist

    logging.basicConfig(
        filename=os.path.join(log_dir, f"pipeline.log"),
        filemode="w",  # Overwrite log file each run
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
    )
    logging.info(f"Logging initialized. Logs will be written to {os.path.join(log_dir, 'pipeline.log')}")

def parallel_step_1_and_2(species, home, fastqdir, haploid, threads, igdetective_home=None):
    output_bam = os.path.join(home, "aligned_bam", species, f"{species}_merged_sorted.bam")
    output_index = f"{output_bam}.csi"
    run_flag = False
    if not os.path.exists(output_bam) or not os.path.exists(output_index):
        run_flag = True
    else:
        run_flag = False

    def step_1():
        if not os.path.exists(output_bam) or not os.path.exists(output_index):
            try:
                logging.info(f"Step 1: Data Preparation for {species}")
                data_prep(species, home, fastqdir, haploid, str(threads))
            except Exception as e:
                logging.error(f"Step 1 failed: {e}")
        else:
            logging.info(f"Step 1: Output exists, skipping Data Preparation for {species} because it's already completed.")        

    def step_2():
        if igdetective_home:  # Only run Step 2 if igdetective_home is provided
            try:
                logging.info(f"Step 2: Loci Location for {species}")
                loci_location(species, home, haploid, igdetective_home)
            except Exception as e:
                logging.error(f"Step 2 failed: {e}")
        else:
            logging.info(f"Skipping Step 2: No IGDetective home provided.")

    # Execute steps in parallel
    with ThreadPoolExecutor(max_workers=2) as executor:
        executor.submit(step_1)
        executor.submit(step_2)
    return run_flag

def process_steps_in_parallel(species, home, annotation):
    with ThreadPoolExecutor(max_workers=2) as executor:
        # Submit Step 5 and Step 6 as tasks to the thread pool
        try:
            # Submit Step 5 as a task to the thread pool
            future_coverage = executor.submit(coverage_analysis, species, home, annotation)
        except Exception as e:
            logging.error(f"Step 5 (coverage analysis) submission failed: {e}")
        
        try:
            # Submit Step 6 as a task to the thread pool
            future_cigar = executor.submit(cigar_processing, species, home, annotation)
        except Exception as e:
            logging.error(f"Step 6 (CIGAR processing) submission failed: {e}")
        
        try:
            # Wait for Step 5 to complete and handle any errors
            future_coverage.result()
        except Exception as e:
            logging.error(f"Step 5 (coverage analysis) execution failed: {e}")
        
        try:
            # Wait for Step 6 to complete and handle any errors
            future_cigar.result()
        except Exception as e:
            logging.error(f"Step 6 (CIGAR processing) execution failed: {e}")


def run_pipeline(args):
    """Run the pipeline with command-line parameters."""

    # Extract arguments
    species_list = args.species.split(",")
    home = args.home
    haploid = args.haploid
    fastqdir = args.fastqdir
    threads = args.t
    if args.igdetective_home:
        igdetective_home = args.igdetective_home
        customIG = None  # Set customIG to None because igdetective_home is provided
    else:
        igdetective_home = None  # Set igdetective_home to None because customIG is provided
        customIG = os.path.join(home, args.customIG) 

    configure_logging(home)

    # Validate input directories
    for dir_path, dir_name in [
        (home, "Home directory"),
        (f"{home}/{fastqdir}", "FASTQ directory"),
    ]:
        if not os.path.exists(dir_path):
            logging.error(f"{dir_name} not found: {dir_path}")
            raise FileNotFoundError(f"{dir_name} not found: {dir_path}")

    logging.info("Starting the pipeline...")

    for species in species_list:
        logging.info(f"Processing species: {species}")
        try:
            # Step 1 and 2: Data Preparation and Loci Location
            logging.info(f"Step 1 and Step 2: Running in parallel for {species}")
            run_flag = parallel_step_1_and_2(species, home, fastqdir, haploid, threads, igdetective_home)

            # Step 3: Convert Primary BAM
            output_primary_bam = os.path.join(home, "aligned_bam", species, f"{species}_merged_sorted_primary.bam")
            output_primary_index = f"{output_primary_bam}.csi"
            if not os.path.exists(output_primary_bam) or not os.path.exists(output_primary_index) or run_flag:
                logging.info(f"Step 3: Convert Primary BAM for {species}")
                convert_primary_bam(species, home, threads)
            else:
                logging.info(f"Step 3: Output exists, skipping Convert Primary BAM for {species} bc already completed")

            # Step 4: Final IG Loci
            if args.igdetective_home:
                logging.info(f"Step 4: Final IG Loci for {species}")
                final_ig_loci(species, home)

            # Step 5 + 6: Coverage Analysis and CIGAR Processing
            logging.info(f"Step 5 and Step 6: Running in parallel for {species}")
            if args.igdetective_home:
                annotation = os.path.join(home, "gene_position", f"{species}.final.Ig_loci.txt")
            else:
                annotation = os.path.join(customIG, f"{species}.customIG.txt")

            process_steps_in_parallel(species, home, annotation)

            logging.info(f"Species {species} processed successfully!")

        except Exception as e:
            logging.error(f"Error processing species {species}: {e}", exc_info=True)
            sys.exit(1)  # Exit with error code

    logging.info("Pipeline completed!")




if __name__ == "__main__":
    import argparse

    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Run the pipeline.")
    parser.add_argument("--species", required=True, help="Comma-separated list of species.")
    parser.add_argument("--home", required=True, help="Path to the home directory.")
    parser.add_argument("--haploid", required=True, help="Haploid status (True/False).")
    parser.add_argument("--fastqdir", required=True, help="Path to the FASTQ directory.")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--igdetective_home', type=str, help="Path to the IGDetective directory.")
    group.add_argument('--customIG', type=str, help="Path to directory containing ${species_name}.customIG.txt.")
    parser.add_argument("--t", required=False, type=int, default=32, help="# of threads to use (default: 32).")

    args = parser.parse_args()

    # Run the pipeline
    run_pipeline(args)
