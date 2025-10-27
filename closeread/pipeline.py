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
from . import logging_config
from .logging_config import setup_global_logging

logger = logging.getLogger(__name__)

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


def parallel_step_1_and_2(species, home, fastqdir, haploid, threads, igdetective_home=None):
    output_bam = os.path.join(home, "aligned_bam", species, f"{species}_merged_sorted.bam")
    output_index = f"{output_bam}.csi"
    igdetective_output = os.path.join(home, "igGene", f"{species}.pri.igdetective", "combined_genes_IGL.txt")
    run_flag = False
    if not os.path.exists(output_bam) or not os.path.exists(output_index):
        run_flag = True
    else:
        run_flag = False

    def step_1():
        if not os.path.exists(output_bam) or not os.path.exists(output_index):
            try:
                logger.info(f"Step 1: Data Preparation for {species}")
                data_prep(species, home, fastqdir, haploid, str(threads))
            except Exception as e:
                logger.error(f"Step 1 failed: {e}")
                sys.exit(1)
        else:
            logger.info(f"Skipping Step 1: Output exists for {species}. If rerun is needed, delete {output_bam} and {output_index}.")        

    def step_2():
        if igdetective_home and not os.path.exists(igdetective_output):  # Only run Step 2 if igdetective_home is provided and output files incomplete
            try:
                logger.info(f"Step 2: IgDetective for {species}")
                loci_location(species, home, haploid, igdetective_home)
            except Exception as e:
                logger.error(f"Step 2 failed: {e}")
                sys.exit(1)
        else:
            logger.info(f"Skipping Step 2 IgDetective: output exist or custom IG file provided. If rerun is needed, delete {home}/igGene/{species}.pri.igdetective")

    # Execute steps in parallel
    with ThreadPoolExecutor(max_workers=2) as executor:
        executor.submit(step_1)
        executor.submit(step_2)
    return run_flag

def process_steps_in_parallel(species, home, annotation):
    with ThreadPoolExecutor(max_workers=2) as executor:
        # Submit Step 5 and Step 6 as tasks to the thread pool
        logger.info(f"Step 5: (base-oriented analysis) for {species}")
        future_coverage = executor.submit(coverage_analysis, species, home, annotation)

        # Submit Step 6 as a task to the thread pool
        logger.info(f"Step 6: (read-oriented analysis) for {species}")
        future_cigar = executor.submit(cigar_processing, species, home, annotation)

        try:
            future_coverage.result()
            logger.info(f"Step 5 (base-oriented analysis) completed successfully for {species}.")
        except Exception as e:
            logger.error(f"Step 5 (base-oriented analysis) failed: {e}")
            sys.exit(1)
        
        try:
            future_cigar.result()
            logger.info(f"Step 6 (read-oriented analysis) completed successfully for {species}.")
        except Exception as e:
            logger.error(f"Step 6 (read-oriented analysis) failed: {e}")
            sys.exit(1)


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

    log_dir = os.path.join(home, 'closeread_logs')
    setup_global_logging(log_dir, level="INFO", console=True)

    # Validate input directories
    for dir_path, dir_name in [
        (home, "Home directory"),
        (f"{home}/{fastqdir}", "FASTQ directory"),
    ]:
        if not os.path.exists(dir_path):
            print(f"{dir_name} not found: {dir_path}")
            raise FileNotFoundError(f"{dir_name} not found: {dir_path}")

    for species in species_list:
        logger.info(f'=== Starting pipeline for {species} ===')
        try:
            # Step 1 and 2: Data Preparation and Loci Location
            run_flag = parallel_step_1_and_2(species, home, fastqdir, haploid, threads, igdetective_home)

            # Step 3: Convert Primary BAM
            output_primary_bam = os.path.join(home, "aligned_bam", species, f"{species}_merged_sorted_primary.bam")
            output_primary_index = f"{output_primary_bam}.csi"
            if not os.path.exists(output_primary_bam) or not os.path.exists(output_primary_index) or run_flag:
                logger.info(f"Step 3: Select only primary alignments for {species}")
                convert_primary_bam(species, home, threads)
                logger.info(f"Step 3 [Select Primary alignments] completed successfully for {species}.")
            else:
                logger.info(f"Skipping Step 3: Output exists for {species}. If rerun is needed, delete files in aligned_bam/{species} folder.")

            # Step 4: Final IG Loci
            if args.igdetective_home:
                logger.info(f"Step 4: Compute IG Loci position file for {species}")
                final_ig_loci(species, home)
                logger.info(f"Step 4 [Compute IG Loci] completed successfully for {species}.")
            else:
                logger.info(f"Skipping Step 4: Custom IG file provided for {species}")

            # Step 5 + 6: Coverage Analysis and CIGAR Processing
            if args.igdetective_home:
                logger.info(f"Using IG loci found by Igdetetive for {species}")
                annotation = os.path.join(home, "gene_position", f"{species}.final.Ig_loci.txt")
            else:
                logger.info(f"Using provided custom IG loci for {species}")
                annotation = os.path.join(customIG, f"{species}.customIG.txt")

            process_steps_in_parallel(species, home, annotation)

            logger.info(f"{species} processed successfully!")

        except Exception as e:
            logger.error(f"Error processing species {species}: {e}", exc_info=True)
            sys.exit(1)  # Exit with error code

    logger.info("Pipeline completed!")




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
