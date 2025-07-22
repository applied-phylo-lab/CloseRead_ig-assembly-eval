import subprocess
import os
from datetime import datetime
import logging
logger = logging.getLogger(__name__)

def coverage_analysis(species, home, annotation):
    """Run coverage analysis."""

    # Define input paths and directories
    current_dir = os.path.dirname(os.path.abspath(__file__))
    parent_dir = os.path.dirname(current_dir)
    script = os.path.join(parent_dir, "scripts/coverage.sh")
    assembly = os.path.join(home, "assemblies", f"{species}.merged.fasta")
    bam = os.path.join(home, "aligned_bam", species, f"{species}_merged_sorted_primary.bam")
    error_dir = os.path.join(home, "errorStats", species)
    os.makedirs(error_dir, exist_ok=True)

    try:
        # Log the start of the process
        logger.info(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} - Starting base-oriented analysis for species: {species}")
        logger.info(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} - Annotation file: {annotation}")
        logger.info(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} - errorStats directory: {error_dir}")
        # Run the script
        res = subprocess.run(
            f"{script} -s {species} -a {assembly} -b {bam} -f {annotation} -d {home}",
            text=True,
            capture_output=True,
            shell=True,
            check=True
        )
        if res.stdout:
            logger.info("[base-oriented analysis stdout]\n%s", res.stdout.strip())
        if res.stderr:
            logger.info("[base-oriented analysis stderr]\n%s", res.stderr.strip())
        logger.info(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} - base-oriented analysis completed successfully for species: {species}")
    except subprocess.CalledProcessError as e:
        logger.info(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} - Error during base-oriented analysis: {e}")
        if e.stdout:
            logger.error("[base-oriented analysis stdout]\n%s", e.stdout.strip())
        if e.stderr:
            logger.error("[base-oriented analysisstderr]\n%s", e.stderr.strip())
        raise RuntimeError(f"Failed to complete base-oriented analysis for species {species}.") from e


if __name__ == "__main__":
    import argparse

    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Run coverage analysis.")
    parser.add_argument("--species", required=True, help="Species name.")
    parser.add_argument("--home", required=True, help="Path to the home directory.")
    parser.add_argument("--annotation", required=True, help="Path to the annotation file.")

    args = parser.parse_args()

    # Call the function
    coverage_analysis(args.species, args.home, args.annotation)
