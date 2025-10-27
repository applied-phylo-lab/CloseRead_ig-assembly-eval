import subprocess
import os
from datetime import datetime
import logging
logger = logging.getLogger(__name__)

def cigar_processing(species, home, annotation):
    """Process CIGAR data."""

    # Define paths for script and input files
    current_dir = os.path.dirname(os.path.abspath(__file__))
    parent_dir = os.path.dirname(current_dir)
    script = os.path.join(parent_dir, "scripts/cigar.py")
    bam = os.path.join(home, "aligned_bam", species, f"{species}_merged_sorted_primary.bam")
    error_dir = os.path.join(home, "errorStats", species)
    os.makedirs(error_dir, exist_ok=True)

    try:
        logger.info(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} - Starting CIGAR processing for species: {species}")
        # Run the script
        res = subprocess.run(
            [
                "python",
                script,
                bam,
                annotation,
                species,
                error_dir,
            ],
            text=True,
            capture_output=True,
            check=True
        )
        if res.stdout:
            logger.info("[read-oriented analysis stdout]\n%s", res.stdout.strip())
        if res.stderr:
            logger.info("[read-oriented analysis stderr]\n%s", res.stderr.strip())
        logger.info(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} - read-oriented analysis completed successfully for species: {species}")

    except subprocess.CalledProcessError as e:
        logger.info(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} - Error during read-oriented analysis: {e}")
        if e.stdout:
            logger.error("[read-oriented analysis stdout]\n%s", e.stdout.strip())
        if e.stderr:
            logger.error("[read-oriented analysis stderr]\n%s", e.stderr.strip())
        raise RuntimeError(f"Failed read-oriented analysis for species {species}.") from e


if __name__ == "__main__":
    import argparse

    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Process CIGAR data.")
    parser.add_argument("--species", required=True, help="Species name.")
    parser.add_argument("--home", required=True, help="Path to the home directory.")
    parser.add_argument("--annotation", required=True, help="Path to the annotation file.")

    args = parser.parse_args()

    # Call the function
    cigar_processing(args.species, args.home, args.annotation)
