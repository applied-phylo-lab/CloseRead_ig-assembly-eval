import subprocess
import os
from datetime import datetime
import logging
logger = logging.getLogger(__name__)

def final_ig_loci(species, home):
    """Process loci into final IG loci."""
    # Define the log file path

    output_dir = os.path.join(home, "gene_position")
    os.makedirs(output_dir, exist_ok=True)

    # Define script and output paths
    current_dir = os.path.dirname(os.path.abspath(__file__))
    parent_dir = os.path.dirname(current_dir)
    script = os.path.join(parent_dir, "scripts/finalGene.py")
    output = os.path.join(home, "gene_position", f"{species}.final.Ig_loci.txt")

    try:
        # Log start of the process
        logger.info(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} - Starting final IG loci processing for species: {species}")
        logger.info(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} - Output file: {output}")

        subprocess.run(
            ["python", script, species, home],
            text=True,
            capture_output=True,
            check=True
        )
        logger.info(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} - Final IG loci generated for {species}.")


    except subprocess.CalledProcessError as e:
        logger.info(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} - Error during IG loci processing: {e}")
        if e.stdout:
            logger.error("[IG loci position compute stdout]\n%s", e.stdout.strip())
        if e.stderr:
            logger.error("[IG loci position compute stderr]\n%s", e.stderr.strip())
        raise RuntimeError(f"Failed to process final IG loci for species {species}.") from e


if __name__ == "__main__":
    import argparse

    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Process loci into final IG loci.")
    parser.add_argument("--species", required=True, help="Species name.")
    parser.add_argument("--home", required=True, help="Path to the home directory.")

    args = parser.parse_args()

    # Call the function
    final_ig_loci(args.species, args.home)