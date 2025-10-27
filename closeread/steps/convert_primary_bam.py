import subprocess
import os
from datetime import datetime
import logging
logger = logging.getLogger(__name__)

def convert_primary_bam(species, home, threads):
    """Convert merged BAM to primary BAM."""

    # Define input and output file paths
    bam = os.path.join(home, "aligned_bam", species, f"{species}_merged_sorted.bam")
    output_bam = os.path.join(home, "aligned_bam", species, f"{species}_merged_sorted_primary.bam")
    output_index = f"{output_bam}.csi"

    # Skip step if output already exists
    if os.path.exists(output_bam) and os.path.exists(output_index):
        logger.info(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} - Primary BAM and index already exist: {output_bam}, {output_index}")
        return

    # Run samtools commands and log output
    try:
        logger.info(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} - Converting BAM to primary BAM for species: {species}")
        logger.info(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} - Input BAM: {bam}")
        logger.info(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} - Output BAM: {output_bam}")
        # Run samtools view
        res = subprocess.run(
            ["samtools", "view", "-b", "-F", "0x800", "-F", "0x100", "-@", str(threads), bam, "-o", output_bam],
            text=True,
            capture_output=True,
            check=True
        )
        if res.stdout:
            logger.info("[samtools view stdout]\n%s", res.stdout.strip())
        if res.stderr:
            logger.info("[samtools view stderr]\n%s", res.stderr.strip())
        logger.info(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} - samtools view completed successfully.")
        # Run samtools index
        res_index = subprocess.run(
            ["samtools", "index", "-c", output_bam],
            text=True,
            capture_output=True,
            check=True
        )
        if res_index.stdout:
            logger.info("[samtools index stdout]\n%s", res_index.stdout.strip())
        if res_index.stderr:
            logger.info("[samtools index stderr]\n%s", res_index.stderr.strip())
        logger.info(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} - samtools index completed successfully.")

    except subprocess.CalledProcessError as e:
        logger.info(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} - Error during BAM conversion: {e}")
        if e.stdout:
            logger.error("[samtools stdout]\n%s", e.stdout.strip())
        if e.stderr:
            logger.error("[samtools stderr]\n%s", e.stderr.strip())
        raise RuntimeError(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} - Failed to convert BAM for species {species}.") from e


if __name__ == "__main__":
    import argparse

    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Convert merged BAM to primary BAM.")
    parser.add_argument("--species", required=True, help="Species name.")
    parser.add_argument("--home", required=True, help="Path to the home directory.")
    parser.add_argument("--t", required=False, type=int, default=32, help="# of threads to use (default: 32).")

    args = parser.parse_args()

    # Call the function
    convert_primary_bam(args.species, args.home)
