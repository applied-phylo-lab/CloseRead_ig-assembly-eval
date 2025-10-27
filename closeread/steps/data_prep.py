import os
import subprocess
from datetime import datetime
import logging

logger = logging.getLogger(__name__)

def data_prep(species, home, fastqdir, haploid, threads):
    """Run data preparation."""

    # Define the script path
    current_dir = os.path.dirname(os.path.abspath(__file__))
    parent_dir = os.path.dirname(current_dir)
    script = os.path.join(parent_dir, "scripts/dataPrepAutomated.sh")
    logger.info(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} - Running data preparation...")
    # Build the command
    cmd = [
        script,
        "-s", species,
        "-w", fastqdir,
        "-h", haploid,
        "-d", home,
        "-t", threads
    ]

    # Run the command, redirecting stdout and stderr to the log file
    try:
        res = subprocess.run(cmd,
                            text=True,
                            capture_output=True,
                            check=True)
        if res.stdout:
            logger.info("[minimap stdout]\n%s", res.stdout.strip())
        if res.stderr:
            logger.info("[minimap stderr]\n%s", res.stderr.strip())
        logger.info("Data preparation completed successfully.")
    except subprocess.CalledProcessError as e:
        # Log at ERROR so it actually goes to pipeline.log
        logger.error("Data preparation failed (exit code %s)", e.returncode)
        if e.stdout:
            logger.error("[minimap stdout]\n%s", e.stdout.strip())
        if e.stderr:
            logger.error("[minimap stderr]\n%s", e.stderr.strip())
        raise


if __name__ == "__main__":
    import argparse

    # Command-line argument parsing
    parser = argparse.ArgumentParser(description="Run the data preparation step.")
    parser.add_argument("--species", required=True, help="Species name.")
    parser.add_argument("--home", required=True, help="Path to the home directory.")
    parser.add_argument("--fastqdir", required=True, help="Path to the FASTQ directory.")
    parser.add_argument("--haploid", required=True, help="Haploid status (True or False).")
    parser.add_argument("--t", required=False, type=int, default=32, help="# of threads to use (default: 32).")

    args = parser.parse_args()

    # Run the data preparation function with parsed arguments
    data_prep(args.species, args.home, args.fastqdir, args.haploid)
