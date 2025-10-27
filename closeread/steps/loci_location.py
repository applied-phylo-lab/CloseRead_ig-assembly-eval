import subprocess
import os
from datetime import datetime
import logging
logger = logging.getLogger(__name__)

def loci_location(species, home, haploid, igdetective_home):
    """Run loci location detection for primary and alternate genomes."""

    # Define input and output paths
    pri_genome = os.path.join(home, "assemblies", f"{species}.pri.fasta")
    alt_genome = os.path.join(home, "assemblies", f"{species}.alt.fasta")
    pri_outdir = os.path.join(home, "igGene", f"{species}.pri.igdetective")
    alt_outdir = os.path.join(home, "igGene", f"{species}.alt.igdetective")

    try:
        # Run for primary genome
        os.makedirs(pri_outdir, exist_ok=True)
        logger.info(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} - Running IgDetective for primary assembly: {pri_genome}")
        res = subprocess.run(
            [
                "python",
                os.path.join(igdetective_home, "run_iterative_igdetective.py"),
                pri_genome,
                pri_outdir,
            ],
            text=True,
            capture_output=True,
            check=True
        )
        if res.stdout:
            logger.info("[IgDetective stdout]\n%s", res.stdout.strip())
        if res.stderr:
            logger.info("[IgDetective stderr]\n%s", res.stderr.strip())
        logger.info(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} - Primary assembly IgDetective completed successfully. Output: {pri_outdir}")

        # Run for alternate genome if haploid is False
        if haploid == "False":
            os.makedirs(alt_outdir, exist_ok=True)
            logger.info(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} - Running IgDetective for alternate assembly: {alt_genome}")
            alt_res = subprocess.run(
                [
                    "python",
                    os.path.join(igdetective_home, "run_iterative_igdetective.py"),
                    alt_genome,
                    alt_outdir,
                ],
                text=True,
                capture_output=True,
                check=True
            )
            if alt_res.stdout:
                logger.info("[IgDetective stdout]\n%s", alt_res.stdout.strip())
            if alt_res.stderr:
                logger.info("[IgDetective stderr]\n%s", alt_res.stderr.strip())
            logger.info(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} - Alternate assembly IgDetective completed successfully. Output: {alt_outdir}")

    except subprocess.CalledProcessError as e:
        logger.info(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} - IgDetective failed: {e}")
        if e.stdout:
            logger.error("[IgDetective stdout]\n%s", e.stdout.strip())
        if e.stderr:
            logger.error("[IgDetective stderr]\n%s", e.stderr.strip())
        raise RuntimeError(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} - IgDetective failed for species {species}.") from e



if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Run loci location detection.")
    parser.add_argument("--species", required=True, help="Species name.")
    parser.add_argument("--home", required=True, help="Path to the home directory.")
    parser.add_argument("--haploid", required=True, help="Haploid status (True or False).")
    parser.add_argument("--igdetective_home", required=True, help="Path to the igDetective home directory.")

    args = parser.parse_args()

    loci_location(args.species, args.home, args.haploid, args.igdetective_home)
