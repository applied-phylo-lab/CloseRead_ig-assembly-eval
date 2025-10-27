"""
CloseRead: A Python package for IG assembly evaluation.
"""

from .pipeline import run_pipeline
from .steps.data_prep import data_prep
from .steps.convert_primary_bam import convert_primary_bam
from .steps.loci_location import loci_location
from .steps.final_ig_loci import final_ig_loci
from .steps.cigar_processing import cigar_processing
from .steps.coverage_analysis import coverage_analysis

__all__ = [
    "run_pipeline",
    "data_prep",
    "convert_primary_bam",
    "loci_location",
    "final_ig_loci",
    "cigar_processing",
    "coverage_analysis",
]


__version__ = "1.2.1"
