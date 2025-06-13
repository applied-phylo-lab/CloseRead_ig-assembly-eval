import os
import logging
import argparse
import re
import sys
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
import pandas as pd
from fpdf import FPDF
import math
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
from PIL import Image
from intervaltree import Interval, IntervalTree
# Import necessary functions from CloseRead_plot.py
from .scripts.CloseRead_plot import plot_locus_length, plot_summary, plot_coverage, plot_mismatch_coverage
# Import necessary functions from CloseRead_pdf.py
from .scripts.CloseRead_pdf import make_pdf, make_pdfdi
# Import necessary functions from CloseRead_Stats.py
from .scripts.CloseRead_stats import read_species_from_file, create_directories, process_gene_data
from .scripts.CloseRead_stats import process_pileup, process_low_coverage_regions, process_read_file, coverage
from .scripts.CloseRead_stats import calculate_bin_counts, write_pileup, find_overlapping_mismatch_regions

def configure_logging(home):
    """Configure logging for plotting."""
    log_dir = os.path.join(home, "logs")
    os.makedirs(log_dir, exist_ok=True)

    logging.basicConfig(
        filename=os.path.join(log_dir, "plot.log"),
        filemode="w",
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
    )
    logging.info(f"Logging initialized. Logs will be written to {os.path.join(log_dir, 'plot.log')}")

def run_plot(args):
    """Run the plotting functionality."""
    home = args.home
    configure_logging(home)

    try:
        # Determine if single species or species file is provided
        if args.sf:
            species_list = read_species_from_file(args.sf)
        else:
            species_list = [args.s]  # Convert single species into a list
    
        gene = args.g
        haploid = False
        single_read_error = args.re
        readview_correct_threshold = args.rc
        chr1_color = "#6AABD7"
        chr2_color = "#F0DDB8"

        # Iterate through each species
        for species in species_list:
            logging.info(f"Generating plots for species: {species} - {gene}")

            haploid = False
            dirOut = f"{args.home}/errorPlots/{species}"
            dirStat = f"{args.home}/errorStats/{species}"
            create_directories(species, dirStat, dirOut)

            #process basepair-view mpileup file
            pri_pileup_file = f'{dirStat}/{gene}_pri_pileup.txt'
            alt_pileup_file = f'{dirStat}/{gene}_alt_pileup.txt'
            pri_pileup = process_pileup(pri_pileup_file)
            # Check if alternate pileup file exists
            if os.path.exists(alt_pileup_file):
                alt_pileup = process_pileup(alt_pileup_file)
                haploid = False
            else:
                print(f"Haploit = True : Alternate mpileup file not found: {alt_pileup_file}")
                alt_pileup = None
                haploid = True
                chr2 = None
                read_alt = None

            # Process read-oriented input
            read = process_read_file(f"{dirStat}/{gene}.txt", dirStat)
            chr1 = pri_pileup['Chrom'].value_counts().index[0]
            read_pri = read[(read['chromosome'] == chr1)]
            pri_pileup = pri_pileup[pri_pileup['Chrom']==chr1]
            if not haploid:
                # skipping alternate IG if more than 2 present
                if len(alt_pileup['Chrom'].value_counts()) > 2:
                    haploid = True
                    chr2 = None
                    read_alt = None
                else:
                    haploid = False
                    chr2 = alt_pileup['Chrom'].value_counts().index[0]
                    read_alt = read[(read['chromosome'] == chr2)]
                    alt_pileup = alt_pileup[alt_pileup['Chrom']==chr2]

            # Concatenate primary and alternate pileup data if diploid
            if alt_pileup is not None:
                merged_pileup = pd.concat([pri_pileup, alt_pileup])
            else:
                merged_pileup = pri_pileup
            merged_pileup.sort_values(by=['Chrom', 'Pos'], inplace=True)
            
            # Find and output break locations 
            break_regions_list, break_regions = process_low_coverage_regions(merged_pileup, lowCov_threshold=args.cov, padding=args.p)
            with open(f"{dirOut}/{gene}.break.txt", 'w') as file:
                for row in break_regions_list:
                    s = str(row)
                    file.write(s+'\n')
            start_break_pri= break_regions[break_regions['Chrom'] == chr1]['Start'].tolist()
            end_break_pri = break_regions[break_regions['Chrom'] == chr1]['End'].tolist()
            if not haploid:
                start_break_alt= break_regions[break_regions['Chrom'] == chr2]['Start'].tolist()
                end_break_alt = break_regions[break_regions['Chrom'] == chr2]['End'].tolist()

            
            #Compute read-view stats
            start_indices_pri, end_indices_pri, high_mismatch_bool_pri, positions_pri, coverage_counts_pri, zero_counts_pri, min_position_pri, max_position_pri, mid_counts_pri, high_mismatch_coverage_pri = coverage(read_pri, single_read_error, readview_correct_threshold)
            min_position_pri = min(pri_pileup['Pos'])
            max_position_pri = max(pri_pileup['Pos'])
            if not haploid:
                start_indices_alt, end_indices_alt, high_mismatch_bool_alt, positions_alt, coverage_counts_alt, zero_counts_alt, min_position_alt, max_position_alt, mid_counts_alt, high_mismatch_coverage_alt = coverage(read_alt, single_read_error, readview_correct_threshold)
                min_position_alt = min(alt_pileup['Pos'])
                max_position_alt = max(alt_pileup['Pos'])
    
            if not args.so:
                #Plot loci length
                plot_locus_length(pri_pileup, alt_pileup, gene, "#6AABD7", "#F0DDB8", dirOut, haploid, chr1, chr2)
                #Plot summary read info
                plot_summary(read_pri, read_alt, "#6AABD7", "#F0DDB8", haploid, dirOut, gene, chr1, chr2)
            
            #Write read-view mismatch to file
            with open(f"{dirOut}/{gene}.read.mismatch.txt", 'w') as file:
                for start, end in zip(start_indices_pri, end_indices_pri):
                    file.write(f"{chr1}:{start}-{end}\n")
                if not haploid:
                    for start, end in zip(start_indices_alt, end_indices_alt):
                        file.write(f"{chr2}:{start}-{end}\n")
    
            if not args.so:
                #Plot read coverage across loci
                plot_coverage(positions_pri, coverage_counts_pri, mid_counts_pri, zero_counts_pri, start_indices_pri,
                            end_indices_pri, high_mismatch_bool_pri.size, start_break_pri, end_break_pri, min_position_pri, max_position_pri, chr1, gene, dirOut)
                if not haploid:
                    plot_coverage(positions_alt, coverage_counts_alt, mid_counts_alt, zero_counts_alt, start_indices_alt, 
                                end_indices_alt, high_mismatch_bool_alt.size, start_break_alt, end_break_alt, min_position_alt, max_position_alt, chr2, gene, dirOut)
    
            #Compute basepair mismatch
            bin_count = calculate_bin_counts(pri_pileup, baseview_correct_threshold=80, bin_size=1000)
            if not haploid:
                alt_bin_count = calculate_bin_counts(alt_pileup, baseview_correct_threshold=80, bin_size=1000)
    
            #Write base-view mismatch exact/average rough position to file - removes the content if the files already exist
            with open(f"{dirOut}/{gene}.base.exactmismatch.csv", 'w') as f:
                f.write("Chrom,Pos,Correct,PercentCorrect,Depth,Indel\n")
            with open(f"{dirOut}/{gene}.base.avgmismatch.csv", 'w') as f:
                f.write("Chrom,Start,End,AvgPercentMismatch,AvgDepth,AvgIndel\n")
            grouped_baseMis_pri = write_pileup(pri_pileup, gene, dirOut)
            if not haploid:
                grouped_baseMis_alt = write_pileup(alt_pileup, gene, dirOut)
    
            #Compute overlapping mismatch regions between read-view and base-view
            overlaps_pri = find_overlapping_mismatch_regions(pri_pileup, start_indices_pri, end_indices_pri)
            with open(f"{dirOut}/{species}.{gene}.{chr1}.finalMismatch.csv", 'w') as file:
                pass
            overlaps_pri.to_csv(f"{dirOut}/{species}.{gene}.{chr1}.finalMismatch.csv", index=False)
            if not haploid:
                overlaps_alt = find_overlapping_mismatch_regions(alt_pileup, start_indices_alt, end_indices_alt)
                with open(f"{dirOut}/{species}.{gene}.{chr1}.finalMismatch.csv", 'w') as file:
                    pass
                overlaps_alt.to_csv(f"{dirOut}/{species}.{gene}.{chr2}.finalMismatch.csv", index=False)
    
            if not args.so:
                #Plot basepair mismatch across loci
                #Define the colors for the heatmap
                colors = [(1, 1, 1), (0.6, 0.6, 0.6), (0.2, 0.2, 0.2), (0, 0, 0)]  
                n_bins = 100  
                cmap_name = 'custom'
                cm = LinearSegmentedColormap.from_list(cmap_name, colors, N=n_bins)
                plot_mismatch_coverage(pri_pileup, bin_count, positions_pri, start_indices_pri, end_indices_pri, 
                                    high_mismatch_bool_pri.size, start_break_pri, end_break_pri, 
                                    chr1_color, chr1, gene, dirOut, cm, min_position_pri, max_position_pri)
                if not haploid:
                    plot_mismatch_coverage(alt_pileup, bin_count, positions_alt, start_indices_alt, end_indices_alt, 
                                    high_mismatch_bool_alt.size, start_break_alt, end_break_alt, 
                                    chr2_color, chr2, gene, dirOut, cm, min_position_alt, max_position_alt)
    
                #Generate PDF - without dotplot
                if args.m:
                    meta = pd.read_csv(args.m, sep=",")
                    LatinName = meta[meta['IndividualID'] == species]['LatinName'].item()
                    CommonName = meta[meta['IndividualID'] == species]['CommonName'].item()
                    Source = meta[meta['IndividualID'] == species]['Source'].item()
                    SourceLink = meta[meta['IndividualID'] == species]['SourceLink'].item()
                    Haplotype = meta[meta['IndividualID'] == species]['Haplotype Resolved'].item()
                    if Haplotype == "No":
                        hapkind = "Not Haplotype Resolved"
                    else:
                        hapkind = "Haplotype Resolved"
                        
                    output_filename = f"{dirOut}/{species}_{gene}_result.pdf"
                    if not haploid:
                        image_files = [
                            f'{dirOut}/{gene}.summary.allreads.png', 
                            f'{dirOut}/{gene}.length.png', 
                            f'{dirOut}/{gene}.{chr1}.readcoverage.all.png',
                            f'{dirOut}/{gene}.{chr1}.basecoverage.PerCorrect.png',
                            f'{dirOut}/{gene}.{chr2}.readcoverage.all.png',
                            f'{dirOut}/{gene}.{chr2}.basecoverage.PerCorrect.png',
                        ]
                        make_pdfdi(image_files, output_filename, species, CommonName, LatinName, hapkind, Source, overlapx=0.7, overlapy=0.01, scale_top=0.6, scale_bottom=0.3)
                    else:
                        image_files = [
                            f'{dirOut}/{gene}.summary.allreads.png', 
                            f'{dirOut}/{gene}.length.png', 
                            f'{dirOut}/{gene}.{chr1}.readcoverage.all.png',
                            f'{dirOut}/{gene}.{chr1}.basecoverage.PerCorrect.png',
                        ]
                        make_pdf(image_files, output_filename, species, CommonName, LatinName, hapkind, Source, overlapx=0.7, overlapy=0.01, scale_top=0.6, scale_bottom=0.3)
            if args.pg:
                genelevel_result = process_gene_data(gene_file=args.pg, merged_pileup=merged_pileup, read=read)
                outFile = f"{dirOut}/{species}.{gene}.genelevel.csv"
                genelevel_result.to_csv(outFile, sep='\t')
                print("gene level read support computed!")


    except Exception as e:
        logging.error(f"Error generating plots: {e}", exc_info=True)
        raise

    logging.info("Plot generation completed successfully!")


def resolve_path(path):
    return os.path.abspath(os.path.expanduser(path))


def run_plot_cli():
    """Command-line interface for the plotting script."""
    parser = argparse.ArgumentParser(description="Run CloseRead Final Evalutation and Plotting.")

    # Add mutually exclusive group for either species or species file input
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--s', type=str, metavar="species", help="Single species identifier (use this if providing one species)")
    group.add_argument('--sf', type=resolve_path, metavar="species_file", help="Path to file containing a list of species (use this if providing multiple species)")

    # Required arguments
    parser.add_argument('--g', type=str, required=True, metavar="gene", help="Gene identifier (IGH/IGK/IGL)")
    parser.add_argument('--home', type=resolve_path, required=True, metavar="home", help="Path to the home directory")

    # Optional arguments
    parser.add_argument('--cov', type=int, default=2, metavar="lowCov_threshold", help="Threshold for low coverage (default: 2)")
    parser.add_argument('--p', type=int, default=2000, metavar="padding", help="Padding around low coverage regions (default: 2000bps)")
    parser.add_argument('--re', type=float, default=0.01, metavar="single_read_error", help="Threshold for high mismatch consideration (default: 0.01)")
    parser.add_argument('--rc', type=int, default=5, metavar="readview_correct_threshold", help="Number of high mismatch reads needed to consider a position as high mismatch (default: 5)")
    parser.add_argument('--bc', type=int, default=5, metavar="baseview_correct_threshold", help="Threshold for exact match percent at a position for heatmap (default: 80 percent)")
    parser.add_argument('--m', type=resolve_path, metavar="meta", help="Path to the meta information .csv file, used for generating PDF.")
    parser.add_argument('--so', action='store_true', help="Output .txt and .csv files only, skip visualization.")
    parser.add_argument('--pg', type=resolve_path, metavar="gene_level assessment", help="Path for gene level annotation file, to generate gene-level read support information.")

    args = parser.parse_args()
    if not args.s and not args.sf:
        parser.error("You must provide either --s (single species) or --sf (species file).")

    run_plot(args)

