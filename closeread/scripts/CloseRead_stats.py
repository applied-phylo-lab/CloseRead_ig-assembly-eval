import re
import sys
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import os
import warnings
import pandas as pd
from intervaltree import Interval, IntervalTree
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
warnings.filterwarnings('ignore')

# Function to read species from a file
def read_species_from_file(file_path):
    with open(file_path, 'r') as f:
        species_list = [line.strip() for line in f.readlines()]
    return species_list

# Create output directories
def create_directories(species, dirStat, dirOut):
    # Create output directory if it doesn't exist
    if not os.path.exists(dirOut):
        os.makedirs(dirOut)
        print(f"Created directory: {dirOut}")
    else:
        print(f"Directory already exists: {dirOut}, will overwrite existing files")

    # Create stats directory if it doesn't exist
    if not os.path.exists(dirStat):
        os.makedirs(dirStat)
        print(f"Created directory: {dirStat}")



def parse_read_bases(read_bases):
    """
    Parse read bases to count correct matches (., and ,) while handling indels,
    read starts (^), read ends ($).
    """
    # Remove read start markers (^) along with the following character indicating mapping quality
    cleaned_bases = re.sub(r'\^.', '', read_bases)

    # Remove read end markers ($)
    cleaned_bases = cleaned_bases.replace('$', '')

    # handle indels: remove sequences following + or - indicating the length and actual indel
    cleaned_bases = re.sub(r'[\+\-](\d+)([ACGTNacgtn]+)', '', cleaned_bases)

    # handle pads and reverse pads
    cleaned_bases = cleaned_bases.replace('*', '').replace('#', '')

    # After removing all special cases, count '.' and ',' as correct matches
    correct = cleaned_bases.count('.') + cleaned_bases.count(',')

    return correct

def count_indels(read_bases):
    """
    Count the number of indels in the read bases, handling both insertions and deletions.
    Indels are indicated by + or - followed by the length of the indel and the sequence.
    """

    # Find all occurrences of indels, both insertions and deletions
    indels = re.findall(r'[\+\-](\d+)([ACGTNacgtn]+)', read_bases)

    # The total count of indels is the sum of the occurrences
    indel_count = sum(int(length) for length, _ in indels)

    return indel_count


def process_pileup(pileup_file):
    """
    Process a pileup file to extract relevant metrics (correct matches, indels, percent correct).
    
    Args:
    pileup_file (str): Path to the pileup file.
    
    Returns:
    Processed data as a pandas DataFrame.
    """
    results = {}  # To store results, keyed by (chromosome, position)
    with open(pileup_file, 'r') as f:
        for line in f:
            try:
                chrom, pos, ref_base, depth, read_bases, _ = line.split()[:6]
                correct = parse_read_bases(read_bases)
                indel = count_indels(read_bases)
                depth = int(depth)
                if depth > 0:
                    percent_correct = (correct / depth) * 100
                else:
                    percent_correct = 0
                # Store the result
                results[(chrom, pos)] = (correct, percent_correct, depth, indel)
            except Exception as e:
                print(f"Error processing line: {line.strip()}, Error: {e}")
                continue

    # Convert to a pandas DataFrame
    pileup_data = [(*key, *value) for key, value in results.items()]
    df_pileup = pd.DataFrame(pileup_data, columns=['Chrom', 'Pos', 'Correct', 'PercentCorrect', 'Depth', 'Indel'])
    
    # Convert necessary columns to numeric
    df_pileup['Pos'] = pd.to_numeric(df_pileup['Pos'])
    df_pileup['Correct'] = pd.to_numeric(df_pileup['Correct'])
    df_pileup['PercentCorrect'] = pd.to_numeric(df_pileup['PercentCorrect'])
    df_pileup['Depth'] = pd.to_numeric(df_pileup['Depth'])
    df_pileup['Indel'] = pd.to_numeric(df_pileup['Indel'])

    # Sort values by Chrom and Pos
    df_pileup.sort_values(by=['Chrom', 'Pos'], inplace=True)
    
    return df_pileup

def process_low_coverage_regions(pileup_df, lowCov_threshold=2, padding=2000):
    """
    Process low coverage regions from a pileup DataFrame.

    Args:
        pileup_df (pd.DataFrame): The pileup data containing 'Chrom', 'Pos', and 'Depth'.
        lowCov_threshold (int): Threshold for low coverage (default: 2).
        padding (int): Number of bases to pad the start and end of low coverage regions (default: 2000).
    
    Returns:
        list: List of low-coverage regions with padding applied, formatted as 'Chrom:Start-End'.
    """
    # Filter for low coverage positions
    low_coverage = pileup_df[pileup_df['Depth'] <= lowCov_threshold].copy()

    # Detect consecutive positions by checking where the difference in positions is greater than 2
    low_coverage['group'] = low_coverage.groupby('Chrom')['Pos'].transform(lambda x: (x.diff() > 2).cumsum())

    # Group by Chrom and group to find the start and end of continuous low-coverage regions
    break_regions = low_coverage.groupby(['Chrom', 'group']).agg({'Pos': ['min', 'max']}).reset_index()

    # Rename columns for clarity
    break_regions.columns = ['Chrom', 'group', 'Start', 'End']

    # Drop the 'group' column
    break_regions.drop('group', axis=1, inplace=True)

    # Create a list of break regions with padding
    break_regions_list = [
        f"{row['Chrom']}:{row['Start']-padding}-{row['End']+padding}"
        for index, row in break_regions.iterrows()
    ]

    return break_regions_list, break_regions



def process_read_file(file_path, dirStat):
    """
    Process the read-oriented file, converting necessary columns to numeric types and returning DataFrame
    
    Args:
    file_path (str): Path to the read file.
    dirStat (str): Directory containing the statistics files.
    
    Returns:
    list: Unique chromosomes found in the data.
    """
    # Define column names
    read_file_columns = [
        'read_name', 'chromosome', 'start', 'read_length', 'mapping_quality', 
        'mismatches', 'mismatch_rate', 'longindels', 'total_indel_length', 
        'indel_rate', 'soft_clipping', 'hard_clipping'
    ]
    
    # Read the file into a DataFrame
    read_file = pd.read_csv(file_path, sep="\t", names=read_file_columns)

    # List of numeric columns to convert
    numeric_columns = [
        'start', 'read_length', 'mapping_quality', 'mismatches', 'mismatch_rate', 
        'longindels', 'total_indel_length', 'indel_rate', 'soft_clipping', 'hard_clipping'
    ]
    
    # Convert all numeric columns at once
    read_file[numeric_columns] = read_file[numeric_columns].apply(pd.to_numeric)
    read_file['end'] = read_file['start'] + read_file['read_length'] - 1
    
    return read_file

def coverage(read_data, single_read_error, readview_correct_threshold):
    """
    Calculate coverage and high mismatch regions for a single haplotype.

    Args:
    read_data (pd.DataFrame): Read data containing 'start', 'read_length', 'mapping_quality', and 'mismatch_rate'.
    single_read_error (float): Mismatch rate threshold. Defaults to 0.01.
    readview_correct_threshold (int): Threshold for high mismatch counts. Defaults to 5.

    Returns:
    Various arrays containing coverage and high mismatch data for the haplotype.
    """
    
    min_position = read_data['start'].min()
    max_position = (read_data['start'] + read_data['read_length']).max()

    """
    Initialize lists for coverage tracking
    `coverage_counts` stores the read coverage with mapping quality 60 for each position between min/max_position of the loci
    `zero_counts` stores the read coverage with mapping quality 0 for each position between min/max_position of the loci
    `mid_counts` stores the read coverage with mapping quality 1-59 for each position between min/max_position of the loci
    `high_mismatch_coverage` stores the number of high mismatch rate reads covering each position between min/max_position of the loci
    """
    coverage_counts = np.zeros(max_position - min_position + 1, dtype=int)
    zero_counts = np.zeros(max_position - min_position + 1, dtype=int)
    mid_counts = np.zeros(max_position - min_position + 1, dtype=int)
    high_mismatch_coverage = np.zeros(max_position - min_position + 1, dtype=int)

    #iterate all reads
    for _, row in read_data.iterrows():
        start_index = row['start'] - min_position
        end_index = start_index + row['read_length']
        # for all position covered by this read, coverage += 1
        if row['mapping_quality'] == 60:
            coverage_counts[start_index:end_index] += 1
        elif row['mapping_quality'] == 0:
            zero_counts[start_index:end_index] += 1
        else:
            mid_counts[start_index:end_index] += 1
        # if this read's mismatch rate exceed single_read_error(default 0.01), all position's high_mismatch_coverage += 1
        if row['mismatch_rate'] > single_read_error:
            high_mismatch_coverage[start_index:end_index] += 1

    # Get the positions with high mismatch rate from read-view
    positions = np.arange(min_position, max_position + 1)
    high_mismatch_positions = positions[high_mismatch_coverage > readview_correct_threshold]

    high_mismatch_bool = high_mismatch_coverage > readview_correct_threshold
    diff = np.diff(high_mismatch_bool.astype(int))
    start_indices = np.where(diff == 1)[0] + 1 + min_position
    end_indices = np.where(diff == -1)[0] + 1 + min_position
    
    # Handle edge cases
    # If the first position is a high mismatch, add 0 at the beginning
    if high_mismatch_bool[0]:
        start_indices = np.insert(start_indices, 0, 0)
    # If the last position is a high mismatch, add the last index at the end
    if high_mismatch_bool[-1]:
        end_indices = np.append(end_indices, high_mismatch_bool.size)

    return start_indices, end_indices, high_mismatch_bool, positions, coverage_counts, zero_counts, min_position, max_position, mid_counts, high_mismatch_coverage



def calculate_bin_counts(pileup, baseview_correct_threshold=80, bin_size=1000):
    """
    Filters the `pileup` data based on a percent correct threshold, bins the positions, 
    and calculates the percentage of well supported reads for each bin.

    Args:
    pileup (pd.DataFrame): The input DataFrame containing the pileup data.
    baseview_correct_threshold (int, optional): The threshold for filtering based on PercentCorrect. Defaults to 80.
    bin_size (int, optional): The size of each bin (in genomic positions). Defaults to 1000.

    Returns:
    pd.DataFrame: DataFrame containing bin start positions, counts, and percentages of well supported bases in bin_size.
    """
    # Filter the data based on the percent correct threshold
    pileup_filtered = pileup[pileup['PercentCorrect'] > baseview_correct_threshold]

    # Create bin start positions
    pileup_filtered['Bin_Start_Pos'] = (pileup_filtered['Pos'] // bin_size) * bin_size

    # Group by bin start position and count occurrences
    bin_count = pileup_filtered.groupby('Bin_Start_Pos').size()

    # Convert the bin counts into a DataFrame and calculate the well supported bases percentage
    bin_count = pd.DataFrame(bin_count, columns=["wellCount"])
    bin_count["wellCount percent"] = bin_count["wellCount"] / bin_size * 100

    # Reset index for cleaner DataFrame
    bin_count = bin_count.reset_index()

    return bin_count

def write_pileup(pileup, gene, dirOut):
    # Write base exact mismatch to CSV
    pileup.to_csv(f"{dirOut}/{gene}.base.exactmismatch.csv", mode='a', header=False, index=False)

    # Filter rows where PercentCorrect is less than 80
    baseMis = pileup[pileup['PercentCorrect'] < 80]
    
    # Calculate PercentMismatch
    baseMis["PercentMismatch"] = 100 - baseMis['PercentCorrect']

    # Create regions for positions where the difference is >= 1000
    baseMis['PosDiff'] = baseMis['Pos'].diff().fillna(0)
    baseMis['Region'] = (baseMis['PosDiff'] >= 1000).cumsum()

    # Group by regions and calculate averages
    grouped_baseMis = baseMis.groupby('Region').agg(
        Chrom=('Chrom', 'first'),
        Start=('Pos', 'min'),
        End=('Pos', 'max'),
        AvgPercentMismatch=('PercentMismatch', 'mean'),
        AvgDepth=('Depth', 'mean'),
        AvgIndel=('Indel', 'mean')
    ).reset_index(drop=True)

    # Write average mismatch regions to CSV
    grouped_baseMis.to_csv(f"{dirOut}/{gene}.base.avgmismatch.csv", mode='a', header=False, index=False)
    return grouped_baseMis

def find_overlapping_mismatch_regions(pileup, start_indices, end_indices, percent_threshold=80, pos_diff_threshold=1000):
    """
    Function to find overlapping regions between pileup intervals and regions with high mismatch.

    Args:
    pileup (pd.DataFrame): DataFrame with pileup data containing 'Chrom', 'Pos', 'PercentCorrect', 'Depth', and 'Indel' columns.
    start_indices (list): List of start indices from read-view high mismatch.
    end_indices (list): List of end indices from read-view high mismatch.
    percent_threshold (int): Threshold for PercentCorrect to define mismatch regions (default: 80).
    pos_diff_threshold (int): Position difference threshold to define non-consecutive regions (default: 1000).
    
    Returns:
    pd.DataFrame: DataFrame containing overlapping high mismatch regions with chromosome name, start, and end positions.
    """
    
    # Filter the pileup to get the regions with PercentCorrect below the threshold (high mismatch)
    baseMis = pileup[pileup['PercentCorrect'] < percent_threshold]
    baseMis["PercentMismatch"] = 100 - baseMis['PercentCorrect']

    # Group consecutive positions (within pos_diff_threshold) into regions
    baseMis['PosDiff'] = baseMis['Pos'].diff().fillna(0)
    baseMis['Region'] = (baseMis['PosDiff'] >= pos_diff_threshold).cumsum()

    # Group the base mismatch regions and calculate average stats
    grouped_baseMis = baseMis.groupby('Region').agg(
        Chrom=('Chrom', 'first'),
        Start=('Pos', 'min'),
        End=('Pos', 'max'),
        AvgPercentMismatch=('PercentMismatch', 'mean'),
        AvgDepth=('Depth', 'mean'),
        AvgIndel=('Indel', 'mean')
    ).reset_index(drop=True)

    # Build an IntervalTree from the grouped_baseMis regions
    baseMis_tree = IntervalTree()
    for idx, row in grouped_baseMis.iterrows():
        baseMis_tree.add(Interval(row['Start'], row['End'] + 1))  # Add 1 for half-open interval

    # Find overlaps between read-view intervals and base-view intervals
    overlapping_intervals = []
    for pileup_start, pileup_end in zip(start_indices, end_indices):
        overlapping = baseMis_tree.overlap(pileup_start, pileup_end + 1)
        
        # For each overlap, compute the actual overlapping region
        for interval in overlapping:
            overlap_start = max(pileup_start, interval.begin)  # Find the max of the start positions
            overlap_end = min(pileup_end, interval.end - 1)    # Find the min of the end positions
            overlapping_intervals.append({
                'chr': grouped_baseMis['Chrom'][0],  # Get the chromosome for the overlap
                'overlap_start': overlap_start,
                'overlap_end': overlap_end
            })

    overlaps_df = pd.DataFrame(overlapping_intervals)

    return overlaps_df


def process_gene_data(gene_file, merged_pileup, read):
    """
    Processes the gene file and performs calculations such as coverage, mismatch rates, 
    and fully spanning reads for each gene, handling two types of input formats.

    Parameters:
    gene_file (str): The path to the gene file (CSV format).
    merged_pileup (pd.DataFrame): DataFrame containing pileup data.
    read (pd.DataFrame): DataFrame containing read data.

    Returns:
    pd.DataFrame: Updated DataFrame with additional columns for each gene.
    """

    # Load the gene file
    genes = pd.read_csv(gene_file, sep='\t')

    # Check if the input file is precomputed format (Gene, Chromosome, Strand, Start, End)
    if 'Start' in genes.columns and 'End' in genes.columns:
        # Detected Type 2 format
        print("Detected precomputed format: Gene, Chromosome, Strand, Start, End")
        
        # Rename columns to match the target format
        genes.rename(columns={'Gene': 'Gene', 'Chromosome': 'Contig', 'Start': 'Pos', 'End': 'EndPos'}, inplace=True)
        # Calculate sequence length from start and end positions
        genes['SeqLength'] = genes['EndPos'] - genes['Pos'] + 1

    # Check if the input file is IgDetected format (GeneType, Contig, Pos, Sequence, etc.)
    elif 'Pos' in genes.columns and 'Sequence' in genes.columns:
        # Detected Type 1 format
        print("Detected IgDetected format: GeneType Contig Pos Strand Sequence Productive Locus")
        
        # Calculate the sequence length and end position
        genes['SeqLength'] = genes['Sequence'].apply(len)
        genes['EndPos'] = genes['Pos'] + genes['SeqLength'] - 1
        # Create Gene column by concatenating 'Locus' and 'GeneType'
        genes['Gene'] = genes['Locus'] + '_' + genes['GeneType']

    else:
        raise ValueError("Unsupported file format. Ensure the input file is in one of the supported formats. \n 1. Gene, Chromosome, Strand, Start, End \n 2. GeneType Contig Pos Strand Sequence Productive Locus")

    # Iterate over each gene in the DataFrame
    for index, row in genes.iterrows():
        # Extract gene information from the row
        chrom = row['Contig']  # Chromosome/contig name
        start = row['Pos']  # Starting position of the gene
        length = row['SeqLength']  # Length of the gene
        end = row['EndPos']  # Ending position of the gene
        gene_name = row['Gene']  # Gene name (constructed for Type 1, direct for Type 2)

        # Filter the pileup data for the exact chromosome and gene region
        gene_df = merged_pileup[(merged_pileup['Chrom'] == chrom) & (merged_pileup['Pos'] >= start) & (merged_pileup['Pos'] <= end)]

        # Calculate positions with at least 10x coverage
        positions_with_10x = (gene_df['Depth'] >= 10).sum()

        # Calculate mismatch rate and classify positions
        gene_df['mismatch_rate'] = 100 - gene_df['PercentCorrect']
        mismatch_positions = (gene_df['mismatch_rate'] > 20).sum()
        match_positions = (gene_df['mismatch_rate'] <= 20).sum()

        # Filter reads that fully span the gene region on the correct chromosome
        reads_spanning_region = read[(read['chromosome'] == chrom) & (read['start'] <= start) & (read['end'] >= end)].shape[0]
        fully_spanning_reads_100 = read[(read['chromosome'] == chrom) & (read['start'] <= start) & (read['end'] >= end) & (read['mismatches'] == 0)].shape[0]

        # Calculate average coverage and percent accuracy
        if length > 0:
            average_coverage = reads_spanning_region / length
            percent_accuracy = (match_positions / length) * 100
        else:
            average_coverage = 0
            percent_accuracy = 0

        # Compile mismatch and match positions as strings for reporting
        position_mismatches = ';'.join(map(str, gene_df['Depth'] - gene_df['Correct']))
        position_matches = ';'.join(map(str, gene_df['Correct']))

        # Update the gene DataFrame with calculated metrics
        genes.loc[index, 'Average_Coverage'] = average_coverage
        genes.loc[index, 'Mismatched_Positions'] = mismatch_positions
        genes.loc[index, 'Matched_Positions'] = match_positions
        genes.loc[index, 'Position_Mismatches'] = position_mismatches
        genes.loc[index, 'Position_Matches'] = position_matches
        genes.loc[index, 'Percent_Accuracy'] = percent_accuracy
        genes.loc[index, 'Positions_With_At_Least_10x_Coverage'] = positions_with_10x
        genes.loc[index, 'Fully_Spanning_Reads'] = reads_spanning_region
        genes.loc[index, 'Fully_Spanning_Reads_100%_Match'] = fully_spanning_reads_100

    # Return the updated DataFrame
    return genes

