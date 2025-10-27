import re
import sys
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import os
import warnings
import pandas as pd
from fpdf import FPDF
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
warnings.filterwarnings('ignore')

#Visualization Functions


def plot_locus_length(pri_pileup, alt_pileup, gene, chr1_color, chr2_color, dirOut, haploid, chr1, chr2=None):
    """
    Plot the locus length for the primary and alternate haplotype.
    
    Args:
    pri_pileup (pd.DataFrame): Primary pileup data.
    alt_pileup (pd.DataFrame): Alternate pileup data (None if haploid).
    gene (str): Gene name for labeling the plot.
    chr1_color (str): Color for chromosome 1.
    chr2_color (str): Color for chromosome 2 (only used if not haploid).
    dirOut (str): Output directory to save the plot.
    haploid (bool): True if haploid, False otherwise.
    chr1 (str): Chromosome 1.
    chr2 (str, optional): Chromosome 2, only used if not haploid.
    """
    # Plot setup
    if not haploid:
        lengthPlt, ax = plt.subplots(figsize=(1.5, 4))
        chrom = [chr1, chr2]
        length = [pri_pileup['Pos'].max() - pri_pileup['Pos'].min(), alt_pileup['Pos'].max() - alt_pileup['Pos'].min()]
        colors = [chr1_color, chr2_color]
    else:
        lengthPlt, ax = plt.subplots(figsize=(1, 4))
        chrom = [chr1]
        length = [pri_pileup['Pos'].max() - pri_pileup['Pos'].min()]
        colors = [chr1_color]

    # Create the bar plot
    y_pos = np.arange(len(chrom))
    ax.bar(y_pos, length, align='center', color=colors, width=0.9 if not haploid else 0.1)
    ax.set_xticks(y_pos)
    ax.set_xticklabels(chrom)
    ax.set_xlabel(f"{gene} locus\nlength", weight='bold', size=10)

    # Customize the plot
    for spine in ['right', 'top', 'left', 'bottom']:
        ax.spines[spine].set_visible(False)
    
    ax.tick_params(axis="both", which="both", bottom=False, top=False, labelbottom=True, left=True, right=False, labelleft=True)
    ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))

    # Add dashed lines for y-axis ticks
    for tick in ax.get_yticks():
        ax.axhline(y=tick, linestyle='dashed', alpha=0.4, color='#eeeeee', zorder=1)
    
    # Rotate x-axis labels
    for tick in ax.get_xticklabels():
        tick.set_rotation(90)

    # Save the plot
    plt.savefig(f'{dirOut}/{gene}.length.png', format="png", dpi=300, bbox_inches='tight')
    plt.show()


def plot_summary(pri_data, alt_data, chr1_color, chr2_color, haploid, dirOut, gene, chr1, chr2):
    """
    Plots a summary evaluation
    Args:
    pri_data (pd.DataFrame): Primary haplotype data.
    alt_data (pd.DataFrame, optional): Alternate haplotype data. Defaults to None for haploid.
    chr1_color (str): Color for chromosome 1.
    chr2_color (str): Color for chromosome 2. Only used if not haploid.
    haploid (bool): Whether the data is haploid or not.
    dirOut (str): Output directory for saving the plot.
    gene (str): Gene name for labeling the plot.
    chr1 (str): Chromosome 1 label.
    chr2 (str, optional): Chromosome 2 label. Only used if not haploid.
    """
    # Font size settings
    plt.rc('font', size=9)
    plt.rc('axes', titlesize=9, labelsize=9)
    plt.rc('xtick', labelsize=9)
    plt.rc('ytick', labelsize=9)
    plt.rc('legend', fontsize=10)
    plt.rc('figure', titlesize=10)
    plt.rcParams["font.weight"] = "normal"
    plt.rcParams["axes.labelweight"] = "normal"

    # Set up subfigures and axes
    fig_size = (10, 5) if not haploid else (5, 5)
    SummPlt = plt.figure(figsize=fig_size)
    w = 0 if not haploid else 0.2
    subfigs = SummPlt.subfigures(3, 2, wspace=w, hspace=0.3)

    if not haploid:
        subfigsnest1 = subfigs[0][0].subplots(1, 2, sharey=True)
        subfigsnest2 = subfigs[0][1].subplots(1, 2, sharey=True)
        subfigsnest3 = subfigs[1][0].subplots(1, 2, sharey=True)
        subfigsnest4 = subfigs[1][1].subplots(1, 2, sharey=True)
        subfigsnest5 = subfigs[2][0].subplots(1, 2, sharey=True)
        subfigsnest6 = subfigs[2][1].subplots(1, 2, sharey=True)

        sns.scatterplot(x='start', y='mapping_quality', data=pri_data, color=chr1_color, s=10, ax=subfigsnest1[0], linewidth=0)
        sns.scatterplot(x='start', y='mismatch_rate', data=pri_data, color=chr1_color, s=10, ax=subfigsnest3[0], linewidth=0)
        sns.scatterplot(x='start', y='longindels', data=pri_data, color=chr1_color, s=10, ax=subfigsnest4[0], linewidth=0)
        sns.scatterplot(x='start', y='soft_clipping', data=pri_data, color=chr1_color, s=10, ax=subfigsnest5[0], linewidth=0)
        sns.scatterplot(x='start', y='hard_clipping', data=pri_data, color=chr1_color, s=10, ax=subfigsnest6[0], linewidth=0)
        
        sns.scatterplot(x='start', y='mapping_quality', data=alt_data, color=chr2_color, s=10, ax=subfigsnest1[1], linewidth=0)
        sns.scatterplot(x='start', y='mismatch_rate', data=alt_data, color=chr2_color, s=10, ax=subfigsnest3[1], linewidth=0)
        sns.scatterplot(x='start', y='longindels', data=alt_data, color=chr2_color, s=10, ax=subfigsnest4[1], linewidth=0)
        sns.scatterplot(x='start', y='soft_clipping', data=alt_data, color=chr2_color, s=10, ax=subfigsnest5[1], linewidth=0)
        sns.scatterplot(x='start', y='hard_clipping', data=alt_data, color=chr2_color, s=10, ax=subfigsnest6[1], linewidth=0)
        
        for f in (subfigsnest1, subfigsnest2, subfigsnest3, subfigsnest4, subfigsnest5, subfigsnest6):
            f[0].spines['right'].set_visible(False)
            f[0].spines['top'].set_visible(False)
            f[1].spines['right'].set_visible(False)
            f[1].spines['top'].set_visible(False)
            f[0].set_xlabel('Genomic Position')   
            f[1].set_xlabel('Genomic Position') 
    
        subfigsnest2[0].hist(pri_data['mapping_quality'], color=chr1_color)
        subfigsnest2[0].set_xlabel('Mapping Quality')
        subfigsnest2[0].set_ylabel('Frequency')
        subfigsnest2[1].hist(alt_data['mapping_quality'], color=chr2_color)
        subfigsnest2[1].set_xlabel('Mapping Quality')    
    
        subfigsnest1[0].set_ylabel('Mapping Quality') 
        subfigsnest3[0].set_ylabel('Mismatch Rate') 
        subfigsnest4[0].set_ylabel('Indels (>2)\nCounts') 
        subfigsnest5[0].set_ylabel('Soft Clipped\nBases Counts') 
        subfigsnest6[0].set_ylabel('Hard Clipped\nBases Counts') 
    
        # Bold font for tick labels
        plt.xticks(fontweight='normal')
        plt.yticks(fontweight='normal')
        custom = [Line2D([], [], marker='.', color=chr1_color, linestyle='None'),
              Line2D([], [], marker='.', color=chr2_color, linestyle='None')]
        plt.legend(custom, [f'{chr1}', f'{chr2}'], bbox_to_anchor=[2.3, 0], loc='lower right', frameon=False, markerscale=3)
    else:
        subfigsnest1 = subfigs[0][0].subplots(1, 1, sharey=True)
        subfigsnest2 = subfigs[0][1].subplots(1, 1, sharey=True)
        subfigsnest3 = subfigs[1][0].subplots(1, 1, sharey=True)
        subfigsnest4 = subfigs[1][1].subplots(1, 1, sharey=True)
        subfigsnest5 = subfigs[2][0].subplots(1, 1, sharey=True)
        subfigsnest6 = subfigs[2][1].subplots(1, 1, sharey=True)
    
        sns.scatterplot(x='start', y='mapping_quality', data=pri_data, color=chr1_color, s=10, ax=subfigsnest1, linewidth=0)
        sns.scatterplot(x='start', y='mismatch_rate', data=pri_data, color=chr1_color, s=10, ax=subfigsnest3, linewidth=0)
        sns.scatterplot(x='start', y='longindels', data=pri_data, color=chr1_color, s=10, ax=subfigsnest4, linewidth=0)
        sns.scatterplot(x='start', y='soft_clipping', data=pri_data, color=chr1_color, s=10, ax=subfigsnest5, linewidth=0)
        sns.scatterplot(x='start', y='hard_clipping', data=pri_data, color=chr1_color, s=10, ax=subfigsnest6, linewidth=0)
        
        for f in (subfigsnest1, subfigsnest2, subfigsnest3, subfigsnest4, subfigsnest5, subfigsnest6):
            f.spines['right'].set_visible(False)
            f.spines['top'].set_visible(False)
            f.set_xlabel('Genomic Position')   
    
        subfigsnest2.hist(pri_data['mapping_quality'], color=chr1_color)
        subfigsnest2.set_xlabel('Mapping Quality')
        subfigsnest2.set_ylabel('Frequency')    
        subfigsnest1.set_ylabel('Mapping Quality') 
        subfigsnest3.set_ylabel('Mismatch Rate') 
        subfigsnest4.set_ylabel('Indels (>2)\nCounts') 
        subfigsnest5.set_ylabel('Soft Clipped\nBases Counts') 
        subfigsnest6.set_ylabel('Hard Clipped\nBases Counts') 
    
        # Bold font for tick labels
        plt.xticks(fontweight='normal')
        plt.yticks(fontweight='normal')
        custom = [Line2D([], [], marker='.', color=chr1_color, linestyle='None'),
              Line2D([], [], marker='.', color=chr2_color, linestyle='None')]
        plt.legend(custom, [f'{chr1}'], bbox_to_anchor=[2.5, 0], loc='lower right', frameon=False, markerscale=3)
        
    #plt.tight_layout()
    plt.savefig(f'{dirOut}/{gene}.summary.allreads.png', format="png", dpi=300, bbox_inches='tight')
    plt.show()


def plot_coverage(
    positions, 
    coverage_counts, 
    mid_counts, 
    zero_counts, 
    start_indices, 
    end_indices, 
    high_mismatch_bool_size, 
    start_breaks, 
    end_breaks, 
    min_position, 
    max_position, 
    chr_label, 
    gene, 
    dirOut, 
    color60="#BBBBBB", 
    colormid="#EECC66", 
    color0="#BB5566"
):
    """
    Plots coverage by read mapping quality across a genomic region.

    Args:
    positions (array): Genomic positions.
    coverage_counts (array): Coverage counts for MapQ = 60.
    mid_counts (array): Coverage counts for MapQ = 1~59.
    zero_counts (array): Coverage counts for MapQ = 0.
    start_indices (array): Start indices for high mismatch regions.
    end_indices (array): End indices for high mismatch regions.
    high_mismatch_bool_size (int): Size of the high mismatch boolean array.
    start_breaks (array): Start indices of break regions.
    end_breaks (array): End indices of break regions.
    min_position (int): Minimum position in the region.
    max_position (int): Maximum position in the region.
    chr_label (str): Chromosome label.
    gene (str): Gene name for file naming and title.
    dirOut (str): Output directory for saving the plot.
    color60 (str): Color for MapQ = 60 reads. Default is "#BBBBBB".
    colormid (str): Color for MapQ = 1~59 reads. Default is "#EECC66".
    color0 (str): Color for MapQ = 0 reads. Default is "#BB5566".
    """

    # Plot setup
    plt.rcParams["font.weight"] = "bold"
    plt.rcParams["axes.labelweight"] = "bold"

    CovPlt, axes = plt.subplots(nrows=1, ncols=1, figsize=(20, 3), dpi=300)

    # Plot read-view high mismatch regions
    for start, end in zip(start_indices, end_indices):
        axes.axvspan(start, end if end != high_mismatch_bool_size else positions[-1], facecolor='rosybrown', alpha=0.2, zorder=0)

    # Plot break regions
    for start, end in zip(start_breaks, end_breaks):
        # if the break is too short, expand it so it is visible in figure
        if end - start < 5000:
            axes.axvspan(start - 4000, end + 4000, facecolor='purple', alpha=0.7)
        else:
            axes.axvspan(start, end, facecolor='purple', alpha=0.7)

    # Plot coverage for different mapping qualities
    axes.fill_between(positions, coverage_counts, step="pre", alpha=0.8, label='MapQ = 60', facecolor=color60, zorder=0.5)
    axes.fill_between(positions, mid_counts, step="pre", alpha=0.9, label='MapQ = 1~59', facecolor=colormid, zorder=0.5)
    axes.fill_between(positions, zero_counts, step="pre", alpha=0.5, label='MapQ = 0', facecolor=color0, zorder=0.5)

    # Plot details
    axes.set_title(f'Coverage by Reads in {chr_label}', fontweight='bold', size=12)
    axes.set_xlabel('Genomic Position')
    axes.set_ylabel('Coverage')
    axes.set_xlim(min_position, max_position)

    # Legend
    CovPlt.legend(loc="center right", bbox_to_anchor=(1.1, 0.5), frameon=False)

    # Save the plot
    plt.tight_layout()
    plt.savefig(f'{dirOut}/{gene}.{chr_label}.readcoverage.all.png', format="png", dpi=300, bbox_inches='tight')
    plt.show()



def plot_mismatch_coverage(pileup, bin_count, positions, start_indices, end_indices, 
                           high_mismatch_size, start_breaks, end_breaks, chr_color, 
                           chr_label, gene, dirOut, cmap, min_position, max_position):
    """
    Plots basepair level coverage (% mismatch per position) and a heatmap of poorly supported positions.

    Args:
    pileup (pd.DataFrame): DataFrame containing pileup data for positions.
    bin_count (pd.DataFrame): DataFrame with poorly supported bin counts.
    positions (array): Array of genomic positions for the data.
    start_indices (array): Start indices for high mismatch regions.
    end_indices (array): End indices for high mismatch regions.
    high_mismatch_size (int): Size of the high mismatch boolean array.
    start_breaks (array): Start indices of break regions.
    end_breaks (array): End indices of break regions.
    chr_color (str): Color for the read mismatch plot.
    chr_label (str): Chromosome label for the title.
    gene (str): Gene name for file naming and plot title.
    dirOut (str): Output directory for saving the plot.
    cmap (str): Colormap for the heatmap.
    """

# Set plotting options
    mpl.rcParams['agg.path.chunksize'] = 1000000000
    plt.rcParams["font.weight"] = "bold"
    plt.rcParams["axes.labelweight"] = "bold"

    # --- DOWNSAMPLING ---
    # Let's subsample the data if it's too large.
    
    num_points = len(pileup)
    max_plot_points = 1500000  # A safe upper limit 

    if num_points > max_plot_points:
        step = max(1, num_points // max_plot_points)
        binner = np.arange(num_points) // step
        print(f"Warning: Data has {num_points} points. Binning into "
              f"{len(np.unique(binner))} groups with median to avoid render errors.")
        plot_pileup = pileup.groupby(binner).agg(
            # For X-axis: Use the average position in the bin
            Pos=('Pos', 'mean'),  
            # For Y-axis: Use the median value in the bin
            PercentCorrect=('PercentCorrect', 'median') 
        ).reset_index(drop=True)
    else:
        plot_pileup = pileup

    # Set up subplots
    fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(20, 3), gridspec_kw={'height_ratios': [1, 3]})

    # Plot read-view high mismatch regions
    for start, end in zip(start_indices, end_indices):
        axes[1].axvspan(start, end if end != high_mismatch_size else positions[-1], facecolor='rosybrown', alpha=0.3, zorder=3)

    # Plot break regions
    for start, end in zip(start_breaks, end_breaks):
        if end - start < 5000:
            axes[1].axvspan(start - 4000, end + 4000, facecolor='purple', alpha=0.7, zorder=3)
        else:
            axes[1].axvspan(start, end, facecolor='purple', alpha=0.7, zorder=3)

    # Plot percentage of mismatches per position
    try:
        # Try the high-quality vector-heavy method first
        axes[1].fill_between(plot_pileup['Pos'], 100 - plot_pileup['PercentCorrect'], step="mid", 
                             alpha=1, color=chr_color, zorder=2, rasterized=True)
    except Exception as e:
        # If it fails (e.g., too many data points), fall back to the rasterized version
        print(f"Warning: 'fill_between' failed with error: {e}. "
              "Falling back to rasterized 'plot'.")
        axes[1].plot(
            plot_pileup['Pos'], 
            100 - plot_pileup['PercentCorrect'], 
            alpha=1, 
            color=chr_color, 
            zorder=2,
            rasterized=True,  # Use the rasterized option
            marker=None
        )

    # Create heatmap for poorly supported percentage
    data_matrix = np.array(100 - bin_count["wellCount percent"])[np.newaxis]
    try:
        sns.heatmap(data_matrix, ax=axes[0], annot=False, cbar=False, yticklabels=False, cmap=cmap, 
                    vmin=0, vmax=50, xticklabels=False)
    except Exception as e:
        axes[0].imshow(
                        data_matrix,
                        aspect='auto',
                        origin='lower',
                        vmin=0, vmax=50,
                        cmap=cmap,
                        interpolation='nearest'
                    )

    for spine in axes[0].spines.values():
        spine.set_visible(True)

    # Set titles and labels
    axes[0].set_title(f'Basepair view coverage (% of mismatch per position) in {chr_label}', 
                      fontweight='bold', size=12)
    axes[1].set_xlabel('Genomic Position')
    axes[1].set_ylabel('% of mismatch\nper position')
    axes[1].set_ylim(0, 101)
    axes[1].set_xlim(min_position, max_position)

    fig.tight_layout()
    axes[1].margins(x=0)

    # Save the plot
    try:
        fig.savefig(f'{dirOut}/{gene}.{chr_label}.basecoverage.PerCorrect.png', format="png", dpi=300, bbox_inches='tight')
    except Exception as e:
        print(f"Warning: 'savefig' with bbox_inches='tight' failed: {e}. Saving without it.")
        fig.savefig(f'{dirOut}/{gene}.{chr_label}.basecoverage.PerCorrect.png', format="png", dpi=300)
    
    plt.close(fig) # Close the figure to free up memory