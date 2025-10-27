# IG Assembly Evaluation with CloseRead

**Brief description:** This tool aims to evaluate the assembly of immunoglobulin (IG) sequences. The evaluation process includes various scripts and workflows to analyze, compare, and validate IG sequence assemblies.

![pipeline overview](plots/overview.png)

## Table of Contents

- [Installation](#installation)
- [Requirements](#requirements)
- [Test Case](#test-case)
- [Input Files](#input-files)
- [Usage](#usage)
- [Outputs](#outputs)
- [Directory Structure](#directory-structure)
- [Citation](#citation)
- [License](#license)
- [Contributing](#contributing)

---

## Installation

Follow these steps to set up the development environment:

```bash
# Clone the repository
git clone https://github.com/phylo-lab-usc/CloseRead_ig-assembly-eval.git

# Navigate to the project directory
cd CloseRead_ig-assembly-eval

# Create and activate the conda environments
conda env create -f ig-assembly-eval.yml
conda activate ig-assembly-eval

# Install closeread
pip install .

```
> **Tip:** Optionally, update your system PATH if needed.

## Requirements
Before running CloseRead, ensure the following tools are installed and configured correctly:

- **IgDetective:**  
  - [GitHub Repository](https://github.com/Immunotools/IgDetective.git)  
  - **Important:** To avoid errors, modify the IgDetective source code so that any relative paths are replaced with absolute paths. Also make sure you are running IgDetective with the python in the above `ig-assembly-eval` enviroment. 
- **samtools:** Required for various downstream analyses.
- **Python 3.10:** This version is included in the provided conda environment. Confirm by running:
  ```bash
  which python
  ```

---


## Test Case
Before analyzing your own data, run the test case to verify that CloseRead is installed correctly. Replace the placeholders (`{PATH_to_IgDetective}`) with your actual paths:
```
cd CloseRead_ig-assembly-eval
closeread-pipeline --species test --home test/ --haploid True --fastqdir hifi_fastq --igdetective_home {PATH to IgDetective}
closeread-plot --s test --g IGH --home test/ --pg test/igGene/test.pri.igdetective/combined_genes_IGH.txt
```
> **Note:** On a 32-core system, the entire process should finish in about 5 minutes. Compare your outputs with the expected results available [here](https://figshare.com/s/9b8db110af1511871669). For guidance on interpreting the results, refer to [this document](https://docs.google.com/document/d/1QOh3Z6noqZ7x-u70QhQv4VhQOrCJJ1hz_NEmvBDaT_A/pub).

## Input Files

#### Required Files:
Place the following files in the specified directories (replace placeholders with your species name):
- **HiFi Fastq/Fastq.gz/Fq/BAM Files:**  
  Located at:  
  ```
  $HOME/$fastqdir/$species_name/
  ```

- **Merged Assembly FASTA:**  
  ```
  $HOME/assemblies/${species_name}.merged.fasta
  ```  
  *Use the entire assembly even if you’re focusing on specific loci.*
  *If you have a phased assembly, concatenate the primary/hapotype1 assembly with the alternate/haplotype2 assembly into 1 file*

- **Primary/Haplotype1 Assembly FASTA:**  
  ```
  $HOME/assemblies/${species_name}.pri.fasta
  ```
  *If you only have a unphased consensus assembly, this is just a copy of the merged.fasta file*

- **Alternate/Haplotype2 Assembly FASTA:**  
  ```
  $HOME/assemblies/${species_name}.alt.fasta
  ```  
  *If you only have a unphased consensus assembly, you can ignore this file.*

- **Assembly Index Files:**  
  Ensure that each assembly FASTA file has a corresponding `.fai` index file.

> **Optional:** If IG loci positions are already known (allowing you to skip running IgDetective), include the custom annotation file named `${species_name}.customIG.txt`. (See the example file in `example/Emax.customIG.txt`. And more on file format see [`docs/format.md`](docs/format.md).)


#### Optional Files

- **Species Meta-information:**  
   This is for the final pdf generation only, if individual png file is preferred, this can be skipped. This is a CSV file (`species_metainfo.csv`) containing meta-information about your species. See the `example` folder for a sample file. And more on file format see [`docs/format.md`](docs/format.md). 

- **Gene-level Annotation:**  
  This can be either in the IgDetective-generated format **or** as a CSV with the columns: `Gene, Chromosome, Strand, Start, End` (include a header). See the `example` folder for a sample file. And more on file format see [`docs/format.md`](docs/format.md).

---

## Usage

#### Running read-to-assembly Pipeline

Use the `closeread-pipeline` command to perform the initial read-to-assembly evaluation. Below is the basic usage:
```bash
closeread-pipeline [OPTIONS]
usage: closeread-pipeline [-h] --species SPECIES --home HOME --haploid HAPLOID --fastqdir FASTQDIR
                          (--igdetective_home IGDETECTIVE_HOME | --customIG CUSTOMIG) [--t T]
```

#### Key Options:
- `--species`: Comma-separated list of species (e.g., `species1,species2`).
- `--home`: Path to the home directory.
- `--haploid`: `True` for consensus assembly or `False` for phased assembly.
- `--fastqdir`: Name of FASTQ directory.
- Either:
  - `--igdetective_home`: Path to the IgDetective directory, **or**
  - `--customIG`: Path to the DIRECTORY containing `${species_name}.customIG.txt` (NOT path directly to the file).
- `--t`: Number of threads (default: 32).

#### Running the Evalutaion and Plotting Pipeline
After running the evaluation pipeline, generate visualizations with `closeread-plot`. Run separately for each locus. For gene-level assessments, include the `--pg` parameter.

```bash
closeread-plot [OPTIONS]
usage: closeread-plot [-h] (--s species | --sf species_file) --g gene --home home [--cov lowCov_threshold]
                      [--p padding] [--re single_read_error] [--rc readview_correct_threshold] [--bc baseview_correct_threshold] [--m meta]
                      [--so stats_only] [--pg gene_level assessment]
```
#### Key Options:
- `--s`: Single species identifier (use this if providing one species) or `--sf` for a file containing a list of species.
- `--g`: Gene identifier (e.g., IGH, IGK, IGL).
- `--home`: Path to the home directory.
- Additional optional parameters:
  - `--cov`: Low coverage threshold (default: 2)
  - `--p`: Purple bar padding around low coverage regions (default: 2000 bp)
  - `--re`: Single-read mismatch threshold (default: 0.01)
  - `--rc`: Number of high-mismatch reads required to flag a position from the read-view (default: 5)
  - `--bc`: Base-view threshold (default: 80% exact match)
  - `--m`: Path to the meta-information CSV file (for generating a PDF report)
  - `--so`: Generate statistics only (skip visualization)
  - `--pg`: Path for a gene-level annotation file (to produce gene-level assessments)


## Outputs
#### From `closeread-pipeline`
The pipeline produces intermediate statistics files in `{HOME}/errorStats/`. Typical outputs include:

- **IG Loci Files:**
  - `IGH.txt`: Read-oriented statistics for both primary and alternate assemblies at the IGH locus.
  - `IGH_alt_pileup.txt`: Basepair-level (mpileup) statistics for the alternate assembly at IGH.
  - `IGH_pri_pileup.txt`: Basepair-level (mpileup) statistics for the primary assembly at IGH.
  - Similar files for `IGK` and `IGL`.
  
- **Non-IG Loci:**
  - `nonIG.txt`: Read-oriented statistics for non-IG loci.
  
- **Flag Files:**
  - `cigar.end` and `pileup.end`: Empty files signaling the completion of certain steps.

Other outputs include:
- **Aligned Files:**
  - BAM files in `{HOME}/aligned_bam/{species}/` (both full and filtered (primary-only) alignments).
- **IgDetective Results:**
  - Stored in `{HOME}/igGene/{species}.{pri/alt}.igdetective/`.
- **IG Loci Positions:**
  - IG loci positions in `{HOME}/gene_position/{species}.final.Ig_loci.txt` (or NA if using `--customIG`).


#### From `closeread-plot`
Final visualization files are saved in `{HOME}/errorPlots/` and typically include:

- `summary.allreads.png`, summary evaluation of the general stats
- `length.png`, locus length plot
- `readcoverage.all.png`, coverage by read mapping quality across a genomic region, for both haplotype if diploid
- `basecoverage.PerCorrect.png`, basepair level coverage (% mismatch per position) and a heatmap of poorly supported positions, for both haplotype if diploid
- `break.txt` containing information for all low coverage positions
- `read.mismatch.txt` containing positional information for all high-mismatch rate region from read-view
- `base.exactmismatch.csv` containing detailed basepair-view mismatch information for every basepair with < bc (80% default) base support
- `base.avgmismatch.csv` containing averaged basepair-view mismatch information, grouping consecutive basepairs (< 1000bps spaced), best for identifying regions with lots of single base mismatches
- `finalMismatch.csv` containing positions marked highly mismatched from both basepair-view and read-view, for both haplotype if diploid
- `result.pdf`, single pdf with all figures and meta information, if meta info given
- `{species}.{gene}.genelevel.csv`, gene level quality assessment, if flag `--pg` given

#### For more information on how to interpret the result please refer to this [document](https://docs.google.com/document/d/1QOh3Z6noqZ7x-u70QhQv4VhQOrCJJ1hz_NEmvBDaT_A/pub) 

#### Logs

All log files are stored in:
```
$HOME/logs/
```

---

## Directory Structure

A typical working directory might look like this:

```plaintext
$HOME/
├── assemblies/          # Assembly input files
│   ├── mEubGla1.merged.fasta
│   ├── mEubGla1.merged.fasta.fai
│   ├── mEubGla1.pri.fasta
│   ├── mEubGla1.pri.fasta.fai
│   ├── mEubGla1.alt.fasta
│   ├── mEubGla1.alt.fasta.fai
├── $fastqdir/           # FASTQ/BAM files directory
│   └── mEubGla1/
├── aligned_bam/         # Aligned BAM files
│   └── mEubGla1/
│       ├── mEubGla1_merged_sorted.bam
│       ├── mEubGla1_merged_sorted.bam.csi
│       ├── mEubGla1_merged_sorted_primary.bam
│       └── mEubGla1_merged_sorted_primary.bam.csi
├── igGene/              # IgDetective results
│   ├── mEubGla1.pri.igdetective/
│   └── mEubGla1.alt.igdetective/    # For diploid assemblies
├── gene_position/       # IG loci positions
│   └── mEubGla1.final.Ig_loci.txt   # or mEubGla1.customIG.txt
├── errorStats/          # Mismatch and coverage statistics
│   └── mEubGla1/
├── errorPlots/          # Final visualizations
│   └── mEubGla1/
└── logs/                # Log files
```

---

## Citation

If you use CloseRead in your research, please cite:

**Assessing Assembly Errors in Immunoglobulin Loci: A Comprehensive Evaluation of Long-read Genome Assemblies Across Vertebrates**  
Yixin Zhu, Corey Watson, Yana Safonova, Matt Pennell, Anton Bankevich  
bioRxiv 2024.07.19.604360; doi: [10.1101/2024.07.19.604360](https://doi.org/10.1101/2024.07.19.604360)

---

## License

This project is licensed under the [GNU General Public License v3 (GPLv3)](LICENSE).

---

## Contributing

We loved to hear feedback on this project! Before starting, please check out our [CONTRIBUTING.md](./CONTRIBUTING.md) file for guidelines on how to help improve CloseRead effectively.

> **Note:** We prefer issues to be submitted first before creating pull requests. This ensures alignment with the project goals and avoids duplicate efforts.
