# Modulation of the host cell transcriptome and epigenome by *Fusobacterium nucleatum*

May 10, 2021

Citation: [BioRxiv]



-----------------------------------------------------------

The included files contain code required to process the RNA-seq and ChIP-seq data.

## Software requirements

```
dependencies:
  - snakemake-minimal=5.4.5
  - samtools=1.9
  - rsem=1.3.2
  - star=2.5.2b
  - bam2fastq=1.1.0
  - sra-=2.10.8
  - openjdk=10.0.2
  - ChromHMM=1.21
  - python=3.7.8
```



## RNA-seq analysis

### Step 1:

[`rnaseq/01_rnaseq_to_gene_quantification.snakemake`](rnaseq/01_rnaseq_to_gene_quantification.snakemake) is a snakemake pipeline that takes the RNA-seq sample names (of compressed `.cram` files), decompresses them, extracts the fastq files, and runs RSEM to get gene expression quantification.

### Step 2:

The output of Step 1 is read into [`rnaseq/02_run_DESeq_and_annotate.Rmd`](rnaseq/02_run_DESeq_and_annotate.Rmd) to prepare and run DESeq2 on the different conditions. This Rmarkdown document also loads in the ensembl gene annotations and assigns gene names to the gene expression data.

### Step 3:

The output from Step 2 is read into [`rnaseq/03_plotting_and_GO_terms.Rmd`](rnaseq/03_plotting_and_GO_terms.Rmd) to generate the expression plots and perform the GO analysis.

## ChIP-seq analysis

### Step 1:

[`chipseq/01_chipseq_to_chromatin_states.snakemake`](chipseq/01_chipseq_to_chromatin_states.snakemake) is a snakemake pipeline that takes the deduplicated and aligned ChIP-seq datasets and runs the ChromHMM process, generating files for analysis and visualization in UCSC Genome Browser.

Note that between `rule compare_models` and `rule rename_states` (line ~230), manual intervention was performed to look at the results up to that point and select the 10-state model as the optimal, which was used for subsequent analysis and rules.

### Step 2:

The outputs from the above are read into [`chipseq/02_analyze_histone_and_chromatin_state_genomewide.R`](chipseq/02_analyze_histone_and_chromatin_state_genomewide.R) to perform the genome-wide analysis of levels of histone marks and chromatin states.

### Step 3:

The Fuso-consistent state files are read into [`chipseq/03_analyze_state_changes.R`](chipseq/03_analyze_state_changes.R) to perform the analysis on the state changes for each segment of the genome. These regions are written to a file (line 237), and [`chipseq/03a_breakup_chromhmm_regions_into_200bp_beds.py`](chipseq/03a_breakup_chromhmm_regions_into_200bp_beds.py) is used to split this bed file into bed files with a line for each 200bp segment, and at max 400,000 lines per file. This is a requirement for loading this data into GREAT to annotate the regions to genes. These files are loaded into http://great.stanford.edu/public/html/, and the gene-region associations are saved and loaded back into [`chipseq/03_analyze_state_changes.R`](chipseq/03_analyze_state_changes.R) (line 249). These are further filtered to just the basal regulatory region, and the Net Epigenomic Score is calculated for each gene. For the selected subsets of genes, the `random_test()` is performed to test how often a correlation between Net Epi Score and log2FC would be expected to be observed by chance, providing a random-resampling p-value.