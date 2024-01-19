# PAPA - Pipeline Alternative polyadenylation (APA)

Snakemake pipeline for detection & quantification of novel last exons/polyadenylation events from bulk RNA-sequencing data

The workflow (in brief) is as follows, but can be toggled depending on your use case:

- Use [StringTie](https://github.com/gpertea/stringtie) to assemble transcripts from aligned RNA-seq reads
- Filter assembled transcripts for minimum mean expression across samples of the same condition ([Biorxiv, Swamy et al., 2020](https://doi.org/10.1101/2020.08.21.261644))
- Filter for novel last exons that extend a known exon (optionally checking that the 5'ends match) or contain a novel last splice junction (optionally checking that the novel 5'ss matches a known 5'ss)
- Filter assembled transcripts for last exons with nearby reference polyA site (e.g. PolyASite) or conserved polyA signal motif (e.g. 18 defined in [Gruber et al., 2016](https://doi.org/10.1101/gr.202432.115))
- Merge novel last exon isoforms with reference annotation, quantify transcripts with Salmon
- Output count/TPM matrices via [tximport](https://doi.org/doi:10.18129/B9.bioc.tximport) for use with downstream differential transcript usage packages
- Perform differential usage analysis between two conditions using [DEXSeq](https://doi.org/doi:10.18129/B9.bioc.DEXSeq)

Please consult manuscript (*coming soon to biorxiv*) for a detailed description of the workflow

## Input files & dependencies

The pipeline requires as input:

- **BAM files** of RNA-seq reads aligned to the genome (and corresponding BAI index files)
- **FASTQ files** of reads aligned to genome in corresponding BAM files (single-end/paired-end)
  - If you need to extract FASTQs from your BAM files, we have a [snakemake pipeline](https://github.com/frattalab/rna_seq_single_steps#Pull-FASTQs-from-BAM-files) that uses samtools extract FASTQs from a directory of BAM files
- **GTF file** of reference transcript models (sourced from **Gencode** or **Ensembl**, though Gencode is recommended as it has been more extensively tested)
- **FASTA file** of genome sequence (and a corresponding FAI index file which can be generated using `samtools faidx`)
- **BED file** of reference poly(A) sites (e.g. [PolyASite 2.0](https://doi.org/10.1093/nar/gkz918))

Please see notes on [data compatibility](#compatible-data) for further information.

The pipeline makes use of Snakemake & conda environments to install pipeline dependencies. If you do not have a conda installation on your system, head to the [conda website for installation instructions](https://conda.io/projects/conda/en/latest/user-guide/install/index.html). [mamba](https://mamba.readthedocs.io/en/latest/installation.html) can also be used as a drop in replacement (and is generally recommended as it's much faster than conda!).

Assuming you have conda/mamba available on your system, the recommended way to run the pipeline is to use the 'execution' conda environment, which will install the required Snakemake version and other dependencies (including *mamba* for faster conda environment installation):

```bash
<conda/mamba> env create -f envs/snakemake_papa.yaml
```

Once installation is complete, you can activate the environment with the following command:

```bash
conda activate papa_snakemake
```

## Configuration

Note: the default config file is set up ready to run with test data packaged with the repository. If you'd like to see example outputs, you can skip straight to [running the pipeline]( #running-the-pipeline)

### Config YAML file

All pipeline parameters and input are declared in `config/config.yaml`, which needs to be customised for each run. All parameters are described in comments in the config file. The first section defines input file locations described above and the output directory destination (see comments in `config/config.yaml` for further details). See comments above each parameter for details.

The pipeline is modular and has several different run modes. Please see [workflow control section](#workflow-control) at the end of the README for further details.

### Sample sheet CSV

Sample information and input file relationships are defined in a CSV sample sheet with the following columns (an example can be found at `config/test_data_sample_tbl.csv`):

- `sample_name` - unique identifier for sample. Identifiers **must not** contain hyphen characters (i.e. '-')
- `condition` - key to group samples of the same condition of interest (e.g. 'control' or 'knockdown'). Note that first value in the first row is considered the 'base' condition/denominator for differential usage analysis.
- `path` - path to BAM file containing aligned reads for given sample. BAM files should be **coordinate sorted** and have a corresponding BAI index file at the same location (suffixed with '.bai')
- `fastq1` - path to FASTQ file storing 1st mates of read pairs aligned for same sample
- `fastq2` - path to FASTQ file storing 2nd mates of read pairs aligned for same sample. If dataset is single-end, leave this field blank

## Running the pipeline

Once you are satisfied with the configuration, you should perform a dry run to check that Snakemake can find all the input files and declare a run according to your parameterisation. Execute the following command (with `papa_snakemake` conda environment active as described above):

```bash
snakemake -n -p --use-conda
```

If you are happy with the dry run, you can execute the pipeline locally with the following command, replacing <integer> with the number of cores you wish to use for parallel execution:

```bash
snakemake -p --cores <integer> --use-conda
```

Note: Conda environments will need to be installed the first time you run the pipeline, which can take a while

## Output

Output is organised into the following subdirectories:

```bash
$ tree -d -L 1 test_data_output
test_data_output
├── benchmarks
├── differential_apa
├── logs
├── salmon
├── stringtie
└── tx_filtering
````

- *benchmarks* - stores the output of snakemake's `benchmark` directive (wall clock time, memory usage) for each job.
- *differential_apa* - stores summarised quantification matrices for last exon isoforms across samples, along with DEXSeq's results table if differential usage analysis is performed.
  - DEXSeq results table - `test_data_output/differential_apa/dexseq_apa.results.processed.tsv`
  - count matrices - `test_data_output/differential_apa/summarised_pas_quantification.<counts/gene_tpm/ppau/tpm>.tsv`
- *logs* - stores re-directed STDERR & STDOUT files for each job
- *salmon* - stores Salmon index & per-sample Salmon quantification outputs
- *stringtie* - stores per-sample StringTie output & expression filtered GTFs
- *tx_filtering* - stores outputs (including intermediates) of filtering transcripts for valid novel last exon isoforms, along with combined reference + novel GTF of last exons used for quantification.
  - All identified novel last exons - `test_data_output/tx_filtering/all_conditions.merged_last_exons.3p_end_filtered.gtf`
  - Combined GTF of reference + novel last exons - `test_data_output/tx_filtering/novel_ref_combined.last_exons.gtf`
  - GTF of last exons used for Salmon quantification - `test_data_output/tx_filtering/novel_ref_combined.quant.last_exons.gtf`

For a full description of output files, field descriptions etc. see [output_docs.md](docs/output_docs.md)

## Workflow control

The 'General Pipeline Parameters' declares how to control the workflow steps that are performed. The pipeline is split into modules - 'identification', 'quantification' and 'differential usage':

- Identification - controlled by `run_identification` parameter - controls whether to perform novel last exon discovery with StringTie + filtering.
- Differential usage - controlled by `run_differential` parameter - controls whether to perform differential usage between conditions with DEXSeq
- Quantification - performed regardless of identification or differential status. Last exon isoforms are quantified with Salmon and summarised to last-exon IDs with tximport

### Running modes when identification is not performed

When `run_identification` is set to False, a reference 'last exon-ome' must still be generated/provided as input to Salmon. There are a few possibilities:

#### 1. Just use reference GTF to construct Salmon index

- Set `use_provided_novel_les` and `use_precomputed_salmon_index` to False

#### 2. Use a pre-specified set of novel last exons to combine with reference last exons and construct Salmon index

- Set `use_provided_novel_les` to True
- Provide path to last exons in GTF format to `input_novel_gtf` parameter

#### 3. Use a pre-computed Salmon index (+ last exon metadata) from a previous PAPA run

- Set `use_precomputed_salmon_index` to True
- Provide path to directory containing pre-computed salmon index to `precomputed_salmon_index`
- Provide paths to computed last exon metadata to `tx2le`, `tx2gene`, `le2gene` and `info` parameters

## Compatible data

- If using single end data, then all samples in the sample sheet must be single end
- Unstranded data must not be used for novel last exon discovery, but can be used for quantification and differential usage with last exon reference from annotation or a pre-specified annotation.
- Otherwise, any bulk RNA-sequencing dataset should be compatible
- As long as the reference poly(A) site file is in BED file, in theory any database (e.g. PolyADB, Gencode manual annotation, custom annotation) can be used
  
  - Reference polyA site BED file can contain 'cluster' coordinates (i.e. not be a single nucloetide), but this can cause some duplication in predicted last exons within a condition. This should not affect quantification/differential usage severely (closely spaced polyA sites are annotated to the same last exon isoform), but can be make life more difficult if you want to extract coordinates for downstream analysis.
    - If multiple 3'ends overlap with a cluster interval, their coordinates will not be updated. This can mean you will have multiple closely spaced 3'end intervals predicted for the same event.
    - When clusters are provided, the nearby predicted last exons will have their 3'ends updated to the End coordinate of the cluster (i.e. not necessarily the most expressed/frequent coordinate within a cluster).


## License

Because Salmon is a core dependency and distributed under the [GNU General Public License v3.0 (GPLv3) licence](https://github.com/COMBINE-lab/salmon/blob/master/LICENSE), PAPA is also licensed under GPLv3 ([Open Source Guide](https://opensource.guide/legal/#which-open-source-license-is-appropriate-for-my-project)).
