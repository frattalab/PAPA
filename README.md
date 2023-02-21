# PAPA - Pipeline Alternative polyadenylation (APA)

Snakemake pipeline for detection & quantification of novel last exons/polyadenylation events.

The workflow (in brief) is as follows, but can be toggled depending on your use case (see below):
- Use [StringTie](https://github.com/gpertea/stringtie) to assemble transcripts from aligned RNA-seq reads
- Filter assembled transcripts for minimum mean expression across samples of the same condition ([Bioxriv preprint - Swamy et al., 2020](https://doi.org/10.1101/2020.08.21.261644))
- Filter for novel last exons that extend a known exon (optionally checking that the 5'ends match) or contain a novel last splice junction (optionally checking that the novel 5'ss matches a known 5'ss)
- Filter assembled transcripts for last exons with nearby reference polyA site (e.g. PolyASite) or conserved polyA signal motif (e.g. 18 defined in [Gruber et al., 2016](https://doi.org/10.1101/gr.202432.115))
- Merge novel last exon isoforms with reference annotation, quantify transcripts with Salmon
- Output count/TPM matrices via [tximport](https://doi.org/doi:10.18129/B9.bioc.tximport) for use with downstream differential transcript usage packages
- Perform differential usage analysis between two conditions using [satuRn](https://doi.org/10.18129/B9.bioc.satuRn)

Please consult [`docs/workflow.md`](docs/workflow.md) for a detailed description of the workflow


## Required input files & dependencies

The pipeline requires as input:
- **BAM files** of RNA-seq reads aligned to the genome (and corresponding BAI index files)
- **FASTQ files** of reads aligned to genome in corresponding BAM files (single-end/paired-end)
   - If you need to extract FASTQs from your BAM files, we have a [snakemake pipeline](https://github.com/frattalab/rna_seq_single_steps#Pull-FASTQs-from-BAM-files) that uses samtools extract FASTQs from a directory of BAM files
- **GTF file** of reference transcript models (sourced from **Gencode** or **Ensembl**, though Gencode is recommended as more extensively tested)
- **FASTA file** of genome sequence (and a corresponding FAI index file which can be generated using `samtools faidx`)
- **BED file** of reference poly(A) sites (e.g. [PolyASite 2.0](https://doi.org/10.1093/nar/gkz918))


The pipeline makes use of Snakemake & conda environments to install pipeline dependencies. If you do not have a conda installation on your system, head to the [conda website for installation instructions](https://conda.io/projects/conda/en/latest/user-guide/install/index.html). [mamba](https://mamba.readthedocs.io/en/latest/installation.html) can also be used as a drop in replacement (generally recommend this as it's much faster than conda!).

Assuming you have conda/mamba available on your system, the recommended way to run the pipeline is to use the 'execution' conda environment, which will install the required Snakemake version and other dependencies (including *mamba* for faster conda environment installation):

```
<conda/mamba> env create -f envs/snakemake_papa.yaml
```

Once installation is complete, you can activate the environment with the following command:

```
conda activate papa_snakemake
```


## Configuration

Note: the default config file is set up ready to run with test data packaged with the repository. If you'd like to see example outputs, you 


Sample information and input file relationships are defined in a CSV sample sheet with the following columns (an example can be found at `config/test_data_sample_tbl.csv`):

- `sample_name` - unique identifier for sample. Identifiers **must not** contain hyphen characters (i.e. '-') 
- `condition` - key to group samples of the same condition of interest (e.g. 'control' or 'knockdown')
- `path` - path to BAM file containing aligned reads for given sample. BAM files should be **coordinate sorted** and have a corresponding BAI index file at the same location (suffixed with '.bai')
- `fastq1` - path to FASTQ file storing 1st mates of read pairs aligned for same sample
- `fastq2` - path to FASTQ file storing 2nd mates of read pairs aligned for same sample

All pipeline parameters and input are declared in `config/config.yaml`, which needs to be edited per-run. All parameters are described in comments in the config file. The first section defines input file locations described above and the output directory destination. The 'General Pipeline Parameters' section declares how to modify filtering of assembled transcripts:

- '*min_mean_tpm*' - Minimum mean TPM across samples of the same condition for a transcript to be retained
- '*min_extension_length*' - Minimum length (nt) of 3'end extension relative to reference transcript for a novel extension isoform to be retained
- '*max_atlas_distance*' - Maximum distance between predicted 3'end and nearest reference poly(A) site for a transcript to be retained
- '*motif_search_region_length*' - Length of region from 3'end (to upstream) to search for presence of polyA signal motifs
- '*polya_signal_motifs*' - Which poly(A) signal motifs to use to retain 3'ends with evidence of genuine cleavage events taking place. 'Gruber' uses 18 conserved motifs identified by [Gruber et al., 2016](https://doi.org/10.1101/gr.202432.115), 'Beaudoing' uses 12 human motifs identified in [Beaudoing et al., 2000](https://doi.org/10.1101/gr.10.7.1001) (all of which found in Gruber). A set of custom PAS motifs (DNA, 1 per line) can also be provided via a TXT file.




## Running the pipeline



Once you are satisfied with the configuration, you should perform a dry run to check that Snakemake can find all the input files and declare a run according to your parameterisation. Execute the following command (with `papa_snakemake` conda environment active as described above):

```
snakemake -n -p --use-conda
```

If you are happy with the dry run, you can execute the pipeline locally with the following command, replacing <integer> with the number of cores you wish to use for parallel execution:

```
snakemake -p --cores <integer> --use-conda
```

Note: Conda environments will need to be installed the first time you run the pipeline. This can take a little while...


## Output

Output is organised into the following subdirectories:
```
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
- *differential_apa* - stores summarised quantification matrices for last exon isoforms across samples, along with satuRn's results table if differential usage analysis is performed.
    - satuRn results table - `test_data_output/differential_apa/saturn_apa.results.processed.tsv`
    - count matrices - `test_data_output/differential_apa/summarised_pas_quantification.<counts/gene_tpm/ppau/tpm>.tsv`
- *logs* - stores re-directed STDERR & STDOUT files for each job
- *salmon* - stores Salmon index & per-sample Salmon quantification outputs
- *stringtie* - stores per-sample StringTie output & expression filtered GTFs
- *tx_filtering* - stores outputs (including intermediates) of filtering transcripts for valid novel last exon isoforms, along with combined reference + novel GTF of last exons used for quantification.
    - All identified novel last exons - `test_data_output/tx_filtering/all_conditions.merged_last_exons.3p_end_filtered.gtf`
    - Combined GTF of reference + novel last exons - `test_data_output/tx_filtering/novel_ref_combined.last_exons.gtf`
    - GTF of last exons used for Salmon quantification - `test_data_output/tx_filtering/novel_ref_combined.quant.last_exons.gtf`


For a full description of output files, field descriptions etc. see [output_docs.md](output_docs.md)


## Advanced use-cases/parameters
