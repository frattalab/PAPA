# PAPA - Pipeline Alternative polyadenylation (APA)

Snakemake pipeline for detection & quantification of novel last exons/polyadenylation events.

The workflow is as follows:
- Use [StringTie](https://github.com/gpertea/stringtie) to assemble transcripts from aligned RNA-seq reads
- Filter assembled transcripts for minimum mean expression across samples of the same condition
- Filter assembled transcripts for those with splice junctions matching any reference junction up until the penultimate junction (spliced last exons) or all junctions (extension last exons)
- Filter assembled transcripts for last exons with nearby reference polyA site (e.g. PolyASite) or conserved polyA signal motif (e.g. 18 defined in [Gruber et al., 2016](https://doi.org/10.1101/gr.202432.115))
- Merge novel last exon isoforms with reference annotation, quantify transcripts with Salmon
- Output count/TPM matrices via [tximport](https://doi.org/doi:10.18129/B9.bioc.tximport) for use with downstream differential transcript usage packages


## Required input files & dependencies

The pipeline requires as input:
- BAM files of RNA-seq reads aligned to the genome (and corresponding BAI index files)
- FASTQ files aligned to genome in corresponding BAM files (**MUST BE PAIRED-END**)
- GTF file of reference transcript models (sourced from **Gencode** or **Ensembl**, though Gencode is recommended as more extensively tested)
- FASTA file of genome sequence (and a corresponding FAI index file which can be generated using `samtools faidx`)
- BED file of reference poly(A) sites (e.g. [PolyASite 2.0](https://doi.org/10.1093/nar/gkz918))


The pipeline makes use of Snakemake & conda environments to install pipeline dependencies. If you do not have a conda installation on your system, head to the [conda website for installation instructions](https://conda.io/projects/conda/en/latest/user-guide/install/index.html).

Assuming you have conda available on your system, the recommended way to run the pipeline is to use the 'execution' conda environment, which will install the required Snakemake version and other dependencies (including *mamba* for faster conda environment installation):

```
conda env create -f envs/snakemake_papa.yaml
```

Once installation is complete, you can activate the environment with the following command:

```
conda activate papa_snakemake
```


## Configuration

Sample information and input file relationships are defined in a CSV sample sheet with the following columns (an example can be found at `config/local_ward_ipsc_example_sample_tbl.csv`):

- `sample_name` - unique identifier for sample
- `condition` - key to group samples of the same condition of interest (e.g. 'control' or 'knockdown')
- `path` - path to BAM file containing aligned reads for given sample. BAM files should be **coordinate sorted** and have a corresponding BAI index file at the same location (suffixed with '.bai')
- `fastq1` - path to FASTQ file storing 1st mates of read pairs aligned for same sample
- `fastq2` - path to FASTQ file storing 2nd mates of read pairs aligned for same sample

All pipeline parameters and input are declared in `config/config.yaml`, which needs to be edited per-run. All parameters are described in comments in the config file. The first section defines input file locations described above and the output directory destination. The 'General Pipeline Parameters' section declares how to modify filtering of assembled transcripts:

- '*min_mean_tpm*' - Minimum mean TPM across samples of the same condition for a transcript to be retained
- '*intron_chain_filter_mode*' - Whether assembled transcripts should have reference matching splice junctions sourced from the same transcript ('transcript' - **NOT YET IMPLEMENTED**) or any transcript ('any').
- '*max_terminal_non_match*' - (for novel spliced events) maximum number of consecutive 3'end non-reference matched introns/splice junctions for a transcript to be considered valid
- '*min_extension_length*' - Minimum length (nt) of 3'end extension relative to reference transcript for a novel extension isoform to be retained
- '*max_atlas_distance*' - Maximum distance between predicted 3'end and nearest reference poly(A) site for a transcript to be retained
- '*motif_search_region_length*' - Length of region from 3'end (to upstream) to search for presence of polyA signal motifs
- '*polya_signal_motifs*' - Which poly(A) signal motifs to use to retain 3'ends with evidence of genuine cleavage events taking place. 'Gruber' uses 18 conserved motifs identified by [Gruber et al., 2016](https://doi.org/10.1101/gr.202432.115), 'Beaudoing' uses 12 human motifs identified in [Beaudoing et al., 2000](https://doi.org/10.1101/gr.10.7.1001) (all of which found in Gruber. A set of custom PAS motifs (DNA, 1 per line) can also be provided via a TXT file
- '*expression_merge_by*' - How to summarise transcript abundance to polyA site isoforms for each gene. 'polyA' sums expression of transcripts if they share a poly(A) site, whereas 'last_exon' sums expression of transcripts if they share a last exon (any overlap; 3'UTR extension events will be considered as separate last exons).
- '*pas_merge_window_size*' - Target window size (nt), centred on poly(A) site) for merging closely spaced polyA sites to a single polyA cluster.


StringTie parameters can also be modified under the 'StringTie Parameters' section. For the following options, it is possible to set **multiple values**, in which case StringTie is run for **every possible combination of parameter values**:

- '*min_junction_reads*' - Minimum number of spliced reads aligning across a junction for it to be included in assembly process. Can be a float >= 1, passed as a single value (e.g. '1') or a list of values (e.g. '[1, 3, 5]')
- '*min_isoform_fraction_abundance*' - minimum fractional abundance relative to most abundance transcript at the locus for a transcript to be retained. Can be a single value between 0-1 passed as a single value (e.g. '0.01') or a list of values (e.g. '[0.01, 0.03, 0.05]')
- '*min_txipt_coverage*' - Minimum per-bp coverage for a predicted transcirpt to be retained. an be a float >= 1, passed as a single value (e.g. '1') or a list of values (e.g. '[1, 1.5]')



## Running the pipeline

Once you are satisfied with the configuration, you should perform a dry run to check that Snakemake can find all the input files and declare a run according to your parameterisation. Execute the following command (with `papa_snakemake` conda environment active as described above:

```
snakemake -n -p --use-conda
```

If you are happy with the dry run, you can execute the pipeline locally with the following command, replacing <integer> with the number of cores you wish to use for parallel execution:

```
snakemake -p --cores <integer> --use-conda
```

Note: Conda environments will need to be installed the first time you run the pipeline. This can take a little while...
