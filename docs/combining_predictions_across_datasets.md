# Workflow for combining PAPA predicted novel last exons across datasets

`scripts/combine_predicted_novel_last_exons.py` allows to combine the predicted novel last exons across multiple independent runs into a single dataset. The main benefit here is standardisation, so that a series of datasets can be quantified using a shared transcriptome reference.

`combine_predicted_novel_last_exons.py` takes a directory as input, and assumes that each subdirectory underneath it corresponds to a top-level PAPA output directory for a given experiment. Assuming you already have this directory structure, you can use `scripts/prepare_input_for_combine_predicted.sh` to generate a new directory with symbolic links to the required files (based on the PAPA output directory structure). The subdirectory name in the input is taken as the 'experiment name' identifier for the dataset

Note: as the script generates symbolic links, I recommend to provide an absolute path as the input directory argument to ensure the generated links will be valid

```bash
$ bash scripts/prepare_input_for_combine_predicted.sh -h
Error: Exactly two arguments required.
Usage: scripts/prepare_input_for_combine_predicted.sh [input_dir] [output_dir]

Create symbolic links to required files for combine_predicted_novel_last_exons.py across multiple PAPA identification runs for different datasets.

Arguments:
  input_dir    Path to directory containing experiment subdirectories
  output_dir   Path where links will be created (one subdirectory per experiment/subdirectory under input_dir)

Required files (automatically linked/assumed to exist in experiment subdirectories under input_dir):
  - tx_filtering/all_conditions.merged_last_exons.3p_end_filtered.gtf
  - tx_filtering/novel_ref_combined.tx2le.tsv
  - tx_filtering/novel_ref_combined.le2gene.tsv
  - tx_filtering/novel_ref_combined.le2genename.tsv
  - differential_apa/summarised_pas_quantification.ppau.tsv

Example: scripts/prepare_input_for_combine_predicted.sh /path/to/experiments /path/to/links

  - tx_filtering/all_conditions.merged_last_exons.3p_end_filtered.gtf
  - tx_filtering/novel_ref_combined.tx2le.tsv
  - tx_filtering/novel_ref_combined.le2gene.tsv
  - tx_filtering/novel_ref_combined.le2genename.tsv
  - differential_apa/summarised_pas_quantification.ppau.tsv
```

Once you've generated the directory structure, you can generate the combined GTF with `combine_predicted_novel_last_exons`

**WARNING**: the script requires at least pandas >=1.4. This is not satisfied using the PAPA conda environment. I tested this script using the `pybioinfo` conda environment in the accompanying [manuscript code repository](https://github.com/frattalab/tdp43-apa/tree/submission-02) (python 3.10.11, pandas 2.0.2, pyranges 0.0.127). 

```bash
$ python scripts/combine_predicted_novel_last_exons.py -h
usage: combine_predicted_novel_last_exons.py [-h] -i INPUT_DIR -o OUTPUT_GTF

Combine PAPA predicted novel last exons from multiple datasets into a single GTF

options:
  -h, --help            show this help message and exit
  -i INPUT_DIR, --input-dir INPUT_DIR
                        Path to directory containing experiment-wise subdirectories of GTFs of predicted novel last exons to merge, along with tx2le, le2genename, le2gene & PPAU tables to
                        additional metadata information. If multiple directories, pass paths consecutively and space-separated
  -o OUTPUT_GTF, --output-gtf OUTPUT_GTF
                        Name of output GTF file of combined last exons (default: None)
```

[You can then follow the instructions in the README to configure the pipeline to use this GTF as a set of input novel last exons to combine with a reference GTF](https://github.com/frattalab/PAPA?tab=readme-ov-file#2-use-a-pre-specified-set-of-novel-last-exons-to-combine-with-reference-last-exons-and-construct-salmon-index).
