import pandas as pd
import os
import sys

configfile: "config/config.yaml"

include: "rules/parse_config.py"

sample_tbl = pd.read_csv(config["sample_tbl"], index_col="sample_name")

# Check fastq2, if all empty then dataset is single end
if sample_tbl['fastq2'].isna().all():
    single_end = True
else:
    single_end = False


SAMPLES = sample_tbl.index.tolist()
CONDITIONS = sample_tbl["condition"].drop_duplicates().tolist()
OPTIONS = sample_tbl.to_dict(orient="index")

GTF = config["annotation_gtf"]

# Make sure it has a slash at end of path
OUTPUT_DIR = os.path.join(config["main_output_dir"], "")
STRINGTIE_SUBDIR = os.path.join(OUTPUT_DIR, config["stringtie_subdir_name"], "")
TX_FILT_SUBDIR = os.path.join(OUTPUT_DIR, config["tx_filtering_subdir_name"], "")
SALMON_SUBDIR = os.path.join(OUTPUT_DIR, config["salmon_subdir_name"], "")
LOG_SUBDIR = os.path.join(OUTPUT_DIR, config["logs_subdir_name"], "")
BMARK_SUBDIR = os.path.join(OUTPUT_DIR, config["benchmarks_subdir_name"], "")
DAPA_SUBDIR = os.path.join(OUTPUT_DIR, config["diff_apa_subdir_name"], "")


## Double check switches for workflow control are valid

# First check all are True/False booleans
assert isinstance(config["run_identification"], bool), f"'run_identification' must be True/False boolean, {config['run_identification']} (type {type(config['run_identification'])}) was provided"
assert isinstance(config["run_differential"], bool), f"'run_differential' must be True/False boolean, {config['run_differential']} (type {type(config['run_differential'])}) was provided"
assert isinstance(config["filter_ref_gtf"], bool), f"'filter_ref_gtf' must be True/False boolean, {config['filter_ref_gtf']} (type {type(config['filter_ref_gtf'])}) was provided"
assert isinstance(config["use_provided_novel_les"], bool), f"'use_provided_novel_les' must be True/False boolean, {config['use_provided_novel_les']} (type {type(config['use_provided_novel_les'])}) was provided"

# Now double check that valid combination of use_provided_novel_les & run_identification is provided (if applicable)
if config["use_provided_novel_les"]:
    assert config["use_provided_novel_les"] and config["run_identification"], f"'use_provided_novel_les' is set to True but 'run_identification' is not. To continue using input novel last exons please set 'run_identification' to True"

include: "rules/filter_gtf.smk"
include: "rules/stringtie.smk"
include: "rules/tx_filtering.smk"
include: "rules/salmon.smk"
include: "rules/differential_apa.smk"

# sys.stderr.write(OPTIONS + "\n")

localrules: all, gtf_list_by_condition, gtf_list_all_tpm_filtered, check_per_sample_mean_tpm_filtered

wildcard_constraints:
    sample = "|".join(SAMPLES),
    condition = "|".join(CONDITIONS)

rule all:
    input:
        rules.process_saturn_tbl.output.processed_tbl if config["run_differential"] else rules.tx_to_le_quant.output.ppau,
        rules.tx_to_le_quant.output.counts,
        os.path.join(DAPA_SUBDIR,
                     "summarised_pas_quantification.tpm.tsv"),
        os.path.join(DAPA_SUBDIR,
                     "summarised_pas_quantification.gene_tpm.tsv"),
        # expand(os.path.join(SALMON_SUBDIR, "quant", "{sample}", "quant.sf"),
        #        sample=SAMPLES,
        #        ),



# def get_stringtie_assembled(sample, output_dir):
#     '''
#     Return path to target StringTie transcriptome assembly
#
#     Want functionality to provide a range of parameter values in same pipeline
#     and Snakemake's Paramspace docs aren't quite cutting it right now...
#
#     If provide a list for given parameter, will perform assembly for each combo of values
#     min_isoform_fraction_abundance (-f)
#     min_junction_reads (-j)
#     min_transcript_coverage (-c) (minimum reads per bp coverage)
#     To be added: disable_end_trimming (-t), point-features (--ptf)
#     '''
#
#     if isinstance(list(), config["min_isoform_fraction_abundance"]):
