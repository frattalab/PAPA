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

# Now double check that sample_name fields do not contain hyphens ('-') - will break column selection with R
if sample_tbl.index.str.contains("-", regex=False).any():
    raise Exception(f"Values in 'sample_name' column of sample table must not contain hyphen characters ('-')")

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
assert isinstance(config["use_precomputed_salmon_index"], bool), f"'use_precomputed_salmon_index' must be True/False boolean, {config['use_precomputed_salmon_index']} (type {type(config['use_precomputed_salmon_index'])}) was provided"

# Now double check that valid combination of use_provided_novel_les & run_identification is provided (if applicable)
if config["use_provided_novel_les"]:
    assert config["use_provided_novel_les"] and config["run_identification"], f"'use_provided_novel_les' is set to True but 'run_identification' is not. To continue using input novel last exons please set 'run_identification' to True"

# make a note if run_identification & use_precomputed_salmon_index are both True - run_id will be overrided and provided files used
if config["run_identification"] and config["use_precomputed_salmon_index"]:
    sys.stderr.write("run_identification will be overridden and pre-provided salmon index, id and info tables will be used\n")


# If differential, make sure that sample table only has two conditions & to set a contrast name
if config["run_differential"]:
    assert sample_tbl["condition"].nunique() == 2, f"condition column in sample table must only contain two distinct conditions, following n found - {sample_tbl['condition'].nunique()}"
    # firs key in sample table condition column = base_key
    BASE_KEY = sample_tbl["condition"][0]
    CONTRAST_KEY = sample_tbl.loc[sample_tbl["condition"] != BASE_KEY, "condition"][0]
    CONTRAST_NAME = CONTRAST_KEY + "vs" + BASE_KEY
    sys.stderr.write(f"Inferred base key for condition - {BASE_KEY}\n")
    sys.stderr.write(f"Inferred contrast key for condition - {CONTRAST_KEY}\n")
    sys.stderr.write(f"Constructed contrast name - {CONTRAST_NAME}\n")


include: "rules/filter_gtf.smk"
include: "rules/stringtie.smk"
include: "rules/tx_filtering.smk"
include: "rules/salmon.smk"
include: "rules/differential_apa.smk"

# sys.stderr.write(OPTIONS + "\n")

localrules: all, gtf_list_by_condition, gtf_list_all_tpm_filtered, check_per_sample_mean_tpm_filtered, make_formulas_txt

wildcard_constraints:
    sample = "|".join(SAMPLES),
    condition = "|".join(CONDITIONS)

rule all:
    input:
        rules.process_dexseq_tbl.output if config["run_differential"] else rules.tx_to_le_quant.output.ppau,
        rules.tx_to_le_quant.output.counts,
        os.path.join(DAPA_SUBDIR,
                     "summarised_pas_quantification.tpm.tsv"),
        os.path.join(DAPA_SUBDIR,
                     "summarised_pas_quantification.gene_tpm.tsv"),
        # expand(os.path.join(SALMON_SUBDIR, "quant", "{sample}", "quant.sf"),
        #        sample=SAMPLES,
        #        ),

