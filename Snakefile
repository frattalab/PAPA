import pandas as pd
import os
import sys

configfile: "config/config.yaml"

include: "rules/parse_config.py"

sample_tbl = pd.read_csv(config["sample_tbl"], index_col="sample_name")

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

min_frac_vals = param_list(config["min_isoform_fraction_abundance"])
min_jnc_vals = param_list(config["min_junction_reads"])
min_cov_vals = param_list(config["min_txipt_coverage"])

include: "rules/filter_gtf.smk"
include: "rules/stringtie.smk"
include: "rules/tx_filtering.smk"
include: "rules/salmon.smk"
include: "rules/differential_apa.smk"

# sys.stderr.write(OPTIONS + "\n")
# sys.stderr.write(min_frac_vals + "\n")
# sys.stderr.write(min_jnc_vals + "\n")
# sys.stderr.write(min_cov_vals + "\n")


localrules: all, gtf_list_by_condition, gtf_list_all_tpm_filtered, check_per_sample_mean_tpm_filtered

wildcard_constraints:
    sample = "|".join(SAMPLES),
    condition = "|".join(CONDITIONS)

rule all:
    input:
        os.path.join(DAPA_SUBDIR,
                     "summarised_pas_quantification.counts.tsv"),
        os.path.join(DAPA_SUBDIR,
                     "summarised_pas_quantification.tpm.tsv"),
        os.path.join(DAPA_SUBDIR,
                     "summarised_pas_quantification.ppau.tsv"),
        os.path.join(DAPA_SUBDIR,
                     "summarised_pas_quantification.gene_tpm.tsv")
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
