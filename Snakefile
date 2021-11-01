import os
import pandas as pd

configfile: "config/config.yaml"

def param_list(param):
    '''
    Return list of all param values converted to string
    If param is not a list/iterable, coerced to a single value list
    '''

    try:
        param = list(param)
        out = [str(p) for p in param]

    except TypeError:
        # Not an iterable
        out = [str(param)]

    return out


sample_tbl = pd.read_csv(config["sample_tbl"], index_col="sample_name")

SAMPLES = sample_tbl.index.tolist()
CONDITIONS = sample_tbl["condition"].tolist()
OPTIONS = sample_tbl.to_dict(orient="index")

GTF = config["annotation_gtf"]

# Make sure it has a slash at end of path
OUTPUT_DIR = os.path.join(config["main_output_dir"], "")
STRINGTIE_SUBDIR = os.path.join(OUTPUT_DIR, config["stringtie_subdir_name"], "")
SALMON_SUBDIR = os.path.join(OUTPUT_DIR, config["salmon_subdir_name"], "")
LOG_SUBDIR = os.path.join(OUTPUT_DIR, config["logs_subdir_name"], "")
DAPA_SUBDIR = os.path.join(OUTPUT_DIR, config["diff_apa_subdir_name"], "")

min_frac_vals = param_list(config["min_isoform_fraction_abundance"])
min_jnc_vals = param_list(config["min_junction_reads"])
min_cov_vals = param_list(config["min_txipt_coverage"])

# print(OPTIONS)
# print(min_frac_vals)
# print(min_jnc_vals)
# print(min_cov_vals)

include: "rules/stringtie.smk"
include: "rules/salmon.smk"
include: "rules/tx_filtering.smk"
include: "rules/differential_apa.smk"


localrules: all, gtf_list_by_condition, gtf_list_all_tpm_filtered

wildcard_constraints:
    sample = "|".join(SAMPLES),
    condition = "|".join(CONDITIONS)

rule all:
    input:
        # expand(os.path.join(STRINGTIE_SUBDIR,
        #                     "min_jnc_{min_jnc}",
        #                     "min_frac_{min_frac}",
        #                     "min_cov_{min_cov}",
        #                     # "{condition}",
        #                     "{condition}.min_mean_tpm_filtered.gtf"),
        #        condition=CONDITIONS,
        #        min_jnc=min_jnc_vals,
        #        min_frac=min_frac_vals,
        #        min_cov=min_cov_vals),
        expand(os.path.join(SALMON_SUBDIR,
                             "pas_quant",
                             "min_jnc_{min_jnc}",
                             "min_frac_{min_frac}",
                             "min_cov_{min_cov}",
                             "summarised_pas_quantification.tsv"),
               min_jnc=min_jnc_vals,
               min_frac=min_frac_vals,
               min_cov=min_cov_vals),
        # expand(os.path.join(SALMON_SUBDIR, "quant", "min_jnc_{min_jnc}", "min_frac_{min_frac}", "min_cov_{min_cov}", "{sample}","quant.sf"),
        #        sample=SAMPLES,
        #        min_jnc=min_jnc_vals,
        #        min_frac=min_frac_vals,
        #        min_cov=min_cov_vals
        #        ),
        # expand(os.path.join(STRINGTIE_SUBDIR,
        #                     "min_jnc_{min_jnc}",
        #                     "min_frac_{min_frac}",
        #                     "min_cov_{min_cov}",
        #                     "tpm_filtered.intron_chain_filtered.3p_end_filtered.all_samples.combined.gtf"),
        #        min_jnc=min_jnc_vals,
        #        min_frac=min_frac_vals,
        #        min_cov=min_cov_vals)


def get_bam(sample, options, output_dir):
    '''
    Returns path to input bam file for given sample
    If sample will undergo additional processing (not yet implemented), path will be output_dir/<processing_step>/{sample}.bam
    If sample will not go additional processing, returns the path provided in the sample table/options

    params:
        sample <str>
        name of sample (in pipeline context should usually pass as wildcards.sample)

        options <dict>
        dict of {sample: {param1: x, param2: y}} generated from input sample table

        output_dir <str>
        path to main output directory (results for each sample stored within here)
    '''

    if config["pre_stringtie_processing"] == "none":
        return options[sample]["path"]

    else:
        raise ValueError("{} is invalid value for 'pre_stringtie_processing' option - please use 'none'".format(config["pre_stringtie_processing"]))

def get_sample_condition(sample, options):
    '''
    Return condition for given sample from options dict (sample table)
    '''

    return options[sample]["condition"]

def get_condition_samples(condition, options):

    return [sample for sample in options.keys() if options[sample]["condition"] == condition]



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
