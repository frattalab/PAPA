import os
import pandas as pd

configfile: "config/config.yaml"







sample_tbl = pd.read_csv(config["sample_table"], index_col="sample_name")

SAMPLES = sample_tbl.index.tolist()
OPTIONS = sample_tbl.to_dict(orient="index")

GTF = config["annotation_gtf"]

# Make sure it has a slash at end of path
OUTPUT_DIR = os.path.join(config["main_output_dir"],"")
STRINGTIE_SUBDIR = os.path.join(OUTPUT_DIR, config["stringtie_subdir_name"])


include: "rules/stringtie.smk"


localrules: all, compose_gtf_list_stringtie

rule all:
    input:
        expand(os.path.join(OUTPUT_DIR, config["stringtie_subdir_name"], "{sample}.assembled.gtf"), sample=SAMPLES),
        os.path.join(STRINGTIE_SUBDIR, "all_samples.merged.gtf")



def get_bam(sample, options, output_dir):
    '''
    Returns path to input bam file for given sample
    If sample will undergo additional processing (not yet implemented), path will be output_dir/aligned/{sample}.bam
    If sample will not go additional to not be realigned, returns the path provided in the sample table/options

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

    # if options[sample]["realign"] == 0 and options[sample]["file_type"] == "bam":
    #
    #     return options[sample]["path"]
    #
    # else:
    #
    #     return os.path.join(output_dir, config["bam_outdir_name"], sample + ".Aligned.sortedByCoord.out.bam")
