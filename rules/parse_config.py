'''
Collection of helper functions to parse config/sample table entries
'''

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


def parse_filter_flags(entry_dict, cl_flag):
    '''
    '''

    if entry_dict is None:
        return ""

    assert isinstance(entry_dict, dict), f"entry_dict should be dict - is type {type(entry_dict)}"

    if len(entry_dict) == 0:
        return ""


    optns = [key + "=" + ",".join(val) for key, val in entry_dict.items()]

    return [cl_flag] + optns


def parse_filter_tags(entry_list, cl_flag):
    '''
    '''

    assert isinstance(entry_list, list), f"entry_list should be list - is type {type(entry_list)}"

    if len(entry_list) == 0:
        return ""

    return cl_flag + " tag=" + ",".join(entry_list)
