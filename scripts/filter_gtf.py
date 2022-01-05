import pyranges as pr
import numpy as np
import pandas as pd
import sys
from functools import reduce
from papa_helpers import eprint
import operator
import argparse
from timeit import default_timer as timer

# https://www.geeksforgeeks.org/python-key-value-pair-using-argparse/
# create a keyvalue class
# Args should be --less-than,--less-equals, --equals, --not-equals, --greater-equals, --greater-than
# all but equals & not equals should be coercable to int or float
# Function should do through each key in turn
class keyvalue(argparse.Action):
    # Constructor calling
    def __call__(self,
                 parser,
                 namespace,
                 values,
                 option_string=None):

        setattr(namespace, self.dest, dict())

        for value in values:
            # split it into key and value
            key, value = value.split('=')

            if "," in value:
                # Multiple values, need a list
                value = value.split(",")
            else:
                # just a single value, see if int/str
                try:
                    value = int(value)
                except ValueError:
                    # probably just a string, keep as is
                    pass

            # assign into dictionary
            getattr(namespace, self.dest)[key] = value


## Pass operators as 2nd argument for 'equality' subsets
# i.e. valid args are strings 'isin' or '~isin', or an operator (e.g. operator.eq())
# https://docs.python.org/3/library/operator.html


def _subset_membership(gr, col_name, filter_tuple):
    '''
    Subset PyRanges object over a specific column to retain rows with value in specified list/set

    gr: PyRanges object
    col_name: name of column in gr on which to apply filtering
    filter_tuple: 2-element tuple, 1st being list/set of values to check for membership, 2nd being one of 'isin' (retained if value in list/set) or '~isin' (retained if value not in list/set)
    '''

    assert isinstance(gr, pr.PyRanges)
    assert isinstance(filter_tuple[0], list) or isinstance(filter_tuple[0], set)
    assert filter_tuple[1] in ["isin", "~isin"]

    if filter_tuple[1] == "isin":
        return gr.subset(lambda df: df[col_name].isin(filter_tuple[0]))

    elif filter_tuple[1] == "~isin":
        return gr.subset(lambda df: ~df[col_name].isin(filter_tuple[0]))


def _subset_equality(gr, col_name, filter_tuple):
    '''
    Subset PyRanges object over a specific column for equality, greater than/less than operations

    gr: PyRanges object
    col_name: name of column in gr on which to apply filtering
    filter_tuple: 2-element tuple, first being value to check in column, 2nd being an comparison function from the operator base module
    e.g. _subset_gr(gr, "gene_type", ("protein_coding", operator.eq))
    '''
    assert isinstance(gr, pr.PyRanges)
    assert col_name in gr.columns

    return gr.subset(lambda df: filter_tuple[1](df[col_name], filter_tuple[0]))


def _subset_str_contains(gr, col_name, filter_tuple):
    '''
    Subset PyRanges object over a specific column for equality, greater than/less than operations

    gr: PyRanges object
    col_name: name of column in gr on which to apply filtering
    filter_tuple: 2-element tuple, first being value to check in column, 2nd being an comparison function from the operator base module
    e.g. _subset_gr(gr, "gene_type", ("protein_coding", operator.eq))
    '''

    assert isinstance(gr, pr.PyRanges)
    assert col_name in gr.columns
    assert filter_tuple[1] in ["contains", "~contains"]

    val = filter_tuple[0]

    if isinstance(val, list):
        val = "|".join(val)

    if filter_tuple[1] == "contains":
        return gr.subset(lambda df: df[col_name].str.contains(val, na=False))

    elif filter_tuple[1] == "~contains":
        return gr.subset(lambda df: ~df[col_name].str.contains(val, na=False))


def _subset_gr(gr, col_name, filter_tuple):
    '''
    '''

    assert isinstance(gr, pr.PyRanges)
    assert col_name in gr.columns

    if isinstance(filter_tuple[1], str):
        assert filter_tuple[1] in ["isin", "~isin", "contains", "~contains"]

        if filter_tuple[1] in ["isin", "~isin"]:
            return _subset_membership(gr, col_name, filter_tuple)

        else:
            # str contains method wanted
            return _subset_str_contains(gr, col_name, filter_tuple)

    else:
        # filter_tuple[1] should be an operator
        return _subset_equality(gr, col_name, filter_tuple)


def _subset_gr_and(gr,
                   col_name_1,
                   filter_tuple_1,
                   col_name_2,
                   filter_tuple_2):
    '''
    '''

    assert isinstance(gr, pr.PyRanges)
    assert col_name_1 in gr.columns
    assert col_name_2 in gr.columns

    # Always matching for a single value
    # Filter down to rows matching 'top_level' key
    gr2 = _subset_equality(gr, col_name_1, filter_tuple_1)

    if isinstance(filter_tuple_2[1], str):
        assert filter_tuple_2[1] in ["isin", "~isin", "contains", "~contains"]
        return _subset_membership(gr2, col_name_2, filter_tuple_2)

    else:
        # filter_tuple_2[1] should be an operator
        return _subset_equality(gr2, col_name_2, filter_tuple_2)


def filter_gtf_attributes(gr, gene_types, tr_types=None):
    '''
    Filter PyRanges with two-level, group-specific filtering schemes given column/attribute values

    gr: PyRanges object

    gene_types: dict of {'column_name': (value, return_true)}
        - value - can be a string, numeric (int or float) or list type
        - return_true - specifies how elements in 'column_name' should relate to value in order for them to be returned
            - return_true can be an operator (from operator base module) or
            one of ['isin', '~isin'] (to test for membership)

    tr_types: optional dict of {'subgroup_key': {'column_name': (value, return_true)}} for value-specific filtering of values in gene_types {col: (['a','b'], 'isin')}
        e.g. if my gene_types dict had {col: (['a','b'], 'isin')}
        You can specify additional filters for group 'a' (and/or 'b') in tr_types
        keys of tr_types must be present in a list/set defined in 1st element of gene_types values (the filter_tuple)


    e.g. I want to extract 'protein_coding' & 'lncRNA' gene_types ('gene_type' GTF attribute key)
    and 'transcript' fields only ('Feature' key (3rd field in GTF))
    for 'protein_coding' gene types, I want to extract TSL 1 transcripts ('transcript_support_level' attribute key)
    filter_gtf_attributes(gr,
                          gene_types={'Feature': ('transcript', operator.eq),
                                      'gene_type': (['protein_coding', lncRNA], 'isin'),
                                      },
                          tr_types={'protein_coding': {'transcript_support_level': (1, operator.eq),
                                                       }
                                    }
                          )

    Notes:
    - Not every key in gene_types needs to have additional filters defined
    - tr_types only becomes active if at least 1 of tuple[0] (value) in gene_types.values() is a list/set.
        - If you have multiple filters but no 'group-specific' filters, put all of your filters in gene_types dict
    - Filters are applied iteratively (i.e. consider multiple filters in the same dict as 'AND' conditions)


    '''
    assert isinstance(gr, pr.PyRanges)
    assert isinstance(gene_types, dict)

    # Apply each 'top-level' filter defined in gene_types
    gr2 = reduce(lambda gr, col: _subset_gr(gr, col, gene_types[col]),
                 gene_types.keys(),
                 gr)

    # Optionally apply key-specific filters
    if tr_types is None:
        return gr2

    else:
        assert isinstance(tr_types, dict)

        # Get a set of filters from gene_types that contain lists
        gene_types_l = {key: val for key,val in gene_types.items()
                        if isinstance(val[0], list) or isinstance(val[0], set)
                        }

        # Double check can find group-specific filter in gene_types dict ('top level' filter)
        # Otherwise will be filtering for keys that don't exist in gr
        # Also create dict of {'tr_type': 'gene_type_col'}
        tt2gene_col = {}

        for key in tr_types.keys():
            found = False
            for gt_col, filter_tuple in gene_types_l.items():
                if key in filter_tuple[0]:
                    found = True
                    tt2gene_col[key] = gt_col

                else:
                    continue

            if not found:
                raise Exception(f"tr_types subgroup key - {key} - not found in gene_types dict. To undergo group specific filtering must filter for this group first")

        # print(tr_type_keys)
        gr_list = []

        # 1st - find set of top level keys (gene_types) that are not sent for additional filtering
        # i.e. values in lists of gene_types tuples that aren't in tr_types.keys()
        # Subset gr for these values and add to gr_list to save in final output
        for col, filter_tuple in gene_types_l.items():
            for ele in filter_tuple[0]:
                if ele in tr_types.keys():

                    continue

                else:
                    # won't undergo extra filtering, need to save unmodified
                    save_gr = _subset_gr(gr2, col, (ele, operator.eq))
                    gr_list.append(save_gr)


        # Now perform group-specific filters
        for key, filter_dict in tr_types.items():
            # Subset for group rows defined in top-level (gene_type) filter
            gr_tt = _subset_gr(gr2,
                               tt2gene_col[key],
                               (key, operator.eq)
                               )

            gr_tt_f = reduce(lambda gr, col: _subset_gr(gr, col, filter_dict[col]),
                             filter_dict.keys(),
                             gr_tt)

            gr_list.append(gr_tt_f)


        return pr.concat(gr_list)


def main(gtf_path, filters_dict, out_path):
    '''
    '''

    eprint("Reading in input GTF file, this can take a while...")

    r_start = timer()
    gtf = pr.read_gtf(gtf_path, duplicate_attr=True)
    r_end = timer()

    eprint(f"Reading complete - took {r_end - r_start} s")

    eprint("Checking for numerics in filters_dict & updating corresponding column types...")
    for col, filter_tup in filters_dict.items():
        if isinstance(filter_tup[0], int) or isinstance(filter_tup[0], float):
            eprint(f"Converting {col} to numeric in GTF PyRanges object...")

            gtf = gtf.assign(col,
                             lambda df: df[col].replace("NA", np.nan).astype(float))

    eprint("Filtering GTF file...")
    gtf = filter_gtf_attributes(gtf, filters_dict)

    eprint(f"Writing filtered gtf to file at {out_path}...")
    gtf.to_gtf(out_path)


if __name__ == '__main__':

    start = timer()

    descrpn = """Script to filter reference GTF based on presence/absence of attributes"""

    parser = argparse.ArgumentParser(description=descrpn,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter, # Add defaults to end of help strings
                                     )

    # Args should be --less-than,--less-equals, --equals, --not-equals, --greater-equals, --greater-than

    parser.add_argument("-i",
                        "--input-gtf",
                        required=True,
                        type=str,
                        default=argparse.SUPPRESS,
                        help="Path to input GTF file")

    parser.add_argument("-t",
                        "--min-tsl",
                        type=int,
                        default=None,
                        choices=list(range(1,6,1)),
                        help="Minimum TSL level for a transcript to be retained (1-5 where 1 = highest support). Optional filter, include to activate"
                        )

    # parser.add_argument("-g",
    #                     "--gene-types",
    #                     type=str,
    #                     default=argparse.SUPPRESS,
    #                     action=keyvalue,
    #                     help="Filter to retain all entries of specified gene type (e.g. 'protein_coding'). Pass as <attribute_name>=<value1>,<value2> e.g. gene_type=protein_coding,lncRNA")

    parser.add_argument("--include-flags",
                        default=argparse.SUPPRESS,
                        nargs="*",
                        action=keyvalue,
                        help="Filter to retain entries matching all provided conditions. Pass as <attribute_name>=<value>, separated by space if multiple and values comma separated if multiple")

    parser.add_argument("--exclude-flags",
                        default=argparse.SUPPRESS,
                        nargs="*",
                        action=keyvalue,
                        help="Filter to exclude entries matching all provided conditions. Pass as <attribute_name>=<value>, separated by space if multiple and values comma separated if multiple")

    parser.add_argument("--tag-include",
                        nargs="*",
                        action=keyvalue,
                        help="Filter to include entries based on conditions operating on the 'tag' attribute column. GTF files can have duplicate tag keys, PyRanges will concatenate values (separated by ',') so need bespoke matching approach")

    parser.add_argument("--tag-exclude",
                        nargs="*",
                        action=keyvalue,
                        help="Filter to exclude entries based on conditions operating on the 'tag' attribute column. GTF files can have duplicate tag keys, PyRanges will concatenate values (separated by ',') so need bespoke matching approach")

    parser.add_argument("-o",
                        "--output-gtf",
                        type=str,
                        default="filtered_gtf.gtf",
                        help="Path to/name of output filtered GTF file")

    if len(sys.argv) == 1:
        parser.print_help()
        parser.exit()

    args = parser.parse_args()
    # Get a dict from NameSpace object
    args = vars(args)

    # print(args)
    # Need to format parsed filters into {col: (val, operator)}

    # filter_args = ["min_tsl", "include_flags", "exclude_flags"]
    filters_dict = {}
    for arg_name, val in args.items():

        if arg_name == "min_tsl":
            if val is not None:
                filters_dict["transcript_support_level"] = tuple([val, operator.le])
            else:
                continue

        elif arg_name == "tag_include":
            # This is a dict of {col: val}
            if val is not None:
                for inc_col, inc_val in val.items():
                    filters_dict[inc_col] = tuple([inc_val, "contains"])

            else:
                continue

        elif arg_name == "tag_exclude":
            # This is a dict of {col: val}
            if val is not None:
                for inc_col, inc_val in val.items():
                    filters_dict[inc_col] = tuple([inc_val, "~contains"])

            else:
                continue

        elif arg_name == "include_flags":
            # This is a dict of {col: val}
            if val is not None:
                for inc_col, inc_val in val.items():
                    if isinstance(inc_val, list):
                        # multiple values
                        filters_dict[inc_col] = tuple([inc_val, "isin"])

                    else:
                        # just a single val
                        filters_dict[inc_col] = tuple([inc_val, operator.eq])

            else:
                continue

        elif arg_name == "exclude_flags":
            # This is a dict of {col: val}
            if val is not None:
                for inc_col, inc_val in val.items():
                    if isinstance(inc_val, list):
                        # multiple values
                        filters_dict[inc_col] = tuple([inc_val, "~isin"])

                    else:
                        # just a single val
                        filters_dict[inc_col] = tuple([inc_val, operator.ne])

            else:
                continue

        else:
            # option is not a filter
            continue

    eprint("Formatted filters from command line input")
    eprint(filters_dict)

    main(args["input_gtf"], filters_dict, args["output_gtf"])

    end = timer()

    eprint(f"Complete: took {round(end - start, 3)} s / {round((end - start) / 60, 3)} min (3 dp)")
