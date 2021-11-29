import pyranges as pr
import numpy as np
import pandas as pd
import sys
from functools import reduce
import operator

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


def _subset_gr(gr, col_name, filter_tuple):
    '''
    '''

    assert isinstance(gr, pr.PyRanges)
    assert col_name in gr.columns

    if isinstance(filter_tuple[1], str):
        assert filter_tuple[1] in ["isin", "~isin"]
        return _subset_membership(gr, col_name, filter_tuple)

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
        assert filter_tuple_2[1] in ["isin", "~isin"]
        return _subset_membership(gr2, col_name_2, filter_tuple_2)

    else:
        # filter_tuple_2[1] should be an operator
        return _subset_equality(gr2, col_name_2, filter_tuple_2)


def filter_gtf_attributes(gr, gene_types,tr_types=None):
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


if __name__ == '__main__':
    p_gtf = sys.argv[1]

    print("Reading in GTF...")

    gtf = pr.read_gtf(p_gtf, duplicate_attr=True)

    # print("Trying filtering just for protein-coding genes")

    pc_only_gt = {"gene_type": ("protein_coding", operator.eq)}

    # print(filter_gtf_attributes(gtf, pc_only_gt)[["gene_id","gene_type"]])
    #
    # print("Trying membership filter, protein-coding & lncRNA genes")

    pc_lnc_gt = {"gene_type": (["protein_coding", "lncRNA"], "isin")}

    # filter_gtf_attributes(gtf, pc_lnc_gt)[["gene_id","gene_type"]].print(n=20)

    # print("Trying protein-coding or lncRNA, plus extract 'gene' rows only (i.e. multiple gene_type filters)")
    #
    # pc_lnc_gt["Feature"] = ("gene", operator.eq)
    #
    # filter_gtf_attributes(gtf, pc_lnc_gt)[["Feature", "gene_id","gene_type"]].print(n=20)

    # print("Trying protein-coding or lncRNA, transcript rows only plus extract protein-coding transcripts for protein-coding genes only")
    #
    # pc_only_tt = {"transcript_type": {"protein_coding": ("protein_coding", operator.eq)}}
    #
    # pc_lnc_gt["Feature"] = ("transcript", operator.eq)
    #
    # print("pc, lncRNA, transcript rows only")
    #
    # x = filter_gtf_attributes(gtf, pc_lnc_gt)[["Feature", "transcript_id", "gene_type", "transcript_type"]].print(n=20, chain=True)
    #
    # print(x.as_df()[["gene_type", "transcript_type"]].groupby("gene_type").describe())
    #
    #
    # print(f"N unique txipts\n{x.as_df()[['gene_type', 'transcript_id']].groupby('gene_type').nunique()}")
    #
    #
    # y = filter_gtf_attributes(gtf, pc_lnc_gt, pc_only_tt)[["Feature", "transcript_id", "gene_type", "transcript_type"]].print(n=20, chain=True)
    #
    # print(y.as_df()[["gene_type", "transcript_type"]].groupby("gene_type").describe())
    #
    # print(y[y.gene_type == "protein_coding"].transcript_type.value_counts())
    #
    # print(y[y.gene_type == "lncRNA"].transcript_type.value_counts())
    #
    # print(f"N unique txipts\n{y.as_df()[['gene_type', 'transcript_id']].groupby('gene_type').nunique()}")

    print("Trying protein-coding or lncRNA, transcript rows only plus extract protein-coding transcripts for protein-coding genes, lncRNA txipts for lncRNAs only")

    pc_only_tt = {"protein_coding": {"transcript_type": ("protein_coding", operator.eq)},
                  "lncRNA": {"transcript_type": ("lncRNA", operator.eq)
                             }
                  }

    pc_lnc_gt["Feature"] = ("transcript", operator.eq)

    print("pc, lncRNA, transcript rows only")

    x = filter_gtf_attributes(gtf, pc_lnc_gt)[["Feature", "transcript_id", "gene_type", "transcript_type"]].print(n=20, chain=True)

    y = filter_gtf_attributes(gtf, pc_lnc_gt, pc_only_tt)[["Feature", "transcript_id", "gene_type", "transcript_type"]].print(n=20, chain=True)

    print(x[x.gene_type == "protein_coding"].transcript_type.value_counts())

    print(x[x.gene_type == "lncRNA"].transcript_type.value_counts())


    print(f"N unique txipts\n{x.as_df()[['gene_type', 'transcript_id']].groupby('gene_type').nunique()}")

    print("-------\n AFTER EXTRA FILTERING \n------")

    print(y[y.gene_type == "protein_coding"].transcript_type.value_counts())

    print(y[y.gene_type == "lncRNA"].transcript_type.value_counts())

    print(f"N unique txipts\n{y.as_df()[['gene_type', 'transcript_id']].groupby('gene_type').nunique()}")


    print("\n TRY ABOVE PLUS TSL FILTERING FOR PC TRANSCRIPTS\n")

    pc_only_tt["protein_coding"]["transcript_support_level"] = (3, operator.le)

    gtf.transcript_support_level = gtf.transcript_support_level.replace("NA", np.nan).astype(float)

    print("This is transcript type filtering \n")
    print(pc_only_tt)

    z = filter_gtf_attributes(gtf, pc_lnc_gt, pc_only_tt)[["Feature", "transcript_id", "gene_type", "transcript_type", "transcript_support_level"]].print(n=20, chain=True)

    print(z[z.gene_type == "protein_coding"].transcript_type.value_counts())
    print(z[z.gene_type == "protein_coding"].transcript_support_level.value_counts())

    print(z[z.gene_type == "lncRNA"].transcript_type.value_counts())
