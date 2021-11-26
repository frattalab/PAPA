import pyranges as pr
import pandas as pd
import sys
from functools import reduce
import operator

## Pass operators as 2nd argument for 'equality' subsets
# i.e. valid args are strings 'isin' or '~isin', or an operator (e.g. operator.eq())
# https://docs.python.org/3/library/operator.html


def _subset_membership(gr, col_name, filter_tuple):
    '''
    '''

    assert isinstance(filter_tuple[0], list) or isinstance(filter_tuple[0], set)
    assert filter_tuple[1] in ["isin", "~isin"]

    if filter_tuple[1] == "isin":
        return gr.subset(lambda df: df[col_name].isin(filter_tuple[0]))

    elif filter_tuple[1] == "~isin":
        return gr.subset(lambda df: ~df[col_name].isin(filter_tuple[0]))



def _subset_equality(gr, col_name, filter_tuple):
    '''
    '''

    return gr.subset(lambda df: filter_tuple[1](df[col_name], filter_tuple[0]))


def _subset_gr(gr, col_name, filter_tuple):
    '''
    '''

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
    Filter PyRanges with multi-level, group-specific filtering schemes given column/attribute values

    gr: PyRanges object

    gene_types: dict of {'column_name': (value, return_true)}
        - value - can be a string, numeric (int or float) or list type
        - return_true - specifies how elements in 'column_name' should relate to value in order for them to be returned
            - return_true can be an operator (from operator base module) or
            one of ['isin', '~isin'] (to test for membership)

    tr_types: optional dict of {'column_name': {'gene_type_value_list[0]': (value, return_true)}} for value-specific filtering for values in gene_types {col: (['a','b'], 'isin')}
        e.g. if my gene_types dict had {col: (['a','b'], 'isin')}
        You can specify additional filters for group 'a' (and/or 'b') in tr_types
        using same structure as above

    e.g. I want to extract 'protein_coding' & 'lncRNA' gene_types ('gene_type' GTF attribute key)
    for 'protein_coding' gene types, I want to extract TSL 1 transcripts ('transcript_support_level' attribute key)
    filter_gtf_attributes(gr,
                          gene_types={'gene_type': (['protein_coding', lncRNA], 'isin')
                                      },
                          tr_types={'transcript_support_level': {'protein_coding': (1, '=='),
                                                                 }
                                    }
                          )

    Notes:
    - not every key in gene_types needs to have additional filters defined
    - tr_types only becomes active if at least 1 of tuple[0] (value) in gene_types.values() is a list.
        - If you have multiple filters but no 'group-specific' filters, put all of your filters in gene_types dict


    '''
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

        # Which values from gene_types filters (values, filter_true) are undergoing additional filters?
        tr_type_keys = [tr_key for tr_dict in tr_types.values() for tr_key in tr_dict.keys()]
        # print(tr_type_keys)

        gr_list = []

        # Every column level filter in tr_types
        for tt_col, tt_dict in tr_types.items():
            #
            for tt_key, tt_tuple in tt_dict.items():
                # Check each key: val in gene_types_l -
                # tt_key should be in one of these gene_types_l.values() tuples
                found = False
                for gt_col, gt_tuple in gene_types_l.items():
                    if tt_key in gt_tuple[0]:

                        tt_gr = _subset_gr_and(gr2,
                                               # Name of col containing nested dict key
                                               col_name_1=gt_col,
                                               # Nested dict key
                                               filter_tuple_1=(tt_key, operator.eq),
                                               col_name_2=tt_col,
                                               filter_tuple_2=tt_tuple)

                        # print(f"this is subset gr for {tt_key}\n{tt_col}\n{tt_tuple}")
                        # print(tt_gr[[tt_col, gt_col]])
                        # print("this is the col of interest")
                        # print(tt_gr.as_df()[tt_col].value_counts())

                        gr_list.append(tt_gr)

                        found = True

                        # Need to check if all values in that column are also in tt_dict.keys()
                        # Otherwise will be filtered out at this step even though don't want a second filter
                        # This just keeps all rows satisfying the top level filteriing criteria

                        no_filt_list = [_subset_gr(gr2,
                                                   gt_col,
                                                   (ele, operator.eq)
                                                   )
                                        for ele in gt_tuple[0]
                                        if ele not in tr_type_keys]
                        #
                        # print(f"n without extra filter - {len(no_filt_list)}")
                        # print(no_filt_list)

                        gr_list.extend(no_filt_list)

                    else:
                        # Try next column: (filter tuple) in gene_types
                        continue

                # Have checked every filter defined in gene_types, this key is not present
                if not found:
                    raise Exception(f"{tt_key} must be present in the values (1st element of tuple) of gene_types dict")

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

    pc_only_tt = {"transcript_type": {"protein_coding": ("protein_coding", operator.eq),
                                      "lncRNA": ("lncRNA", operator.eq)
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

    
