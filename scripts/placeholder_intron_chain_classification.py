import pyranges as pr
import numpy as np
import pandas as pd
import os
import sys
import time

# A plus strand transcript
test_ref_tr1 = {"Chromosome": [1]*4,
                "Start": [10,100,200,300],
                "End": [30,120,220,340],
                "Strand": ["+"]*4,
                "Feature": ["exon"]*4,
                "gene_id": ["ref_gene_1"]*4,
                "transcript_id": ["ref_tr_1"] * 4}


pr.from_dict(test_ref_tr1)

# {"Chromosome": [],
#                 "Start": [],
#                 "End": [],
#                 "Strand": [],
#                 "Feature": [],
#                 "gene_id": [],
#                 "transcript_id": []
#                }


# A minus strand transcript
test_ref_tr2 = {"Chromosome": [2]*3,
                "Start": [10,80,100],
                "End": [20,90,120],
                "Strand": ["-"]*3,
                "Feature": ["exon"]*3,
                "gene_id": ["ref_gene_2"]*3,
                "transcript_id": ["ref_tr_2"]*3
               }


test_ref = pr.concat([pr.from_dict(test_ref_tr1), pr.from_dict(test_ref_tr2)])
test_ref

# Now make test novel transcripts to cover my test cases

# 1. Novel last exons in first intron of annotated transcript (e.g. STMN2)
# For these to pass, they should share an identical 3'end with a first exon of a known transcript

test_novel_fi = {"Chromosome": [1]*4,
                "Start": [10,50]*2,
                "End": [30,70,35,70],
                "Strand": ["+"]*4,
                "Feature": ["exon"]*4,
                "gene_id": ["nov_gene_1"]*4,
                "transcript_id": ["nov_tx_fi_p"]*2 + ["nov_tx_fi_f"]*2,
               }

pr.from_dict(test_novel_fi)

# 2. Internal intron, spliced in last exon (fully contained within last exon) (e.g. ONECUT1)
# For these to pass, they should match the intron chain of a known transcript up until the penultimate exon

test_novel_si = {"Chromosome": [1]*6,
                "Start": [10,100,140] + [50,100,140],
                "End": [30,120,160] + [70,120,160],
                "Strand": ["+"]*6,
                "Feature": ["exon"]*6,
                "gene_id": ["nov_gene_1"]*6,
                "transcript_id": ["nov_tx_si_p"]*3 + ["nov_tx_si_f"]*3
               }

pr.from_dict(test_novel_si)

# 3. Internal intron bleedthrough (e.g. SIN3B)
# For these events to pass, they should match the intron chain of a known transcript up until the penultimate exon

test_novel_bl = {"Chromosome": [1]*6,
                "Start": [10,100,200] + [50,100,200],
                "End": [30,120,240] + [70,120,240],
                "Strand": ["+"]*6,
                "Feature": ["exon"]*6,
                "gene_id": ["nov_gene_1"]*6,
                "transcript_id": ["nov_tx_bl_p"]*3 + ["nov_tx_bl_f"]*3
               }

pr.from_dict(test_novel_bl)

# 4. Internal intron with novel internal and terminal exon (e.g.)
# For these events to pass, they should match the intron chain of a known transcript,
# but have a continuous chain of length n of novel events at the 3'end of the transcript
# (n can be varied)
# Event know from NP is fully contained within annotated intron - also set this constraint?

test_novel_mult = {"Chromosome": [1]*4,
                   "Start": [10,100,130,150],
                   "End": [30,120,140,160],
                   "Strand": ["+"]*4,
                   "Feature": ["exon"]*4,
                   "gene_id": ["nov_gene_1"]*4,
                   "transcript_id": ["nov_tx_mult_p"]*4
                    }

pr.from_dict(test_novel_mult)

# 5. 3'UTR intron fully contained within an annotated 3'UTR (e.g. TDP-43)
# For this to pass, they should match the intron chain of a known transcript up until the penultimate exon
# (Annotate as a 3'UTR intron (spliced out) after filtering for intron chain match)
test_novel_3ui = {"Chromosome": [1]*5,
                "Start": [10,100,200,300,330],
                "End": [30,120,220,310,340],
                "Strand": ["+"]*5,
                "Feature": ["exon"]*5,
                "gene_id": ["nov_gene_1"]*5,
                "transcript_id": ["nov_tx_3ui_p"]*5
                 }

pr.from_dict(test_novel_3ui)

#6. Distal last exon spliced from penultimate exon (i.e. a mutually exclusive last exon) (e.g. SMC1A)
# For this to pass, they should match the intron chain of a known transcript up until the penultimate exon
# Will have to annotate more precisely later (i.e. differentiate from 7)

test_novel_exc_dist = {"Chromosome": [1]*4,
                "Start": [10,100,200,360],
                "End": [30,120,220,380],
                "Strand": ["+"]*4,
                "Feature": ["exon"]*4,
                "gene_id": ["nov_gene_1"]*4,
                "transcript_id": ["nov_tx_exc_dist_p"]*4
               }

pr.from_dict(test_novel_exc_dist)

#7. Distal last exon spliced from ann
# For this to pass, they should match the intron chain of a known transcript up until the penultimate exon
test_novel_le_dist = {"Chromosome": [1]*5,
                      "Start": [10,100,200,300,360],
                      "End": [30,120,220,310,380],
                      "Strand": ["+"]*5,
                      "Feature": ["exon"]*5,
                      "gene_id": ["nov_gene_1"]*5,
                      "transcript_id": ["nov_tx_le_dist_p"]*5
                     }

pr.from_dict(test_novel_le_dist)

#8. Distal last exon spliced from annoatrd (minus strand)
# For this to pass, they should match the intron chain of a known transcript up until the penultimate exon
test_novel_le_dist_minus = {"Chromosome": [2]*3,
                      "Start": [40,80,100],
                      "End": [50,90,120],
                      "Strand": ["-"]*3,
                      "Feature": ["exon"]*3,
                      "gene_id": ["nov_gene_2"]*3,
                      "transcript_id": ["nov_tx_le_dist_p_minus"]*3
                     }

pr.from_dict(test_novel_le_dist_minus)

test_novel_gr = pr.concat([pr.from_dict(event) for event in [test_novel_3ui,
                                                            test_novel_bl,
                                                            test_novel_exc_dist,
                                                            test_novel_fi,
                                                            test_novel_le_dist,
                                                            test_novel_mult,
                                                            test_novel_si,
                                                            test_novel_le_dist_minus]
                          ]
                         )

test_novel_gr


def introns_from_df(df):
    '''
    '''

    n_exons = len(df)

    if n_exons < 2:
        return None
        #print(df)
        #raise Exception("at least two exons are required for transcript to have an intron")
    # n exons = n-1 introns

    strand = df["Strand"].drop_duplicates().tolist()[0]
#     print(strand)
    chrom = df["Chromosome"].drop_duplicates().tolist()[0]
    gene_id = df["gene_id"].drop_duplicates().tolist()[0]
    tx_id = df["transcript_id"].drop_duplicates().tolist()[0]
    feature = "intron"
    introns = {}
    for i in range(0, n_exons - 1):
        if strand == "+":
            intron_start = df.iloc[i, lambda x: x.columns.get_loc("End")]
            intron_end = df.iloc[i+1, lambda x: x.columns.get_loc("Start")]
            introns[str(i)] = {"Chromosome": chrom,
                               "Start": intron_start,
                               "End": intron_end,
                               "Strand": strand,
                               "Feature": feature,
                               "gene_id": gene_id,
                               "transcript_id": tx_id}
        elif strand == "-":
            intron_start = df.iloc[i, lambda x: x.columns.get_loc("End")]
            intron_end = df.iloc[i+1, lambda x: x.columns.get_loc("Start")]
            introns[str(i)] = {"Chromosome": chrom,
                               "Start": intron_start,
                               "End": intron_end,
                               "Strand": strand,
                               "Feature": feature,
                               "gene_id": gene_id,
                               "transcript_id": tx_id}
    return pd.DataFrame.from_dict(introns, orient = "index")



def introns_by_tx(gr, by="transcript_id", nb_cpu=1):
    '''
    '''
    # Sort by position (for safety)
    gr = gr.sort()

    return gr.apply(lambda df: df.groupby(by).apply(introns_from_df), nb_cpu=nb_cpu)


test_novel_introns = introns_by_tx(test_novel_gr)
test_ref_introns = introns_by_tx(test_ref)



def rle(inarray):
        """
        run length encoding. Partial credit to R rle function.
        Multi datatype arrays catered for including non Numpy
        returns: tuple (runlengths, startpositions, values)
        https://stackoverflow.com/questions/1066758/find-length-of-sequences-of-identical-values-in-a-numpy-array-run-length-encodi
        Thomas Browne
        """
        ia = np.asarray(inarray)                # force numpy
        n = len(ia)
        if n == 0:
            return (None, None, None)
        else:
            y = ia[1:] != ia[:-1]               # pairwise unequal (string safe)
            i = np.append(np.where(y), n - 1)   # must include last element posi
            z = np.diff(np.append(-1, i))       # run lengths
            p = np.cumsum(np.append(0, z))[:-1] # positions
            return(z, p, ia[i])

def match_intron_chains(novel_gr, ref_gr, id_col = "transcript_id", nb_cpu = 1):
    '''
    grs should contain introns
    '''
    # {novel_id: {matches: [ref_id], chain_match: [n], terminal_non_match: [n]}}
    match_info_dict = {}



    for key, dfs in pr.itergrs([novel_gr, ref_gr], strand=True, keys=True):
        print("----processing chrom and strand pair {0} & {1}".format(key[0], key[1]))
        # dfs = novel_gr & ref_gr matched by key (chromosome & strand)
        # pandas group by objects of (transcript_id, df)
        by_tx = tuple([df.groupby("transcript_id") for df in dfs])

        novel_txipts = by_tx[0]
        ref_txipts = by_tx[1]

        # Pyranges keys are tuples of (chr,strand)
        strand = key[1]


        #Comparing each novel transcript against all ref transcripts
        for novel_id, novel_introns in novel_txipts:

            if strand == "-":
                # Standard PyRanges sort - First in df = last intron (smallest values i.e. leftmost)
                # Reverse so first row in df is always the first intron
                # is reset_index necessary?
                novel_introns = novel_introns[::-1].reset_index(drop=True)

            else:
                pass

            # As a first pass, check for matches in intron chain of first introns between novel and reference Txs
            if strand == "+":
                first_matches = [np.array_equal(novel_introns.head(1)[["Start","End"]],
                                                ref_introns.head(1)[["Start","End"]])
                                 for ref_id,ref_introns in ref_txipts]
            else:
                first_matches = [np.array_equal(novel_introns.head(1)[["Start","End"]],
                                                ref_introns[::-1].reset_index(drop=True) # rev order so 1st row = 1st intron
                                                .head(1)[["Start","End"]]
                                               )
                                 for ref_id,ref_introns in ref_txipts]


            if not sum(first_matches) > 0:
                print("{0} does not match any reference transcripts in its first intron. Skipping".format(novel_id))
                continue

            # Compare full introns chains of novel transcript against each ref transcript with match in first intron.

            for ref_tr, first_match in zip(ref_txipts, first_matches):

                if not first_match:
                    continue


                ref_id = ref_tr[0]
#                 print("ref_id {}".format(ref_id))
                ref_introns = ref_tr[1]

                if strand == "-":
                    ref_introns = ref_introns[::-1].reset_index(drop=True) # reverse so first row = first intron

                # To avoid a slicing error for ref txipts shorter than novel
                n_novel_introns = len(novel_introns)
                n_ref_introns = len(ref_introns)
                max_chain = min(n_novel_introns, n_ref_introns)
#             print(max_chain)
#             print(ref_tr[1].iloc[0,])
#             print(n_novel_introns)
#             print(len(ref_tr[1]))

                # Row-wise, check whether match with corresponding intron of reference transcript
                novel_chain_match = pd.DataFrame([(novel_introns.iloc[i,:][["Start","End"]]
                                                   .eq(ref_introns.iloc[i,:][["Start","End"]]
                                                      )
                                                  )
                                                  for i in range(max_chain)])

                # Collapse to single True/False per row - does intron completely match?
                novel_chain_match = novel_chain_match.apply(np.all, axis=1, raw=True)

#             print(novel_chain_match)

                runs, starts, vals = rle(novel_chain_match)
            #print(runs[0])
            #print(starts)
            #print("\n {0}".format(np.where(vals == False)))

                # Considered a valid match if a intron chain completely identical
                # or matches at beginning of txipt but differs at the end
                if np.all(vals) or np.array_equal(vals, [True,False]):

                    # Don't want to throw away all valid matches (yet)
                    # Possible genuine match, should update dict with info
                    # match (from start of ref Txipt)

                    if vals.size == 1:
                        #i.e. all true/introns match (e.g. bleedthrough)
                        terminal_non_match = 0
                    else:
                        # All valid = consective match & non-match, so non-match = 2nd in array
                        terminal_non_match = runs[1]

#                 print(terminal_non_match)

                    if novel_id not in match_info_dict:

                        match_info_dict[novel_id] = {"matches": [ref_id],
                                                     "chain_match": [runs[0]], #Always starts with true, so take length of true
                                                     "terminal_non_match": [terminal_non_match],
                                                    }

                    else:
                        # Append to dict
                        match_info_dict[novel_id]["matches"].append(ref_id)
                        match_info_dict[novel_id]["chain_match"].append(runs[0])
                        match_info_dict[novel_id]["terminal_non_match"].append(terminal_non_match)


                else:
                    continue

    # Output df for easier parsing
#     match_df = pd.concat({key: pd.DataFrame.from_dict(d, orient = "columns") for key, d in match_info_dict.items()}, axis=0)

    return match_info_dict


path_stie_gtf_chr1_nov = "../two_sample_example_output/stringtie/chr1.no_ref_id.TDP43-F_S6.assembled.gtf"
path_ref_gtf = "../data/annotation/gencode.v34.annotation.gtf"

stie_chr1 = pr.read_gtf(path_stie_gtf_chr1_nov)
stie_chr1

ref = pr.read_gtf(path_ref_gtf)
ref

stie_chr1_exons = stie_chr1.subset(lambda df: df["Feature"] == "exon")
stie_chr1_exons

ref_exons = ref.subset(lambda df: df["Feature"] == "exon")
ref_chr1_exons = ref_exons.subset(lambda df: df["Chromosome"] == "chr1")
ref_exons


ref_chr1_pc_exons = ref_chr1_exons.subset(lambda df: df.gene_type.isin(["protein_coding", "lncRNA"]))
ref_chr1_pc_exons


stie_chr1_introns = introns_by_tx(stie_chr1_exons,nb_cpu=1)
stie_chr1_introns

ref_chr1_introns = introns_by_tx(ref_chr1_pc_exons, nb_cpu=2)
ref_chr1_introns


### Alternative (hopefully more scalable approach)

# Try and use pyranges internals as much as possible, and avoid manual looping and comparisons

# Essentially, Find overlapping introns, then filter for those that are identical


# New idea


def intron_id(gr):

    return gr.assign("intron_id",
                                         lambda df: pd.Series([":".join([tx_id, str(start), str(end)])
                                                     for tx_id, start, end in zip(df["transcript_id"],
                                                                                  df["Start"],
                                                                                  df["End"])
                                                    ])
                                        )

def sort_introns_by_strand(df):
    '''
    '''
    # first reset_index call removes the original index of the group (e.g. row 4005 in df)
    # second reset_index call adds the sorted index as a column to the dataframe (the order along exon in each transcript)
    if (df.Strand == '+').all():
        return df.sort_values(by=['End']).reset_index(drop=True).reset_index()
    elif (df.Strand == '-').all():
        return df.sort_values(by=['Start'], ascending=False).reset_index(drop=True).reset_index()


def filter_transcripts_by_chain(novel_exons,ref_exons, match_type = "transcript",nb_cpu = 1):
    '''
    '''

    novel_cols_to_keep = ["Feature","transcript_id"]
    ref_cols_to_keep = ["Feature", "transcript_id", "gene_id", "gene_name"]

    assert match_type in ["transcript", "any"], "match_type must be one of 'transcript' or 'any'. value passed - {}".format(str(match_type))

    #1. Find introns by transcript & give each intron a unique ID
    novel_introns = introns_by_tx(novel_exons).sort()
    ref_introns = introns_by_tx(ref_exons).sort()

    novel_introns = intron_id(novel_introns)
    ref_introns = intron_id(ref_introns)


    #2. Track number of introns in each novel transcript
    novel_tx_intron_counts = (novel_introns.as_df()
                              .groupby("transcript_id").size())


    # novel_introns, ref_introns

    # 3. Store intron_ids for each transcript, sorted by intron_number (where 1 = first intron regardless of strand) in a df/Series
    novel_intron_ids_ordered = (novel_introns.as_df()
                                .groupby("transcript_id")
                                .apply(sort_introns_by_strand)
                                .reset_index(drop=True)
                                .rename({"index": "intron_number"}, axis="columns")
                               )
    novel_intron_ids_ordered["intron_number"] = novel_intron_ids_ordered["intron_number"].add(1)

    # df of txipt_id | intron_id | intron_number
    novel_intron_ids_ordered = novel_intron_ids_ordered.loc[:,["transcript_id","intron_id","intron_number"]]
#     print(novel_intron_ids_ordered.dtypes)


    #4. Find novel introns with any overlap with reference introns
    # Inner join to add ref_rows to novel gr
    joined = novel_introns.join(ref_introns, strandedness="same",how=None, suffix ="_ref",nb_cpu=nb_cpu)

    #5. Filter for overlaps that exactly match (or differ by given tolerance)
    joined = joined.subset(lambda df: abs(df.Start - df.Start_ref) + abs(df.End - df.End_ref) <= 0, nb_cpu=nb_cpu)

    # Minimal info needed on matches between novel and reference introns
    joined = joined.as_df()[["transcript_id","intron_id","transcript_id_ref","intron_id_ref"]]

#     print(joined.dtypes)

    #6. Join ordered novel introns with match info

    if match_type == "transcript":
        # Looking for intron to match any annotated intron, regardless of reference transcript
        novel_ref_match_info = novel_intron_ids_ordered.merge(joined,
                                                              how="left",
                                                              on="intron_id",
                                                              suffixes=["_novel","_match"]
                                                             )

    else:
        # Looking for introns (except last) to match the same reference transcript
        # merge_ordered can do a 'grouped merge' filling in empty rows (introns) for each transcript_id
        # This is especially useful if want to do transcript-specific intron matching
        # For each reference transcript, all novel introns will be filled with NaN if no overlap for given transcript_id
        # (i.e. novel txipt matches all but last intron of reference transcript)

        novel_ref_match_info = (pd.merge_ordered(novel_intron_ids_ordered,
                                    joined,
                                    how="left",
                                    on="intron_id",
                                    right_by="transcript_id_ref", # group matches by ref tx & join tx by tx
                                    suffixes=["_novel","_match"],
                                    fill_method=None)
                   .sort_values(by=["transcript_id_novel","intron_number"])
                   .dropna(axis="columns", subset="transcript_id_match") # This fills for
              )


    return novel_ref_match_info
