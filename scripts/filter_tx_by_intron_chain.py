#!/usr/bin/env python3

import pyranges as pr
import numpy as np
import pandas as pd
import argparse
import os
import sys
import logging
from timeit import default_timer as timer


'''
Filter assembled transcripts for those that match reference transcripts in their intron chain up until their penultimate intron
This allows identification of novel last exons (conservatively), but limits addition of false positive transcripts
'''


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


def rle(inarray):
        """
        run length encoding. Partial credit to R rle function.
        Multi datatype arrays catered for including non Numpy
        returns: tuple (runlengths, startpositions, values)
        https://stackoverflow.com/questions/1066758/find-length-of-sequences-of-identical-values-in-a-numpy-array-run-length-encoding
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



def validate_matching_chain(df, max_terminal_non_match=2):
    '''
    apply to grouped df
    '''

    runs, starts, vals = rle(df["match"])

    # Valid matches are:
    # All introns match (e.g. bleedthrough event)
    # All but last x introns match (usually one) (i.e. runs = 1,0)

    if np.all(vals):
        # All introns match (e.g. bleedthrough event)
        return True

    elif np.array_equal(vals, [1,0]) and runs[-1] <= max_terminal_non_match:
        # all but last x introns match (x = max_terminal_non_match) (i.e. runs = 1,0)
        return True

    else:
        return False


def filter_transcripts_by_chain(novel_exons, novel_introns, ref_exons, ref_introns, match_type = "transcript", max_terminal_non_match=2, nb_cpu = 1):
    '''
    '''

    novel_cols_to_keep = ["Feature","transcript_id"]
    ref_cols_to_keep = ["Feature", "transcript_id", "gene_id", "gene_name"]

    assert match_type in ["transcript", "any"], "match_type must be one of 'transcript' or 'any'. value passed - {}".format(str(match_type))

    #1. Find introns by transcript & give each intron a unique ID
    # print("finding introns...")
    # t1 = timer()

    #novel_introns = introns_by_tx(novel_exons, nb_cpu=nb_cpu).sort()
    #ref_introns = introns_by_tx(ref_exons, nb_cpu=nb_cpu).sort()

    # t2 = timer()

    # print("took {} (s)".format(t2 - t1))
    print("filtering transcripts by intron chain matching...")

    print("sorting all grs by position for safety...")
    t1 = timer()

    novel_exons = novel_exons.sort()
    novel_introns = novel_introns.sort()
    ref_exons = ref_exons.sort()
    ref_introns = ref_introns.sort()

    t2 = timer()
    print("took {} (s)".format(t2 - t1))

    print("adding intron_id column...")

    t3 = timer()
    novel_introns = intron_id(novel_introns)
    ref_introns = intron_id(ref_introns)
    t4 = timer()

    print("took {} s".format(t4 - t3))

    #2. Track number of introns in each novel transcript
    # novel_tx_intron_counts = (novel_introns.as_df()
                              # .groupby("transcript_id").size())


    # novel_introns, ref_introns

    # 3. Store intron_ids for each transcript, sorted by intron_number (where 1 = first intron regardless of strand) in a df/Series
    print("generating df of novel txipts sorted by intron number...")

    t5 = timer()

    novel_intron_ids_ordered = (novel_introns.as_df()
                                .groupby("transcript_id")
                                .apply(sort_introns_by_strand)
                                .reset_index(drop=True)
                                .rename({"index": "intron_number"}, axis="columns")
                               )
    novel_intron_ids_ordered["intron_number"] = novel_intron_ids_ordered["intron_number"].add(1)

    # df of txipt_id | intron_id | intron_number
    novel_intron_ids_ordered = novel_intron_ids_ordered.loc[:,["transcript_id","intron_id","intron_number"]]

    t6 = timer()
    print("took {} s".format(t6 - t5))
#     print(novel_intron_ids_ordered.dtypes)


    #4. Find novel introns with any overlap with reference introns
    # Inner join to add ref_rows to novel gr
    print("finding overlaps between novel and reference introns...")

    t7 = timer()

    # Have to convert Starts and Ends to np.int64 to prevent left-join error
    # https://github.com/biocore-ntnu/pyranges/issues/170
    # TO DO: re-report
    joined = pr.PyRanges(novel_introns.as_df(), int64=True).join(pr.PyRanges(ref_introns.as_df(), int64=True), strandedness="same", suffix ="_ref", nb_cpu=nb_cpu)

    t8 = timer()

    print("took {} s".format(t8 - t7))

    #5. Filter for overlaps that exactly match (or differ by given tolerance - for now not supported)
    print("filtering overlaps for exact matches...")

    t9 = timer()
    joined = joined.subset(lambda df: abs(df.Start - df.Start_ref) + abs(df.End - df.End_ref) <= 0, nb_cpu=nb_cpu)
    t10 = timer()

    print("took {} s".format(t10 - t9))

    # Minimal info needed on matches between novel and reference introns
    joined = joined.as_df()[["transcript_id","intron_id","transcript_id_ref","intron_id_ref"]]

#     print(joined.dtypes)

    #6. Join ordered novel introns with match info
    #7. Assign a simple tracker column 'match' of True (where intron is matched) and False (where intron is not matched)

    print("preparing for filtering intron matches...")
    t11 = timer()

    if match_type == "any":
        # Looking for intron to match any annotated intron, regardless of reference transcript
        novel_ref_match_info = novel_intron_ids_ordered.merge(joined,
                                                              how="left",
                                                              on="intron_id",
                                                              suffixes=["_novel","_match"]
                                                             )

        # Assign 'match' column for each intron.
        # Since we don't really care which intron it matches, & no matches will mean NaN
        novel_ref_match_info["match"] = novel_ref_match_info["transcript_id_ref"]
        novel_ref_match_info["match"] = novel_ref_match_info["match"].fillna(0)
        novel_ref_match_info["match"] = novel_ref_match_info["match"].replace("\w*", 1, regex=True)

        novel_ref_match_info = novel_ref_match_info.drop_duplicates(subset=["intron_id"])

        # Minimal informative info is novel tx, novel intron_id & number, match column
        novel_ref_match_info = novel_ref_match_info[["transcript_id_novel","intron_id","intron_number","match"]]


    elif match_type == "transcript":
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
                               )

        # merge_ordered fills rows for each intron for each ref tx in df, regardless of whether any overlap
        # .dropna(axis="rows", subset=["intron_id_ref"])
        novel_ref_match_info = (novel_ref_match_info.groupby(["transcript_id_novel", "transcript_id_ref"])
                                .filter(lambda df: (df["intron_id_ref"].notna()).any()) # Retained if ref tx has >=1 matching introns
                                .reset_index(drop=True))

        # Make a match column where 1 = match, 0 = no match for each ref id and novel intron
        novel_ref_match_info["match"] = novel_ref_match_info["intron_id_ref"]
        novel_ref_match_info["match"] = novel_ref_match_info["match"].fillna(0)
        novel_ref_match_info["match"] = novel_ref_match_info["match"].replace("\w*", 1, regex=True)

    t12 = timer()
    print("took {} s".format(t12 - t11))


    # 8. Filter down matching transcripts to those that all ref introns except penultimate or all introns...
    print("filtering for valid intron chain matches...")
    t13 = timer()
    if match_type == "any":
        # Only need to check by novel transcript_id
        filt_novel_ref_match_info = (novel_ref_match_info.groupby("transcript_id_novel")
                                     .filter(lambda x: validate_matching_chain(x, max_terminal_non_match)
                                            )
                                    )

    elif match_type == "transcript":
        # Check novel tx vs each ref tx
        filt_novel_ref_match_info = (novel_ref_match_info.groupby(["transcript_id_novel","transcript_id_ref"])
                                     .filter(lambda x: validate_matching_chain(x, max_terminal_non_match)
                                            )
                                    )
    t14 = timer()
    print("took {} s".format(t14 - t13))


    # Return simplified df of novel transcript_id & matching transcript_id if applicable

    if match_type == "any":
        return filt_novel_ref_match_info["transcript_id_novel"].drop_duplicates()

    elif match_type == "transcript":
        return filt_novel_ref_match_info[["transcript_id_novel","transcript_id_ref"]].drop_duplicates()


def main(novel_path, ref_path, match_by, max_terminal_non_match, out_gtf, nb_cpu):
    '''
    '''
    start = timer()

    print("reading in input gtf files, this can take a while...")
    print("reading gtf containing novel assembled transcripts...")

    s1 = timer()
    novel = pr.read_gtf(novel_path)
    e1 = timer()

    print("took {} s".format(e1 - s1))

    print("reading gtf containing reference transcripts...")

    s2 = timer()
    ref = pr.read_gtf(ref_path)
    e2 = timer()

    print("took {} s".format(e2 - s2))

    print("extracting protein-coding & lncRNA gene types from reference annotation file...")

    start2 = timer()
    ref_pc = ref.subset(lambda df: df["gene_type"].isin(["protein_coding", "lncRNA"]), nb_cpu=nb_cpu)
    end2 = timer()

    print("took {} s".format(end2 - start2))

    print("extracting novel exons from input GTF file...")

    start4 = timer()
    novel_exons = novel.subset(lambda df: df["Feature"] == "exon", nb_cpu=nb_cpu)
    end4 = timer()

    print("took {} s".format(end4 - start4))

    print("extracting exons of protein_coding and lncRNA genes from reference annotation...")

    start5 = timer()
    ref_pc_exons = ref.subset(lambda df: df["Feature"] == "exon", nb_cpu=nb_cpu)
    end5 = timer()

    print("took {} s".format(end5 - start5))

    print("finding introns for each reference transcript...")

    start3 = timer()

    try:
        ref_pc_introns = ref_pc.features.introns(by="transcript", nb_cpu=nb_cpu)
    except KeyError:
        # Specific error with ray when execute introns func for the 2nd time in same script...
        # KeyError: 'by'
        # Whilst avoiding working out what's going on, I can run my super slow intron finding script...
        print("pr.features.introns returned KeyError ('by'), using my hacky intron finding workaround...")
        ref_pc_introns = introns_by_tx(ref_pc_exons, nb_cpu=nb_cpu)

    end3 = timer()

    print("took {} s".format(end3 - start3))

    print("finding introns for each novel transcript...")

    start1 = timer()

    try:
        novel_introns = novel.features.introns(by="transcript", nb_cpu=nb_cpu)
    except KeyError:
        # Specific error with ray when execute introns func for the 2nd time in same script...
        # KeyError: 'by'
        # Whilst avoiding working out what's going on, I can run my super slow intron finding script...
        print("pr.features.introns returned KeyError ('by'), using my hacky intron finding workaround...")
        novel_introns = introns_by_tx(novel_exons, nb_cpu=nb_cpu)

    end1 = timer()

    print("took {} s".format(end1 - start1))

    print("finding novel transcripts with valid matches in their intron chain to reference transcripts...")

    start6 = timer()
    valid_matches = filter_transcripts_by_chain(novel_exons,
                                                novel_introns,
                                                ref_pc_exons,
                                                ref_pc_introns,
                                                match_type=match_by,
                                                max_terminal_non_match=max_terminal_non_match,
                                                nb_cpu=nb_cpu)
    end6 = timer()

    print("took {} s".format(end6 - start6))

    if isinstance(valid_matches, pd.Series):
        # match type was any
        valid_novel = novel.subset(lambda df: df["transcript_id"].isin(set(valid_matches.tolist())), nb_cpu=nb_cpu)
        valid_novel.to_gtf(out_gtf)

    elif isinstance(valid_matches, pd.DataFrame):
        # match_by/match_type was transcript
        valid_novel = novel.subset(lambda df: df["transcript_id"].isin(set(valid_matches["transcript_id_novel"].tolist())), nb_cpu=nb_cpu)
        valid_novel.to_gtf(out_gtf)


    end = timer()

    print("Completed: script took {} s / {} min (3dp) ".format(round(end - start, 3), round((end - start) / 60, 3)))





if __name__ == '__main__':

    descrpn = """Script to filter assembled transcripts for intron chain matches with reference transcripts
                 to identify novel 3'end transcript isoforms"""

    parser = argparse.ArgumentParser(description=descrpn)

    parser.add_argument("-i", "--input-transcripts", default='', dest="novel_gtf", help = "path to GTF file containing novel transcripts assembled by StringTie.", required=True)
    parser.add_argument("-r", "--reference-transcripts", default='', type=str, dest="ref_gtf", help="path to GTF file containing reference transcripts against which to match intron chains of novel transcripts. Should contain same chromosome naming scheme", required=True)
    parser.add_argument("-m", "--match-by", default="any", type=str, choices=["any", "transcript"], dest="match_by", help="Consider novel transcript a valid match if all but penultimate intron(s) match introns of any transcript ('any') or the same transcript ('transcript') (default: %(default)s)")
    parser.add_argument("-n", "--max-terminal-non-match", default=1, type=int, dest="max_terminal_non_match", help="Maximum number of uninterrupted reference-unmatched introns at 3'end of novel transcript for it to be considered a valid match (default: %(default)s)")
    parser.add_argument("-c", "--cores", default=1, type=int, help="number of cpus/threads for parallel processing (default: %(default)s)")
    parser.add_argument("-o", "--output-gtf", type=str, default="intron_chain_matched_transcripts.gtf", dest="output_gtf", help="name of output GTF file containing novel transcripts with intron chain matches (default: %(default)s)")

    if len(sys.argv) == 1:
        parser.print_help()
        parser.exit()

    args = parser.parse_args()



    main(args.novel_gtf, args.ref_gtf, args.match_by, args.max_terminal_non_match, args.output_gtf, args.cores)
