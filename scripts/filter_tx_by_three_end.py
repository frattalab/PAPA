#!/usr/bin/env python3

import pyranges as pr
import pandas as pd
import numpy as np
import pyfaidx
from Bio.Seq import Seq
from Bio import motifs
from papa_helpers import eprint, get_terminal_regions, add_region_number
import argparse
import os
import sys
from timeit import default_timer as timer



'''
As a second step of filtering to look for specificity of called last exon
Check whether:
1. Predicted 3'end of transcript falls within x nt of known 3'ends (e.g. PolyASite)

2. Last x nt of last exon contains any conserved polyadenylation signal motif
    - These motifs are typically positionally enriched ~20-25nt upstream of the cleavage site
    - StringTie is unlikely to get the precise position right (TECTool showed this (Gypas 2018), generally a challenge for all tools that use SR data)
    - Without considering positional specificity, PAS motifs are considerably enriched in 3'UTRs (~80-85 %)
    and are depleted in other regions of the transcripts (< 10 % in 5'UTRs and internal coding exons)
    (Sethi et al.; 2021 (F3UTER preprint)). Therefore, the presence of a PAS motif close to the predicted 3'end
    should provide evidence that a genuine cleavage event is happening in this region
    (even if StringTie doesn't get precise cleavage site right)
    - e.g. 1/18 positionally enriched motifs discovered by Gruber et al., 2016 (PolyASite v1.0))
    or 12 discovered by Beaudong 2000 (all of which captured by Gruber et al., 2016). All of these motifs are
    conserved in human and mouse.

Intron-chain filtered transcripts that pass either of these two filters are retained for further analysis

The script takes as input:
- A GTF file of intron-chain filtered transcripts from 'filter_tx_by_intron_chain.py'
- 'Match stats' TSV file output by filter_tx_by_intron_chain.py
    - Only include event types that involve a polyadenylation event (e.g. exclude 3'UTR introns - these are by default reported in output)
- BED file of reference poly(A) sites (e.g. from PolyASite atlas)
- PAS motifs
    - 'Gruber' (18 from Gruber 2016) and 'Beaudong' (12 from Beaudong 2000) are built into script.
    - A TXT file of (DNA) motifs, one per-line, can also be supplied
- Max distance to nearest polyA site (int)
- Upstream region length for PAS motif searching (int)

The script outputs:
- a filtered GTF containing passed transcripts
- a TSV of 'match information' e.g.
    - filter status,
    - PAS motifs found and their location relative to 3'end
'''

beaudong_pas_motifs = """AATAAA
ATTAAA
TATAAA
AGTAAA
AATACA
CATAAA
AATATA
GATAAA
AATGAA
AAGAAA
ACTAAA
AATAGA""".split("\n")

gruber_pas_motifs = """AATAAA
ATTAAA
TATAAA
AGTAAA
AATACA
CATAAA
AATATA
GATAAA
AATGAA
AAGAAA
ACTAAA
AATAGA
AATAAT
AACAAA
ATTACA
ATTATA
AACAAG
AATAAG""".split("\n")

def region_around_pas(gr,
                      extend_5p=100,
                      extend_3p=0):
    '''
    Return 3'end of each interval of gr extended upstream by extend_5p and downstream by extend_3p
    Note: pr.three_end() returns the final nucleotide of each region, which is a 1 length interval
    If you wish to have a region of specified length (e.g. Final 100nt), you should provide extend_5p
    as your target length -1
    '''

    three = gr.three_end()

    return three.slack({"5": extend_5p, "3": extend_3p})


def nearest_atlas_site(gr_3p,
                       atlas_gr,
                       max_distance,
                       class_outcol="atlas_filter",
                       distance_outcol="nearest_atlas_distance"):
    '''
    Return gr with distance (distance_outcol) to nearest region (poly(A) site) defined in atlas_gr
    and class_outcol specifying whether distance is greater than (0) max_distance or not (1)
    '''

    # Define cols from atlas_gr to drop at the end
    # pr.nearest returns cols from atlas_gr for nearest feature
    # in out gr only want to retain the nearest distance from this object
    # Cols unique to atlas_gr won't have suffix added, so need to ID here
    gr_3p_cols = gr_3p.columns.tolist()
    cols_to_drop = [col for col in atlas_gr.columns if col not in gr_3p_cols]

    gr_3p_nr = gr_3p.nearest(atlas_gr,
                          strandedness="same",
                          overlap=True,
                          how=None # either direction
                          ).drop(cols_to_drop).drop(like="_b$")

    # Rename distance column
    gr_3p_nr = gr_3p_nr.apply(lambda df: df.rename({"Distance": distance_outcol}, axis=1))



    # Add 'class_outcol' - 1/0 if nearest site is <= (1) or > (0) max_distance
    gr_3p_nr = gr_3p_nr.assign(class_outcol,
                            lambda df: pd.Series(np.where(df[distance_outcol] <= max_distance,
                                                 1,
                                                 0)
                                                 )
                            )

    return gr_3p_nr


def find_pas_motifs(gr,
                    pas_motifs,
                    seq_col="seq",
                    class_outcol="motif_filter",
                    motifs_outcol="pas_motifs",
                    ):
    '''
    '''

    assert seq_col in gr.columns
    assert isinstance(pas_motifs, list)

    # Check that objects stored in seq_col are BioPython 'Seq' objects
    # (gr.as_df()[seq_col].apply(lambda x: isinstance(x, Seq))).all()

    # Generate motifs object
    # from generate a list of Seq objects storing motifs
    pas_seqs = [Seq(motif) for motif in pas_motifs]
    pas_motifs = motifs.create(pas_seqs)


    # Search for exact match to any provided motif in each region
    out_gr = gr.assign(motifs_outcol,
                   lambda df: df.seq.apply(lambda region:
                                          ",".join([str(idx) + "_" + str(motif)
                                                    for idx, motif in pas_motifs.instances.search(region)
                                                    ]
                                                   )
                                           ),
                       )

    # If no match found then an empty string is output
    # Add col specifying any match (1) / no match (0)

    out_gr = out_gr.assign(class_outcol,
                           lambda df: pd.Series(np.where(df[motifs_outcol] != "",
                                                         1,
                                                         0)
                                               )
                          )

    return out_gr



def main(gtf_path,
         match_stats_path,
         fasta_path,
         atlas_path,
         pas_motifs,
         max_atlas_dist,
         motif_search_length,
         output_prefix):
    '''
    '''

    assert isinstance(pas_motifs, list)

    eprint("reading in GTF of intron-chain filtered transcripts...")
    gtf = pr.read_gtf(gtf_path)

    eprint("reading in TSV of match status for intron-chain filtered transcripts...")
    match_stats = pd.read_csv(match_stats_path, sep="\t")


    #1. Filter for valid polyadenylation events want to check for 3'end genuineness
    # Since spliced 3'UTR intron events do not have a cleavage event, will not assess here
    # These will be reported without filtering in the output GTF
    eprint("Extracting valid polyadenylation events from match stats TSV & gtf...")
    valid_utr_spliced = set(match_stats.loc[(match_stats["match_class"] == "valid") &
                                            (match_stats["isoform_class"] == "3utr_intron_spliced"),
                                            "transcript_id_novel"])


    valid_to_test = set(match_stats.loc[(match_stats["match_class"] == "valid") &
                                        (match_stats["isoform_class"] != "3utr_intron_spliced"),
                                        "transcript_id_novel"])

    # Length of valid_to_test set check (if 0 output empty...)

    gtf_to_test = gtf.subset(lambda df: df["transcript_id"].isin(valid_to_test))

    #2 Get last exons for each transcript and extract three ends

    le = get_terminal_regions(gtf_to_test.subset(lambda df: df.Feature == "exon"),
                              source="stringtie",
                              filter_single=True)

    le_polya = le.three_end()


    #3. Find distance from predicted 3'end to nearest atlas site
    # Reporting 2 columns - 1 with distance (nt), 1 specifying whether at/below max_distance (1/0)

    eprint("Reading in BED file of poly(A) sites...")
    atlas = pr.read_bed(atlas_path)

    eprint("Finding nearest atlas poly(A) sites to predicted 3'ends...")
    le_polya = nearest_atlas_site(le_polya, atlas, max_atlas_dist)


    #4. Search for PAS motifs within last x nt of each predicted transcript...

    #A - Extend upstream by user-specified nt
    # pr.three_end() returns the final nucleotide of each region (1nt length interval)
    # To get specified length (e.g. Final 100nt), need to subtract 1 from provided value

    le_polya = le_polya.slack({"5": motif_search_length - 1})

    #B - Read in sequence for last x nt from genome fasta
    eprint(f"Reading in sequence for last {str(motif_search_length)} nt of each input transcript from genome fasta...")

    # Have to unstrand as fasta files only key by chromosome
    three_end_seqs = pr.get_fasta(le_polya.unstrand(), fasta_path)

    # Convert to BioPython Seq objects
    three_end_seqs = three_end_seqs.apply(lambda x: Seq(x))

    # Add to gr
    le_polya.seq = three_end_seqs


    #C - Search for presence of any of provided motifs in last x nt of Txs
    # Returns col with binary found/not_found (1/0) & motif-seq + position
    eprint(f"Finding poly(A) signal motifs in last {motif_search_length} nt of each transcript...")

    motif_start = timer()

    le_polya = find_pas_motifs(le_polya, pas_motifs)

    motif_end = timer()

    eprint(f"Finding motifs complete: took {motif_end - motif_start} s")

    #5.Generate 'match_stats' dataframe plus summary counts dfs
    # Subset to df of Tx_id, atlas_filter, motif_filter, atlas_distance & motifs_found
    # Find set of valid transcript IDs that pass either filter

    eprint("Generating 'match class' table and 'summary stats' tables...")

    pas_match_stats = le_polya.as_df()[["transcript_id",
                                        "atlas_filter",
                                        "motif_filter",
                                        "nearest_atlas_distance",
                                        "pas_motifs"
                                        ]]

    # A - add 'match_class' column ('valid/not_valid') if either of atlas_filter/motif_filter are 1
    pas_match_stats["match_class"] = np.where((pas_match_stats["atlas_filter"] == 1) |
                                               (pas_match_stats["motif_filter"] == 1),
                                               "valid",
                                               "not_valid")

    pas_valid_ids = set(pas_match_stats.loc[pas_match_stats["match_class"] == "valid",
                                            "transcript_id"])

    # Create dummy spliced 3'UTR introns match_stats df
    # i.e. all 'valid', but NaN for other columns as can't assess

    utr_spliced_match_stats = pd.DataFrame({"transcript_id": list(valid_utr_spliced),
                                            "match_class": ["valid"]*len(valid_utr_spliced),
                                            "atlas_filter": [np.nan]*len(valid_utr_spliced),
                                            "motif_filter": [np.nan]*len(valid_utr_spliced),
                                            "nearest_atlas_distance": [np.nan]*len(valid_utr_spliced),
                                            "pas_motifs": [np.nan]*len(valid_utr_spliced)
                                            }
                                           )

    pas_match_stats = pd.concat([pas_match_stats, utr_spliced_match_stats],
                                ignore_index=True)

    pas_match_stats = pas_match_stats.rename({"transcript_id": "transcript_id_novel"},
                                             axis=1)

    # B - Return isoform_class to pas_match_stats
    pas_match_stats = pas_match_stats.merge(match_stats[["transcript_id_novel",
                                                         "isoform_class"]],
                                            on="transcript_id_novel",
                                            how="left")

    # C - Generate 'valid & 'not_valid' summary counts dfs
    valid_summary = (pas_match_stats.loc[lambda x: x["match_class"] == "valid", :]
                                    .drop_duplicates(subset=["transcript_id_novel"])
                                    ["isoform_class"]
                                    .value_counts(dropna=False)
                     )

    nv_summary = (pas_match_stats.loc[lambda x: x["match_class"] == "not_valid", :]
                                 .drop_duplicates(subset=["transcript_id_novel"])
                                 ["isoform_class"]
                                 .value_counts(dropna=False)
                  )

    eprint("Transcripts passing atlas/motif filters by isoform class...")
    eprint(valid_summary)

    eprint("Transcripts failing atlas/motif filters by isoform class...")
    eprint(nv_summary)

    #6. Filter input GTF for valid IDs and output GTF to file
    # get set of valid IDs from pas_match_stats (& filter GTF for these)
    eprint("Filtering input GTF for transcripts passing atlas/motif filters...")
    valid_ids = set(pas_match_stats.loc[pas_match_stats["match_class"] == "valid",
                                        "transcript_id_novel"]
                    )


    valid_gtf = gtf.subset(lambda df: df.transcript_id.isin(valid_ids))

    eprint(f"writing filtered GTF to {output_prefix + '.gtf'}")
    valid_gtf.to_gtf(output_prefix + ".gtf")

    eprint("writing 'match class' and 'summary counts' tables to TSV...")

    pd.DataFrame(valid_summary).reset_index().to_csv(output_prefix + ".valid.class_summary_counts.tsv",
                                                     sep="\t",
                                                     header=["isoform_class","count"],
                                                     index=False,
                                                     na_rep="NA")

    pd.DataFrame(nv_summary).reset_index().to_csv(output_prefix + ".not_valid.class_summary_counts.tsv",
                                                  sep="\t",
                                                  header=["isoform_class","count"],
                                                  index=False,
                                                  na_rep="NA")

    pas_match_stats.to_csv(output_prefix + ".match_stats.tsv",
                           sep="\t",
                           header=True,
                           index=False,
                           na_rep="NA")




if __name__ == '__main__':

    start = timer()

    descrpn = """Script to filter intron-chain filtered transcripts (filter_tx_by_intron_chain.py)
    for transcripts containing 'high confidence' predicted cleavage regions, based on proximity to PolyASite PAS
    or presence of poly(A) signal motifs close to predicted 3'end"""

    parser = argparse.ArgumentParser(description=descrpn,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter, # Add defaults to end of help strings
                                     )

    parser.add_argument("-i", "--input-gtf",
                        dest="input_gtf", type=str,
                        default=argparse.SUPPRESS,
                        required=True,
                        help="Path to GTF file containing intron-chain filtered transcripts output by filter_tx_by_intron_chain.py")

    parser.add_argument("-s","--match-status-tsv",
                        dest="match_status_tsv",
                        type=str,required=True,
                        default=argparse.SUPPRESS,
                        help="Path to '<prefix>.match_stats.tsv' file corresponding to input intron-chain filtered transcripts output by filter_tx_by_intron_chain.py")

    parser.add_argument("-f","--fasta",
                        type=str,
                        required=True,
                        default=argparse.SUPPRESS,
                        help="Path to genome FASTA file sequence. It's expected that a .fai file (FASTA index, generated by samtools faidx) is also present at the same location.")

    parser.add_argument("-a","--pas-atlas",
                        dest="atlas", type=str,
                        required=True,
                        default=argparse.SUPPRESS,
                        help="Path to BED file of poly(A) sites defined from orthogonal sequencing (e.g. 3'end sequencing, PolyASite atlas)")

    parser.add_argument("-p","--pas-motifs",
                        dest="motifs",type=str,
                        required=True,
                        default=argparse.SUPPRESS,
                        help="polyA signal motifs defined by Gruber 2016 (PolyASite v1.0, 'Gruber', 18 motifs)" +
                             "or Beaudong 2000 ('Beaudong', 12 motifs (all of which recaptured by Gruber)). " +
                             "Path to TXT file containing PAS motifs (DNA, one per line) can also be supplied")

    parser.add_argument("-m","--max-atlas-distance",
                        dest="max_atlas_dist", type=int,
                        default=100,
                        help="Maximum distance (nt) between predicted 3'end and nearest atlas polyA site for transcript to be retained")

    parser.add_argument("-u","--motif-upstream-region-length",
                        dest="motif_upstream_length", type=int,
                        default=100,
                        help="length (nt) of region from upstream to predicted 3'end in which to search for presence of any of defined poly(A) signal motifs (retained if true)")

    parser.add_argument("-o","--output-prefix",
                        dest="output_prefix",
                        type=str,
                        default="three_end_matched_transcripts",
                        help="Prefix for output files. '.<suffix>' added depending on output file type ('.gtf' for valid transcripts GTF, '.match_stats.tsv' for TSV with match/filtering status etc.)",)

    if len(sys.argv) == 1:
        parser.print_help()
        parser.exit()

    args = parser.parse_args()

    # Define choice of PAS motifs
    if args.motifs.capitalize() == "Gruber":
        pas_motifs = gruber_pas_motifs

    elif args.motifs.capitalize() == "Beaudong":
        pas_motifs = beaudong_pas_motifs

    elif os.path.exists(args.motifs):
        eprint(f"reading in pas motifs from file - {args.motifs}")
        with open(args.motifs) as infile:
            pas_motifs = [line.rstrip("\n") for line in infile]

        eprint(f"pas motifs {', '.join(pas_motifs)}")

    else:
        raise ValueError(f"-p/--pas-motifs argument invalid - must be one of 'Gruber', 'Beaudong' or a valid path to TXT file. You passed {args.motifs}")


    main(args.input_gtf,
         args.match_status_tsv,
         args.fasta,
         args.atlas,
         pas_motifs,
         args.max_atlas_dist,
         args.motif_upstream_length,
         args.output_prefix)

    end = timer()

    eprint(f"Complete: took {round(end - start, 3)} s / {round((end - start) / 60, 3)} min (3 dp)")
