#!/usr/bin/env python3

import pyranges as pr
import pandas as pd
import pyfaidx
from Bio.Seq import Seq
from Bio import motifs
from papa_helpers import eprint, get_terminal_regions, add_region_number
import argparse
import os



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
