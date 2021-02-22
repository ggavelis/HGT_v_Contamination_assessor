# HGT_v_Contamination_assessor

This script takes existing metadata about each DNA/AA sequence, and uses that--in combination with an alien index value--to determine whether each sequence should be flagged as a contaminant. This script expects 2 inputs (1) a fasta file to decontaminate, (2) a 'supertsv' metadata file that contains the following fields  (seq_id | alien_index_value | num_splice_variants | lineage_of_best_BLAST_hit | spliced_leader[True/False] | polyA_tail[True/False]). This script also expects 1 parameter -- the AI cutoff used for screening. (default AI_cutoff = 0.01)

#### Rationale

Alien indices can be used as heuristics to infer whether a sequence is likely to be native or foreign (e.g. a contaminant or HGT). But decontaminating a 'dirty' dataset based on AI alone is inadvisable, since this approach is also likely to remove bona fide HGT. To mitigate this problem of 'overcleaning,' I have broken AI cleaning into two steps.
I. A first-pass "flagging" step that flags alls seqs whose AI excede the cutoff.
II. A second-pass "rescue" step that uses sequence metadata to redeem certain sequences. For example:
  A. Any sequence with a dinoflagellate spliced-leader is unflagged as native. (Inferred from upstream script that builds the supertsv)
  B. Any best-hit to prokaryotes is unflagged if it has:
    i. A poly-A tail (Inferred from upstream script)
    ii. Multiple splice isoforms (inferred from Trinity, via upstream script)

This script also gathers metrics about the frequency of HGT from various groups.
