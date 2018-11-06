"""
Make a short test case

SKH 20180619
"""

import pandas as pd
import subprocess
from Bio import SeqIO
import glob
import os

def main():
    total_sites = 10
    total_number_seqs = 5
    prefs = "HA_hybridDoud_prefs.csv"
    alignment = "HA_hybrid_H1.fasta"

    # make new preference set
    prefs = pd.read_csv(prefs)
    prefs = prefs.head(total_sites)
    assert prefs['site'].iloc[-1] == len(prefs)
    prefs.to_csv("HA_short_prefs_nsites{0}.csv".format(total_sites), index=False)

    # make new alignment
    final = []
    for seq in SeqIO.parse(alignment, "fasta"):
        seq.seq = seq.seq[:total_sites*3]
        assert len(seq.seq) == (3 * total_sites)
        final.append(seq)
    final = final[:total_number_seqs]
    with open("HA_short_nsites{0}_nseqs{1}.fasta".format(total_sites, total_number_seqs), "w") as output_handle:
        SeqIO.write(final, output_handle, "fasta")

    # make new tree
    for fname in glob.glob("RAxML*"):
        os.remove(fname)
    raxml_cmd = ["raxml", "-s", "HA_short_nsites{0}_nseqs{1}.fasta".format(total_sites, total_number_seqs), "-n", "temp", "-m", "GTRCAT", "-p1", "-T", "2"]
    subprocess.check_call(raxml_cmd)
    os.rename("RAxML_bestTree.temp", "HA_short_tree_nsites{0}_nseqs{1}.newick".format(total_sites, total_number_seqs))
    for fname in glob.glob("RAxML*"):
        os.remove(fname)


if __name__ == '__main__':
    main()
