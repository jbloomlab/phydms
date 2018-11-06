# Data for testing model adequacy

This directory contains data to test the model adequacy protocol in `phydms`.

## Original data
This data is from the divergence timing manuscript

* [`HA_hybridDoud_prefs.csv`](HA_hybridDoud_prefs.csv): Mike's WSN preferences with only the shared sites between H1 and H3.  
* [`HA_hybrid_H1.fasta`](HA_hybrid_H1.fasta): 34 H1 sequences.  

## Short data
The data below were created using the[`shorten.py`](shorten.py) script.
They are the preferences and alignment from above but shortened to the number of sites and the number of sequences specified in the script.  

* [`HA_short_prefs_nsites10.csv`](HA_short_prefs_nsites10.csv)  
* [`HA_short_nsites10_nseqs34.fasta`](HA_short_nsites10_nseqs34.fasta)  
* [`HA_short_tree_nsites10_nseqs34.newick`](HA_short_tree_nsites10_nseqs34.newick)  
* [`HA_short_nsites10_nseqs5.fasta`](HA_short_nsites10_nseqs5.fasta)  
* [`HA_short_tree_nsites10_nseqs5.newick`](HA_short_tree_nsites10_nseqs5.newick)  

## `phydms` outputs

The script [`run_phydms.bas`](run_phydms.bash) runs `phydms` with the short data listed above.
The outputs are found in the directory `phydms_nsites*_nseqs*/`.

## Scripts

I make the smaller test set using the script [`shorten.py`](shorten.py).
This script takes in a preference set, an alignment, the target number of sites, and the target number of sequences.
The outputs include a preference set, an alignment, and a tree built by `RAxML`.


I run `phydms` on the small data using the script [`run_phydms.bash`](`run_phydms.bash`).

## Notes

I tried to make a preference set with just *one site* but `RAxML` threw an error.
I tried to make a preference set with just *two sites* but `phydms` threw an error.
