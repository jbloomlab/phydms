# Data for the model adequacy profiling script/test

This directory contains the data to run a profiling test for a model adequacy test with `phydms`.

## Original data
This data is from the divergence timing manuscript

* [`HA_hybridDoud_prefs.csv`](HA_hybridDoud_prefs.csv): Mike's WSN preferences with only the shared sites between H1 and H3.  
* [`HA_hybrid_H1.fasta`](HA_hybrid_H1.fasta): 34 H1 sequences.  

## Short data
The data below were created using the[`_shorten.py`](_shorten.py) script.
They are the preferences and alignment from above but shortened to the number of sites and the number of sequences specified in the script.  

* [`HA_short_prefs.csv`](HA_short_prefs.csv)
* [`HA_short.fasta`](HA_short.fasta)
* [`HA_short_tree.newick`](HA_short_tree.newick)

## `phydms` outputs

The script [`run_phydms.bas`](run_phydms.bash) runs `phydms` with the short data listed above.
The outputs are found in the directory [`phydms/`](phydms/)

## Scripts

I make the smaller test set using the script [`_shorten.py`](_shorten.py).
This script takes in a preference set, an alignment, the target number of sites, and the target number of sequences.
The outputs include a preference set, an alignment, and a tree built by `RAxML`.


I run `phydms` on the small data using the script [`phydms/`](phydms/).

## Notes

I tried to make a preference set with just *one site* but `RAxML` threw an error.
I tried to make a preference set with just *two sites* but `phdyms` threw an error.
