# Work-Flow of the EIC-SIDIS Study

## Generate SIDIS Events with four setting:
 * D2 target, free PDF, by running ./GetSIDIS eic_d2_pion_free.dat N
 * D2 target, nPDF which should be the same as the free one (a consistancy check), by running ./GetSIDIS eic_d2_pion.dat N
 * C12 target, free PDF, by running ./GetSIDIS eic_c12_pion_free.dat N
 * C12 target, nPDF from EPS09, by running ./GetSIDIS eic_c12_pion.dat N

 Note that to accumulate enough statistics, we need to run many files (N up to a large number). 
 I have a farm script "gofarm" to submit jobs to batch-farms, also "run.sh" is usefull. 

## Bin the data and calculate observables, such as counts and average kinematic quantities in each bin
* Using "skim.C" to merge all ROOT files generated (0-N), and then split them into individual files for different bins
* Use this script, "MakeBinHist.C", to read the binned data from "skim" ROOT files for D2, C12 in two different configuration

## Read in binned data of D2 and C12, calculate asymmetry (or super-ratio)
 * A python Jupyter notebook script was written to do this job, named "asym_lo.ipynb"
 
