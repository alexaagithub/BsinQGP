****************************
          README
****************************

area setup:  
source /cvmfs/cms.cern.ch/cmsset_default.sh  
export SCRAM_ARCH=slc6_amd64_gcc700  
cmsrel CMSSW_10_3_1_patch3  
cd CMSSW_10_3_1_patch3/src  
mkdir UserCode  
git clone <repo>  
mkdir prefiltered_trees  

to compile: scram b  

to run: bmesons 

It runs over data and mc trees in lstore. 
To run over your trees, change inputs  

First delete the folders:
old/; results_eff/; results_eff_8_bins/; results_final/; results_new


and then do:  
mkdir results

mkdir results/Bu/mc_validation_plots/ss_mc

mkdir results/Bs/mc_validation_plots/ss_mc

mkdir results/Bu/mc_validation_plots/ss_sp

mkdir results/Bs/mc_validation_plots/ss_sp

mkdir results/Bu/mc_validation_plots/mc_sp

mkdir results/Bs/mc_validation_plots/mc_sp

mkdir results/Bu/mc_validation_plots/ss_mc_sp

mkdir results/Bs/mc_validation_plots/ss_mc_sp


mkdir results/Bu/Bpt

mkdir results/Bs/Bpt

mkdir results/Bu/sideband_sub/

mkdir results/Bs/sideband_sub/

mkdir results/Bu/splot/Bmass

mkdir results/Bs/splot/Bmass

mkdir results/Bu/splot/sig

mkdir results/Bs/splot/sig

mkdir results/Bu/splot/bkg

mkdir results/Bs/splot/bkg

mkdir results/Bu/splot/sig_bkg

mkdir results/Bs/splot/sig_bkg


mkdir results/Bu/pulls

mkdir results/Bs/pulls
