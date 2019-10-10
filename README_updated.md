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

First run selection_skimmer.cc!  

It runs over data and mc samples in lstore  
right now for B+  
to run over other samples, change inputs  
It loops over the data or mc samples, flattening the tree, applying the selection, and writing a new tree with less variables  

running example: selection_skimmer --colsyst PbPb --mc 0 --channel 1  

Then you can run sideban_subtraction.cc over the outputs of selection_skimmer that are in prefiltered_trees  
It applies a sideband subtraction method to the data and then compares the result with the monte carlo distributions  

First do:  
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


To run:  
bmesons 
