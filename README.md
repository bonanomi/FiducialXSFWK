# FiducialXS

Framework for fiducial and differential cross section measurements using CJLST TTrees for Run 2 data.

Before using this package setting up Combine:
```
export SCRAM_ARCH=slc7_amd64_gcc700
cmsrel CMSSW_10_2_13
cd CMSSW_10_2_13/src
cmsenv
git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
cd HiggsAnalysis/CombinedLimit
cd $CMSSW_BASE/src/HiggsAnalysis/CombinedLimit
git fetch origin
git checkout v8.1.0
scramv1 b clean; scramv1 b
```

## Brief presentation of the codes
In this section a brief description of the codes is given, together with the ideal workflow to run the analysis
1. **skim_MC_tree.cp** and **skim_data_tree.cpp**: Starting from CJLST TTrees, the branches we are interested in are selected only, both for data and MC. 
