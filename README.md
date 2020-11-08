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
1. **skim_MC_tree.cpp** and **skim_data_tree.cpp**: Starting from CJLST TTrees, the branches we are interested in are selected only, both for data and signal MC
2. **templates** folder: Templates and normalization coefficients for the backgrounds' PDF are extracted from MC (ggZZ and qqZZ) and data (ZX)
3. **coefficients** folder: All the coefficients of the signal parameterization are calculated
4. **fit** folder: The maximum likelihood fit is performed
5. **LHScans**: Likelihood scans are plotted, best-fit values and the correspodnign uncertainties are calculated
6. **producePlots.py**: Unfolded differential distributions are plotted
