# FiducialXS

Framework for fiducial and differential cross section measurements using CJLST TTrees for Run 2 data.

Before using this package setting up Combine:
'''
export SCRAM\_ARCH=slc7\_amd64\_gcc700
cmsrel CMSSW\_10\_2\_13
cd CMSSW\_10\_2\_13/src
cmsenv
git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
cd HiggsAnalysis/CombinedLimit
cd $CMSSW\_BASE/src/HiggsAnalysis/CombinedLimit
git fetch origin
git checkout v8.1.0
scramv1 b clean; scramv1 b
'''
