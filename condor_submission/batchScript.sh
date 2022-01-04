#!/bin/bash

## ----- Setting pre-compiled CMSSW version with combine ----- ##
cp /home/llr/cms/tarabini/CMSSW_10_2_13_withCombine.tgz .
source /cvmfs/cms.cern.ch/cmsset_default.sh
tar -xf CMSSW_10_2_13_withCombine.tgz
rm CMSSW_10_2_13_withCombine.tgz
cd CMSSW_10_2_13/src/
scramv1 b ProjectRename
eval `scramv1 runtime -sh`

## ----- Setting kerberos to access EOS area ----- ##
export KRB5CCNAME=FILE:/etc/krb5c/$(whoami)
/opt/exp_soft/cms/t3/eos-login -username atarabin

mkdir coefficients coefficients/matrix_eff coefficients/matrix_eff/2016 coefficients/matrix_eff/2017 coefficients/matrix_eff/2018 coefficients/matrix_nonfid
mkdir combine_files
mkdir datacard datacard/datacard_2016 datacard/datacard_2017 datacard/datacard_2018
mkdir fit
mkdir inputs
mkdir LHScans LHScans/plots
mkdir templates
mkdir templates/plots templates/plots/2016 templates/plots/2017 templates/plots/2018
mkdir templates/plots/2016/OBS templates/plots/2017/OBS templates/plots/2018/OBS
mkdir impacts

cp /home/llr/cms/tarabini/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXSFWK/inputs/observables.py inputs/.
cp /home/llr/cms/tarabini/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXSFWK/inputs/higgs_xsbr_13TeV.py inputs/.

cd coefficients
cp /home/llr/cms/tarabini/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXSFWK/coefficients/RunCoefficients.py .
python RunCoefficients.py --obsName 'OBS' --obsBins 'BIN' --year 'Full'

cp /home/llr/cms/tarabini/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXSFWK/coefficients/pdfUncertainty.py .
python pdfUncertainty.py --obsName 'OBS' --obsBins 'BIN' --year 'Full'
# python pdfUncertainty.py --obsName 'OBS' --obsBins 'BIN' --year 'Full' --nnlops

cd ../templates
cp /home/llr/cms/tarabini/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXSFWK/templates/RunTemplates.py .
cp /home/llr/cms/tarabini/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXSFWK/templates/plot_templates.cpp .
cp /home/llr/cms/tarabini/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXSFWK/templates/setTDRStyle.C .
python RunTemplates.py --obsName 'OBS' --obsBins 'BIN' --year 'Full'
c++ -o  plot_templates plot_templates.cpp `root-config --cflags --glibs`
IN="BIN"
boundaries=$(echo $IN | tr "|" "\n")
./plot_templates OBS $boundaries

cd ../fit
cp /home/llr/cms/tarabini/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXSFWK/fit/RunFiducialXS.py .
cp /home/llr/cms/tarabini/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXSFWK/fit/createDatacard.py .
cp /home/llr/cms/tarabini/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXSFWK/fit/createXSworkspace.py .
cp /home/llr/cms/tarabini/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXSFWK/fit/addConstrainedModel.py .
cp /home/llr/cms/tarabini/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXSFWK/fit/impacts.py .
python RunFiducialXS.py --obsName 'OBS' --obsBins 'BIN' --year 'Full' UNBLIND
python impacts.py --obsName 'OBS' --obsBins 'BIN' --year 'Full' UNBLIND

cd ../LHScans
cp /home/llr/cms/tarabini/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXSFWK/LHScans/plot_LLScan.py .
cp /home/llr/cms/tarabini/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXSFWK/LHScans/plotting.py .
python plot_LLScan.py --obsName 'OBS' --obsBins 'BIN' --year 'Full' UNBLIND

cd ..
cp /home/llr/cms/tarabini/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXSFWK/tdrStyle.py .
cp /home/llr/cms/tarabini/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXSFWK/producePlots.py .
cp /home/llr/cms/tarabini/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXSFWK/inputs/accUnc_OBS_NNLOPS.py .
python producePlots.py --obsName 'OBS' --obsBins 'BIN' --year 'Full'

## ----- Moving all outputs ----- ##
# Outputs from RunCoefficients
mv inputs/inputs_sig_* /home/llr/cms/tarabini/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXSFWK/inputs
# Outputs from pdfUncertainty
mv inputs/accUnc_* /home/llr/cms/tarabini/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXSFWK/inputs
# Outputs from RunTemplates
mv inputs/inputs_bkg* /home/llr/cms/tarabini/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXSFWK/inputs
if [ -d /home/llr/cms/tarabini/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXSFWK/templates/2016/OBS ]; then
  rm -rf /home/llr/cms/tarabini/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXSFWK/templates/2016/OBS
  mkdir /home/llr/cms/tarabini/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXSFWK/templates/2016/OBS
else
  mkdir /home/llr/cms/tarabini/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXSFWK/templates/2016/OBS
fi
mv templates/2016/OBS/* /home/llr/cms/tarabini/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXSFWK/templates/2016/OBS/.
if [ -d /home/llr/cms/tarabini/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXSFWK/templates/2017/OBS ]; then
  rm -rf /home/llr/cms/tarabini/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXSFWK/templates/2017/OBS
  mkdir /home/llr/cms/tarabini/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXSFWK/templates/2017/OBS
else
  mkdir /home/llr/cms/tarabini/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXSFWK/templates/2017/OBS
fi
mv templates/2017/OBS /home/llr/cms/tarabini/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXSFWK/templates/2017/OBS/.
if [ -d /home/llr/cms/tarabini/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXSFWK/templates/2018/OBS ]; then
  rm -rf /home/llr/cms/tarabini/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXSFWK/templates/2018/OBS
  mkdir /home/llr/cms/tarabini/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXSFWK/templates/2018/OBS
else
  mkdir /home/llr/cms/tarabini/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXSFWK/templates/2018/OBS
fi
mv templates/2018/OBS /home/llr/cms/tarabini/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXSFWK/templates/2018/OBS/.
# Outputs from plot_templates
if [ -d /home/llr/cms/tarabini/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXSFWK/templates/plots/2016/OBS ]; then
  rm -rf /home/llr/cms/tarabini/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXSFWK/templates/plots/2016/OBS
  mkdir /home/llr/cms/tarabini/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXSFWK/templates/plots/2016/OBS
else
  mkdir /home/llr/cms/tarabini/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXSFWK/templates/plots/2016/OBS
fi
mv templates/plots/2016/OBS/* /home/llr/cms/tarabini/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXSFWK/templates/plots/2016/OBS/.
if [ -d /home/llr/cms/tarabini/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXSFWK/templates/plots/2017/OBS ]; then
  rm -rf /home/llr/cms/tarabini/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXSFWK/templates/plots/2017/OBS
  mkdir /home/llr/cms/tarabini/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXSFWK/templates/plots/2017/OBS
else
  mkdir /home/llr/cms/tarabini/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXSFWK/templates/plots/2017/OBS
fi
mv templates/plots/2017/OBS/* /home/llr/cms/tarabini/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXSFWK/templates/plots/2017/OBS/.
if [ -d /home/llr/cms/tarabini/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXSFWK/templates/plots/2018/OBS ]; then
  rm -rf /home/llr/cms/tarabini/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXSFWK/templates/plots/2018/OBS
  mkdir /home/llr/cms/tarabini/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXSFWK/templates/plots/2018/OBS
else
  mkdir /home/llr/cms/tarabini/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXSFWK/templates/plots/2018/OBS
fi
mv templates/plots/2018/OBS/* /home/llr/cms/tarabini/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXSFWK/templates/plots/2018/OBS/.
# Outputs from RunFiducialXS
mv fit/commands_OBS.py /home/llr/cms/tarabini/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXSFWK/fit/.
mv datacard/hzz4l_all_13TeV_xs_* /home/llr/cms/tarabini/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXSFWK/datacard/.
mv datacard/datacard_2016/* /home/llr/cms/tarabini/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXSFWK/datacard/datacard_2016/.
mv datacard/datacard_2017/* /home/llr/cms/tarabini/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXSFWK/datacard/datacard_2017/.
mv datacard/datacard_2018/* /home/llr/cms/tarabini/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXSFWK/datacard/datacard_2018/.
# Outputs from plot_LLScan.py
mv LHScans/plots/lhscan_compare* /home/llr/cms/tarabini/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXSFWK/LHScans/plots/.
mv LHScans/resultsXS_LHScan_* /home/llr/cms/tarabini/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXSFWK/LHScans/.
# Outputs from producePlots.py
if [ -d /home/llr/cms/tarabini/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXSFWK/plots/OBS ]; then
  rm -rf /home/llr/cms/tarabini/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXSFWK/plots/OBS
  mkdir /home/llr/cms/tarabini/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXSFWK/plots/OBS
else
  mkdir /home/llr/cms/tarabini/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXSFWK/plots/OBS
fi
mv plots/OBS/* /home/llr/cms/tarabini/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXSFWK/plots/OBS/.
# Outputs from impacts.py
mv impacts/* /home/llr/cms/tarabini/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXSFWK/impacts/.

mv combine_files/* /home/llr/cms/tarabini/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXSFWK/combine_files/.

rm -rf *
