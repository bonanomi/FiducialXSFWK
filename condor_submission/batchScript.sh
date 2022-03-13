#!/bin/bash

## ----- Setting pre-compiled CMSSW version with combine ----- ##
cp /home/llr/cms/tarabini/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXSFWK/condor_submission/CMSSW_10_2_13_withCombine.tgz .
source /cvmfs/cms.cern.ch/cmsset_default.sh
tar -xf CMSSW_10_2_13_withCombine.tgz
rm CMSSW_10_2_13_withCombine.tgz
cd CMSSW_10_2_13/src/
scramv1 b ProjectRename
eval `scramv1 runtime -sh`

## ----- Setting kerberos to access EOS area ----- ##
# export KRB5CCNAME=FILE:/etc/krb5c/$(whoami)
# /opt/exp_soft/cms/t3/eos-login -username atarabin

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
mkdir helperstuff

cp /home/llr/cms/tarabini/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXSFWK/helperstuff/binning.py helperstuff/.
cp /home/llr/cms/tarabini/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXSFWK/helperstuff/observables.py helperstuff/.

cp /home/llr/cms/tarabini/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXSFWK/inputs/higgs_xsbr_13TeV.py inputs/.

cd coefficients
cp /home/llr/cms/tarabini/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXSFWK/coefficients/RunCoefficients.py .
python RunCoefficients.py --obsName 'VAR' --year 'Full' ZZFLOATING
FIRST
SECOND

jes=SETTING
if [ $jes == true ];then
  mkdir JES JES/tables JES/tables/OBS
  cd JES
  # cp /home/llr/cms/tarabini/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXSFWK/coefficients/JES/RunJES.py .
  cp /home/llr/cms/tarabini/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXSFWK/coefficients/JES/PrintJES.py .
  # cp /home/llr/cms/tarabini/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXSFWK/coefficients/JES/zx.py .
  cp /home/llr/cms/tarabini/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXSFWK/coefficients/JES/JESNP_OBS.py .
  cp /home/llr/cms/tarabini/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXSFWK/coefficients/JES/JESNP_evts_OBS.py .
  # python RunJES.py --obsName 'VAR' --year 'Full'
  python PrintJES.py --obsName 'VAR' --year '2016'
  python PrintJES.py --obsName 'VAR' --year '2017'
  python PrintJES.py --obsName 'VAR' --year '2018'
  cd ..
fi

cp /home/llr/cms/tarabini/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXSFWK/coefficients/pdfUncertainty.py .
python pdfUncertainty.py --obsName 'VAR' --year 'Full' ZZFLOATING
# python pdfUncertainty.py --obsName 'VAR' --year 'Full' --nnlops
cd ../inputs
sed "s/ggH125/ggH125_NNLOPS/g" accUnc_OBS.py > accUnc_OBS_NNLOPS.py #FIXME: This is temporary!

cd ../templates
cp /home/llr/cms/tarabini/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXSFWK/templates/RunTemplates.py .
cp /home/llr/cms/tarabini/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXSFWK/templates/plot_templates.py .
cp /home/llr/cms/tarabini/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXSFWK/templates/tdrStyle.py .
python RunTemplates.py --obsName 'VAR' --year 'Full' ZZFLOATING
# c++ -o  plot_templates plot_templates.cpp `root-config --cflags --glibs`
# IN="BIN"
# boundaries=$(echo $IN | tr "|" "\n")
# ./plot_templates OBS $boundaries
python plot_templates.py --obsName 'VAR' --year 'Full'


cd ../fit
cp /home/llr/cms/tarabini/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXSFWK/fit/RunFiducialXS.py .
cp /home/llr/cms/tarabini/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXSFWK/fit/createDatacard.py .
cp /home/llr/cms/tarabini/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXSFWK/fit/createXSworkspace.py .
cp /home/llr/cms/tarabini/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXSFWK/fit/addConstrainedModel.py .
cp /home/llr/cms/tarabini/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXSFWK/fit/impacts.py .
cp /home/llr/cms/tarabini/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXSFWK/fit/impacts_klambda.sh .
python RunFiducialXS.py --obsName 'VAR' --year 'Full' UNBLIND ZZFLOATING
if [[ "VAR" == "pT4l_kL" ]]; then
  sh impacts_klambda.sh
else
  python impacts.py --obsName 'VAR' --year 'Full' UNBLIND
  FIFTH UNBLIND
fi

cd ../LHScans
cp /home/llr/cms/tarabini/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXSFWK/LHScans/plot_LLScan.py .
cp /home/llr/cms/tarabini/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXSFWK/LHScans/plotting.py .
python plot_LLScan.py --obsName 'VAR' --year 'Full' UNBLIND
FOURTH UNBLIND

cd ..
if [[ "VAR" != "pT4l_kL" ]]; then
  cp /home/llr/cms/tarabini/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXSFWK/tdrStyle.py .
  cp /home/llr/cms/tarabini/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXSFWK/producePlots.py .
  cp /home/llr/cms/tarabini/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXSFWK/producePlots_v4.py .
  python producePlots.py --obsName 'VAR' --year 'Full' UNBLIND
  THIRD UNBLIND
fi

## ----- Moving all outputs ----- ##
# Outputs from RunCoefficients
mv inputs/inputs_sig_* /home/llr/cms/tarabini/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXSFWK/inputs

# Outputs from JES
if [ $jes == true ]; then

  # mv inputs/JES* /home/llr/cms/tarabini/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXSFWK/inputs

  # if [ -d /home/llr/cms/tarabini/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXSFWK/coefficients/JES/plots/OBS ]; then
  #   rm -r /home/llr/cms/tarabini/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXSFWK/coefficients/JES/plots/OBS
  #   mkdir /home/llr/cms/tarabini/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXSFWK/coefficients/JES/plots/OBS
  # else
  #   mkdir /home/llr/cms/tarabini/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXSFWK/coefficients/JES/plots/OBS
  # fi
  # mv coefficients/JES/plots/OBS/* /home/llr/cms/tarabini/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXSFWK/coefficients/JES/plots/OBS/.

  if [ -d /home/llr/cms/tarabini/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXSFWK/coefficients/JES/tables/OBS ]; then
    rm -r /home/llr/cms/tarabini/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXSFWK/coefficients/JES/tables/OBS
    mkdir /home/llr/cms/tarabini/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXSFWK/coefficients/JES/tables/OBS
  else
    mkdir /home/llr/cms/tarabini/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXSFWK/coefficients/JES/tables/OBS
  fi
  mv coefficients/JES/tables/OBS/* /home/llr/cms/tarabini/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXSFWK/coefficients/JES/tables/OBS/.

fi

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
mv templates/2017/OBS/* /home/llr/cms/tarabini/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXSFWK/templates/2017/OBS/.
if [ -d /home/llr/cms/tarabini/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXSFWK/templates/2018/OBS ]; then
  rm -rf /home/llr/cms/tarabini/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXSFWK/templates/2018/OBS
  mkdir /home/llr/cms/tarabini/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXSFWK/templates/2018/OBS
else
  mkdir /home/llr/cms/tarabini/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXSFWK/templates/2018/OBS
fi
mv templates/2018/OBS/* /home/llr/cms/tarabini/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXSFWK/templates/2018/OBS/.

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
if [ -d /home/llr/cms/tarabini/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXSFWK/impacts/OBS ]; then
  rm -rf /home/llr/cms/tarabini/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXSFWK/impacts/OBS
  mkdir /home/llr/cms/tarabini/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXSFWK/impacts/OBS
else
  mkdir /home/llr/cms/tarabini/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXSFWK/impacts/OBS
fi
mv impacts/* /home/llr/cms/tarabini/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXSFWK/impacts/OBS/.

# Combine files
mv combine_files/* /home/llr/cms/tarabini/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXSFWK/combine_files/.

rm -rf *
