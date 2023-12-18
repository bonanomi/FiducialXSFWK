#!/bin/bash

#Check if the name of the observable is provided in the command line
if [[ $# -eq 0 ]] ; then
    echo 'no name of the observable is provided'
    exit 1
fi

# path="/eos/user/a/atarabin/www/fiducial/run2_cardsValidation" #ReReco
path="/eos/user/a/atarabin/www/fiducial/run2_UL_cardsValidation" #UL

/opt/exp_soft/cms/t3/eos-login -init -username atarabin

cd coefficients
source /opt/exp_soft/llr/root/vv6.20.06-el7-gcc9-py37/etc/init.sh
python RunPlotCoefficients.py --obsName "${1}" --year 'Full'
cd ../fit
python RunPlotCorrelation.py --obsName "${1}" --year 'Full'
cd ..

if [ -d "$path/${1}" ]; then
  ### Take action if $DIR exists ###
  rm -r $path/${1}
  mkdir $path/${1}
  cp $path/index.php $path/${1}
else
  ###  Control will jump here if $DIR does NOT exists ###
  mkdir $path/${1}
  cp $path/index.php $path/${1}
fi

mkdir $path/${1}/datacard_2016 $path/${1}/datacard_2017 $path/${1}/datacard_2018
cp $path/index.php $path/${1}
cp $path/index.php $path/${1}/datacard_2016
cp $path/index.php $path/${1}/datacard_2017
cp $path/index.php $path/${1}/datacard_2018

mkdir $path/${1}/templatesBkgs_2016 $path/${1}/templatesBkgs_2017 $path/${1}/templatesBkgs_2018
cp $path/index.php $path/${1}/templatesBkgs_2016/.
cp $path/index.php $path/${1}/templatesBkgs_2017/.
cp $path/index.php $path/${1}/templatesBkgs_2018/.
cp templates/plots/2016/${1}/* $path/${1}/templatesBkgs_2016/.
cp templates/plots/2017/${1}/* $path/${1}/templatesBkgs_2017/.
cp templates/plots/2018/${1}/* $path/${1}/templatesBkgs_2018/.

cp coefficients/matrix_eff/2016/eff_2016_${1}_* $path/${1}/.
cp coefficients/matrix_eff/2017/eff_2017_${1}_* $path/${1}/.
cp coefficients/matrix_eff/2018/eff_2018_${1}_* $path/${1}/.

cp coefficients/matrix_nonfid/2016/nonFid_2016_${1}_* $path/${1}/.
cp coefficients/matrix_nonfid/2017/nonFid_2017_${1}_* $path/${1}/.
cp coefficients/matrix_nonfid/2018/nonFid_2018_${1}_* $path/${1}/.

# Move datacards
if [ ${1} = "njets_pt30_eta4p7" ]; then
  obs_datacards="NJ"
elif [ ${1} = "pT4l" ]; then
  obs_datacards="PTH"
elif [ ${1} = "rapidity4l" ]; then
  obs_datacards="YH"
elif [ ${1} = "pTj1" ]; then
  obs_datacards="PTJET"
else
  obs_datacards=${1}
fi

cp datacard/datacard_2016/hzz4l_*_13TeV_xs_${obs_datacards}_bin*_* $path/${1}/datacard_2016/.
cp datacard/datacard_2017/hzz4l_*_13TeV_xs_${obs_datacards}_bin*_* $path/${1}/datacard_2017/.
cp datacard/datacard_2018/hzz4l_*_13TeV_xs_${obs_datacards}_bin*_* $path/${1}/datacard_2018/.
# Move ws
if [ ${1} = "pT4l_kL" ]; then
  cp datacard/datacard_2016/hzz4l_*_13TeV_xs_SM_125_${obs_datacards}_kLambda* $path/${1}/datacard_2016/.
  cp datacard/datacard_2017/hzz4l_*_13TeV_xs_SM_125_${obs_datacards}_kLambda* $path/${1}/datacard_2017/.
  cp datacard/datacard_2018/hzz4l_*_13TeV_xs_SM_125_${obs_datacards}_kLambda* $path/${1}/datacard_2018/.
  cp datacard/hzz4l_all_13TeV_xs_${1}_bin_kLambda* $path/${1}/.
else
  cp datacard/datacard_2016/hzz4l_*_13TeV_xs_SM_125_${obs_datacards}_v* $path/${1}/datacard_2016/.
  cp datacard/datacard_2017/hzz4l_*_13TeV_xs_SM_125_${obs_datacards}_v* $path/${1}/datacard_2017/.
  cp datacard/datacard_2018/hzz4l_*_13TeV_xs_SM_125_${obs_datacards}_v* $path/${1}/datacard_2018/.
  cp datacard/hzz4l_all_13TeV_xs_${1}_bin_v* $path/${1}/.
fi

cp plots/${1}/asimov/${1}_unfoldwith* $path/${1}/.
cp plots/${1}/asimov/corr_${1}_*.png $path/${1}/.
cp plots/${1}/asimov/model/* $path/${1}/.
if [ ${1} = "pT4l_kL" ]; then
  cp impacts/${1}/impacts_v3_kappa_lambda* $path/${1}/.
  cp LHScans/plots/lhscan_compare_${1}_kappa* $path/${1}/.
else
  cp impacts/${1}/impacts_*_${1}_*_asimov* $path/${1}/.
  cp LHScans/plots/lhscan_compare_${1}_r* $path/${1}/.
fi

cp fit/commands_${1}.py $path/${1}/.
cp impacts/${1}/commands_impacts_${1}_v* $path/${1}/.
