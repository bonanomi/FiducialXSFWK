#!/bin/bash

#Check if the name of the observable is provided in the command line
if [[ $# -eq 0 ]] ; then
    echo 'no name of the observable is provided'
    exit 1
fi

path="/eos/user/a/atarabin/www/fiducial/run2_cardsValidation"

/opt/exp_soft/cms/t3/eos-login -init -username atarabin

cd coefficients
source /opt/exp_soft/llr/root/vv6.20.06-el7-gcc9-py37/etc/init.sh
python RunPlotCoefficients.py --obsName "$1" --year 'Full'
cd ..

if [ -d "$path/$1" ]; then
  ### Take action if $DIR exists ###
  rm -r $path/$1
  mkdir $path/$1
  cp $path/index.php $path/$1
else
  ###  Control will jump here if $DIR does NOT exists ###
  mkdir $path/$1
  cp $path/index.php $path/$1
fi

mkdir $path/$1/datacard_2016 $path/$1/datacard_2017 $path/$1/datacard_2018
cp $path/index.php $path/$1
cp $path/index.php $path/$1/datacard_2016
cp $path/index.php $path/$1/datacard_2017
cp $path/index.php $path/$1/datacard_2018

mkdir $path/$1/templatesBkgs_2016 $path/$1/templatesBkgs_2017 $path/$1/templatesBkgs_2018
cp $path/index.php $path/$1/templatesBkgs_2016/.
cp $path/index.php $path/$1/templatesBkgs_2017/.
cp $path/index.php $path/$1/templatesBkgs_2018/.
mv templates/plots/2016/$1/* $path/$1/templatesBkgs_2016/.
mv templates/plots/2017/$1/* $path/$1/templatesBkgs_2017/.
mv templates/plots/2018/$1/* $path/$1/templatesBkgs_2018/.

mv coefficients/matrix_eff/2016/eff_2016_$1_* $path/$1/.
mv coefficients/matrix_eff/2017/eff_2017_$1_* $path/$1/.
mv coefficients/matrix_eff/2018/eff_2018_$1_* $path/$1/.

mv coefficients/matrix_nonfid/2016/nonFid_2016_$1_* $path/$1/.
mv coefficients/matrix_nonfid/2017/nonFid_2017_$1_* $path/$1/.
mv coefficients/matrix_nonfid/2018/nonFid_2018_$1_* $path/$1/.

mv datacard/datacard_2016/hzz4l_*_13TeV_xs_$1_bin*_v* $path/$1/datacard_2016/.
mv datacard/datacard_2017/hzz4l_*_13TeV_xs_$1_bin*_v* $path/$1/datacard_2017/.
mv datacard/datacard_2018/hzz4l_*_13TeV_xs_$1_bin*_v* $path/$1/datacard_2018/.
mv datacard/hzz4l_*_13TeV_xs_$1_bin*_v* $path/$1/.

mv plots/$1/asimov/$1_unfoldwith* $path/$1/.
mv impacts/impacts_*_$1_*_asimov* $path/$1/.
mv LHScans/plots/lhscan_compare_$1_* $path/$1/.

mv fit/commands_$1.py $path/$1/.
mv impacts/commands_impacts_$1.py $path/$1/.
