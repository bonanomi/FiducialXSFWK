#!/bin/bash
# $1 -> true or false to unblind

declare -a obs=(
"mass4l noJES"
# "mass4l_zzfloating noJES"
# # "njets_pt30_eta4p7 JES"
# "pT4l noJES"
# "rapidity4l noJES"
# "costhetaZ1 noJES"
# "costhetaZ2 noJES"
# "phi noJES"
# "phistar noJES"
# "costhetastar noJES"
# "massZ1 noJES"
# "massZ2 noJES"
# # "pTj1 JES"
# # "pTHj JES"
# # "mHj JES"
# # "pTj2 JES"
# # "mjj JES"
# # "absdetajj JES"
# # "dphijj JES"
# # "pTHjj JES"
# # "TCjmax JES"
# # "TBjmax JES"
# "D0m noJES"
# "Dcp noJES"
# "D0hp noJES"
# "Dint noJES"
# "DL1 noJES"
# "D0m noJES"
# "Dcp noJES"
# "DL1Zg noJES"
# # "rapidity4l vs pT4l JES"
# # "njets_pt30_eta4p7 vs pT4l JES"
# # "pTj1 vs pTj2 JES"
# # "pT4l vs pTHj JES"
# "massZ1 vs massZ2 noJES"
# "TCjmax vs pT4l JES"
)
source /opt/exp_soft/cms/t3/t3setup

echo
echo "------------------------------------------------------------"

for i in "${!obs[@]}"; do
  name=($(echo ${obs[$i]} | cut -d" " -f1))
  name_folder=$name
  if [ $(echo ${obs[$i]} | cut -d" " -f2) == "vs" ]; then
    name=$name" vs "$(echo ${obs[$i]} | cut -d" " -f3)
    name_folder=$name_folder"_"$(echo ${obs[$i]} | cut -d" " -f3)
  fi
  jes=($(echo ${obs[$i]} | rev | cut -d" " -f1 | rev))
  echo
  echo "name of the observable:" $name
  echo "name of folders:" $name_folder
  echo $jes
  if [ -d $name_folder ]; then
    rm -rf $name_folder
  fi
  mkdir $name_folder
  cd $name_folder
  cp ../batchScript.sh .
  sed -i "s/OBS/$name_folder/g" batchScript.sh
  sed -i "s/VAR/$name/g" batchScript.sh
  sed -i "s/BIN/$bins/g" batchScript.sh

  if [ $1 == true ]; then
    sed -i "s/asimov/data/g" batchScript.sh
    sed -i "s/UNBLIND/--unblind/g" batchScript.sh
  else
    sed -i "s/UNBLIND//g" batchScript.sh
  fi

  #JES flag
  if [ $jes == "JES" ]; then
    sed -i "s/SETTING/true/g" batchScript.sh
  else
    sed -i "s/SETTING/false/g" batchScript.sh
  fi

  #Compute AC acceptances using AC samples
  #-----FIXME:Since AC samples have been processed for 2016 only, little trick to postpone the problem
  #-----FIXME:RunCoefficients is run just for 2016
  if [ $name_folder == D0m ]; then
    CoeffAC="python RunCoefficients.py --obsName 'D0m' --year '2016' --AC_onlyAcc --AC_hypothesis '0M'"
    CoeffACbis=""
    plotsv4="python producePlots_v4.py --obsName 'D0m' --year 'Full'"
    LLscanv4="python plot_LLScan.py --obsName 'D0m' --year 'Full' --v4"
  elif [ $name_folder == Dcp ]; then
    CoeffAC="python RunCoefficients.py --obsName 'Dcp' --year '2016' --AC_onlyAcc --AC_hypothesis '0Mf05ph0'"
    CoeffACbis=""
    plotsv4="python producePlots_v4.py --obsName 'Dcp' --year 'Full'"
    LLscanv4="python plot_LLScan.py --obsName 'Dcp' --year 'Full' --v4"
  elif [ $name_folder == D0hp ]; then
    CoeffAC="python RunCoefficients.py --obsName 'D0hp' --year '2016' --AC_onlyAcc --AC_hypothesis '0PH'"
    CoeffACbis=""
    plotsv4="python producePlots_v4.py --obsName 'D0hp' --year 'Full'"
    LLscanv4="python plot_LLScan.py --obsName 'D0hp' --year 'Full' --v4"
  elif [ $name_folder == Dint ]; then
    CoeffAC="python RunCoefficients.py --obsName 'Dint' --year '2016' --AC_onlyAcc --AC_hypothesis '0PHf05ph0'"
    CoeffACbis=""
    plotsv4="python producePlots_v4.py --obsName 'Dint' --year 'Full'"
    LLscanv4="python plot_LLScan.py --obsName 'Dint' --year 'Full' --v4"
  elif [ $name_folder == DL1 ]; then
    CoeffAC="python RunCoefficients.py --obsName 'DL1' --year '2016' --AC_onlyAcc --AC_hypothesis '0L1'"
    CoeffACbis="python RunCoefficients.py --obsName 'DL1' --year '2016' --AC_onlyAcc --AC_hypothesis '0L1f05ph0'"
    plotsv4="python producePlots_v4.py --obsName 'DL1' --year 'Full'"
    LLscanv4="python plot_LLScan.py --obsName 'DL1' --year 'Full' --v4"
  elif [ $name_folder == massZ1 ]; then
    CoeffAC=""
    CoeffACbis=""
    plotsv4="python producePlots_v4.py --obsName 'massZ1' --year 'Full'"
    LLscanv4="python plot_LLScan.py --obsName 'massZ1' --year 'Full' --v4"
  elif [ $name_folder == massZ2 ]; then
    CoeffAC=""
    CoeffACbis=""
    plotsv4="python producePlots_v4.py --obsName 'massZ1' --year 'Full'"
    LLscanv4="python plot_LLScan.py --obsName 'massZ1' --year 'Full' --v4"
  else
    CoeffAC=""
    CoeffACbis=""
    plotsv4=""
    LLscanv4=""
  fi
  sed -i "s/FIRST/$CoeffAC/g" batchScript.sh
  sed -i "s/SECOND/$CoeffACbis/g" batchScript.sh
  sed -i "s/THIRD/$plotsv4/g" batchScript.sh
  sed -i "s/FOURTH/$LLscanv4/g" batchScript.sh


  if [ $name == mass4l_zzfloating ]; then
    zzfloating="--m4lLower 105 --m4lUpper 160"
  else
    zzfloating=""
  fi
  sed -i "s/ZZFLOATING/$zzfloating/g" batchScript.sh

  /opt/exp_soft/cms/t3/t3submit batchScript.sh
  cd ..
  echo
done