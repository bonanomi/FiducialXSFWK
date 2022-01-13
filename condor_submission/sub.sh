#!/bin/bash
# $1 -> true or false to unblind

declare -a obs=(
"mass4l |105|140| noJES")
#"pT4l |0|10|20|30|45|60|80|120|200|3000| noJES"
#"rapidity4l |0.0|0.15|0.3|0.45|0.6|0.75|0.9|1.2|1.6|2.5| noJES"
#"D0m |0.0|0.4|0.5|0.6|0.7|0.8|0.9|1.0| noJES"
#"Dcp |-0.75|-0.25|-0.1|0.0|0.1|0.25|0.75| noJES"
# "pTj1 |-2|30|55|95|200| JES")

source /opt/exp_soft/cms/t3/t3setup

for i in "${!obs[@]}"; do
  name=($(echo ${obs[$i]} | cut -d" " -f1))
  bins=($(echo ${obs[$i]} | cut -d" " -f2-))
  jes=($(echo ${obs[$i]} | rev | cut -d" " -f1 | rev))
  echo $name
  echo $bins
  echo $jes
  if [ -d $name ]; then
    rm -rf $name
  fi
  mkdir $name
  cd $name
  cp ../batchScript.sh .
  sed -i "s/OBS/$name/g" batchScript.sh
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
  if [ $name == D0m ]; then
    CoeffAC="python RunCoefficients.py --obsName 'D0m' --obsBins '$bins' --year '2016' --AC_onlyAcc --AC_hypothesis '0M'"
    CoeffACbis=""
    plotsAC="python producePlots_v4.py --obsName 'D0m' --obsBins '$bins' --year 'Full'"
    LLscanAC="python plot_LLScan.py --obsName 'D0m' --obsBins '$bins' --year 'Full' --v4"
  elif [ $name == Dcp ]; then
    CoeffAC="python RunCoefficients.py --obsName 'Dcp' --obsBins '$bins' --year '2016' --AC_onlyAcc --AC_hypothesis '0Mf05ph0'"
    CoeffACbis=""
    plotsAC="python producePlots_v4.py --obsName 'Dcp' --obsBins '$bins' --year 'Full'"
    LLscanAC="python plot_LLScan.py --obsName 'Dcp' --obsBins '$bins' --year 'Full' --v4"
  elif [ $name == D0hp ]; then
    CoeffAC="python RunCoefficients.py --obsName 'D0hp' --obsBins '$bins' --year '2016' --AC_onlyAcc --AC_hypothesis '0PH'"
    CoeffACbis=""
    plotsAC="python producePlots_v4.py --obsName 'D0hp' --obsBins '$bins' --year 'Full'"
    LLscanAC="python plot_LLScan.py --obsName 'D0hp' --obsBins '$bins' --year 'Full' --v4"
  elif [ $name == Dint ]; then
    CoeffAC="python RunCoefficients.py --obsName 'Dint' --obsBins '$bins' --year '2016' --AC_onlyAcc --AC_hypothesis '0PHf05ph0'"
    CoeffACbis=""
    plotsAC="python producePlots_v4.py --obsName 'Dint' --obsBins '$bins' --year 'Full'"
    LLscanAC="python plot_LLScan.py --obsName 'Dint' --obsBins '$bins' --year 'Full' --v4"
  elif [ $name == DL1 ]; then
    CoeffAC="python RunCoefficients.py --obsName 'DL1' --obsBins '$bins' --year '2016' --AC_onlyAcc --AC_hypothesis '0L1'"
    CoeffACbis="python RunCoefficients.py --obsName 'DL1' --obsBins '$bins' --year '2016' --AC_onlyAcc --AC_hypothesis '0L1f05ph0'"
    plotsAC="python producePlots_v4.py --obsName 'DL1' --obsBins '$bins' --year 'Full'"
    LLscanAC="python plot_LLScan.py --obsName 'DL1' --obsBins '$bins' --year 'Full' --v4"
  else
    CoeffAC=""
    CoeffACbis=""
    plotsAC=""
    LLscanAC=""
  fi

  sed -i "s/FIRST/$CoeffAC/g" batchScript.sh
  sed -i "s/SECOND/$CoeffACbis/g" batchScript.sh
  sed -i "s/THIRD/$plotsAC/g" batchScript.sh
  sed -i "s/FOURTH/$LLscanAC/g" batchScript.sh

  /opt/exp_soft/cms/t3/t3submit batchScript.sh
  cd ..
  echo
done
