#!/bin/bash
# $1 -> true or false to unblind

declare -a obs=(
"mass4l |105|140|"
"pT4l |0|10|20|30|45|60|80|120|200|3000|"
"rapidity4l |0.0|0.15|0.3|0.45|0.6|0.75|0.9|1.2|1.6|2.5|"
"D0m |0.0|0.4|0.5|0.6|0.7|0.8|0.9|1.0|"
"Dcp |-0.75|-0.25|-0.1|0|0.1|0.25|0.75|")

source /opt/exp_soft/cms/t3/t3setup

for i in "${!obs[@]}"; do
  bins=($(echo ${obs[$i]} | rev | cut -d" " -f1  | rev))
  name=($(echo ${obs[$i]} | rev | cut -d"_" -f2-  | rev))
  echo $name
  echo $bins
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
  /opt/exp_soft/cms/t3/t3submit batchScript.sh
  cd ..
  echo
done
