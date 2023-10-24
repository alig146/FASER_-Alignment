#!/bin/sh

for i in {0..24}; do
  i1=$(( i ))
  i2=$(( i+25 ))
  i3=$(( i+50 ))
  i4=$(( i+75 ))
  echo $i  $i2 $i3 $i4
  cat faserkf_misalignmc_template.py | sed   -e "s/666NAME666/${i}/g" -e "s/666NAME1666/${i1}/g"  -e "s/666NAME2666/${i2}/g" -e "s/666NAME3666/${i3}/g" -e "s/666NAME4666/${i4}/g" > faserkf_misalignmc_${i}.py
  cat run_faserkf_misalignmc_template.sh | sed -e "s/666NAME666/${i}/g" > run_faserkf_misalignmc_${i}.sh
  chmod 755 run_faserkf_misalignmc_${i}.sh
done
# for i in {0..24}; do
#   i1=$(( i+100 ))
#   i2=$(( i+25+100 ))
#   i3=$(( i+50+100 ))
#   i4=$(( i+75+100 ))
#   j=$(( i+25 ))
#   echo $j $i2 $i3 $i4
#   cat faserkf_misalignmc_template.py | sed   -e "s/666NAME666/${j}/g" -e "s/666NAME1666/${i1}/g"  -e "s/666NAME2666/${i2}/g" -e "s/666NAME3666/${i3}/g" -e "s/666NAME4666/${i4}/g" > faserkf_misalignmc_${j}.py
#   cat run_faserkf_misalignmc_template.sh | sed -e "s/666NAME666/${j}/g" > run_faserkf_misalignmc_${j}.sh
#   chmod 755 run_faserkf_misalignmc_${j}.sh
# done
