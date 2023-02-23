#!/bin/bash

for scriptname in script_*.sh
do
  aux=${scriptname%%.sh}
  id_string=${aux##script_}
 
  bsub -q longparallel -o out_${id_string} -e err_${id_string} -J run_${id_string} -n32 -G sft ${PWD}/${scriptname}
  sleep 2
done

