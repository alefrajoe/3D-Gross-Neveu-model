#!/bin/bash

for acc_file in `ls dati_simulation*/acc*/*22l*`
do
	echo $acc_file
	tail -n 1 $acc_file | awk '{print $NF}'
done
