#!/bin/bash

NCOLOR=1
NFLAV=4
NLAT=8
DIM=3

rpn="sun"
together="&"

num_sosia=4
num_sosia=`echo "${num_sosia}-1" | bc -l`

num_beta=1
num_beta=`echo "${num_beta}-1" | bc -l`

num_beta_scalar=8
num_beta_scalar=`echo "${num_beta_scalar}-1" | bc -l`

num_u=1
num_u=`echo "${num_u}-1" | bc -l`


first_beta=0.00
step=0.05

first_beta_scalar=0.80
step_beta_scalar=0.02

first_u=0.00
step_u=0.00

mkdir -p dati_simulation${DIM}d/acceptance_fermion${rpn}${NCOLOR}c${NFLAV}f${NLAT}l
mkdir -p dati_simulation${DIM}d/data_fermion${rpn}${NCOLOR}c${NFLAV}f${NLAT}l
mkdir -p dati_simulation${DIM}d/terma_fermion${rpn}${NCOLOR}c${NFLAV}f${NLAT}l

scriptname=script_fermion${rpn}${NCOLOR}c${NFLAV}f${NLAT}l${first_beta}${first_beta_scalar}${first_u}${num_sosia}.sh

#create the template to execute 
echo "#!/bin/bash" > $scriptname
echo " " >> $scriptname
echo "cd $PWD" >> $scriptname

for aux in `seq 0 1 ${num_beta}`
do
	for gau in `seq 0 1 ${num_beta_scalar}`
	do
                for u_val in `seq 0 1 ${num_u}`
                do	
                	for sosia in `seq 0 1 ${num_sosia}`
                        do
				auxl=`echo "${first_beta}+${step}*${aux}" | bc -l`
				gaul=`echo "${first_beta_scalar}+${step_beta_scalar}*${gau}" | bc -l`
				ul=`echo "${first_u}+${step_u}*${u_val}" | bc -l`
				BETAH=`printf %.5f $aux1`
				GAMMAH=`printf %.5f $gau1`
				echo "./QED_${rpn}${NCOLOR}${NFLAV}${NLAT} ${auxl} ${gaul} 0.0 0.0 0.34 0.0 0 4 0 $RANDOM ${together}" >> $scriptname
			done
		done
	done
done

echo "wait" >> $scriptname
chmod +x $scriptname

icc -O3 *.c -o QED_${rpn}${NCOLOR}${NFLAV}${NLAT} -std=c11 -pedantic-errors -no-multibyte_chars -lm


