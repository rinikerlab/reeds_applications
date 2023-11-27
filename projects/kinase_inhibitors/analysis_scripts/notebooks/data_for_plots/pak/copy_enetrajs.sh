#!/bin/bash

ffs="gaff complex_gaff gaff_rest complex_gaff_rest openff complex_openff"

for ff in $ffs; 
do
	echo $ff
	mkdir -p $ff
	
	for i in {1..5};do
		scp euler:/cluster/home/cchampion/work/REEDS/PAK_rerun/${ff}/g_prod_seed${i}/analysis/data/*energies_s1.dat ${ff}/energies_seed${i}.dat 
	done
	scp euler:/cluster/home/cchampion/work/REEDS/PAK_rerun/${ff}/g_prod_seed1/input/*.imd ${ff}/input.imd

done
