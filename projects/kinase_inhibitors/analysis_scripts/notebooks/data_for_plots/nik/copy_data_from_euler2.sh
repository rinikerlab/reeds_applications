#!/bin/bash

ffs="complex_openff"

for ff in $ffs; 
do
	echo $ff
	for i in {1..5};do
		scp euler:/cluster/home/cchampion/work/REEDS/NIK_fixed/${ff}/g3_prod_seed${i}/analysis/free_energy/deltaGs_mbar.npy deltaGs_${ff}_seed${i}.npy
	done
done
