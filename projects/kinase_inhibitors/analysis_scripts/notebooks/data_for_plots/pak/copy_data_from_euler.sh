#!/bin/bash

ffs="gaff complex_gaff gaff_rest complex_gaff_rest openff complex_openff"

for ff in $ffs; 
do
	echo $ff
	for i in {1..5};do
		scp euler:/cluster/home/cchampion/work/REEDS/PAK_rerun/${ff}/g_prod_seed${i}/analysis/free_energy/deltaGs_mbar.npy deltaGs_${ff}_seed${i}.npy
	done
done
