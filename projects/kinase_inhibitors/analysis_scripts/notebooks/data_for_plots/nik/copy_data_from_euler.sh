#!/bin/bash

ffs="gaff complex_gaff openff complex_openff hybrid_openff complex_hybrid_openff"

#ffs="gaff openff hybrid_openff"

ffs="complex_hybrid_openff"

for ff in $ffs; 
do
	echo $ff
	for i in {1..5};do
		scp euler:/cluster/home/cchampion/work/REEDS/NIK_fixed/${ff}/g_prod_seed${i}/analysis/free_energy/deltaGs_mbar.npy deltaGs_${ff}_seed${i}.npy
	done
done
