#!/bin/bash


# /cluster/home/cchampion/work/REEDS/NIK_fixed/TI_comparison/openff_repexTI/protein/leg_1_2/seed1/results

mkdir -p solvent complex

for leg in 1_2 1_3 1_4 1_5 1_6;
do
	mkdir -p complex/leg_$leg/
	#mkdir -p solvent/leg_$leg/
	
	for s in 1 2 3 4 5;
	do
		scp -r euler:/cluster/home/cchampion/work/REEDS/NIK_fixed/TI_comparison/openff_repexTI/protein/leg_$leg/seed$s/results complex/leg_$leg/seed$s
		#scp -r euler:/cluster/home/cchampion/work/REEDS/NIK_fixed/TI_comparison/openff_repexTI/solvent/leg_$leg/seed$s/results solvent/leg_$leg/seed$s
	done
done



