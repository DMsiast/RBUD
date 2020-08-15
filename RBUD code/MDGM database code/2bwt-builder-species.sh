#!/bin/bash
#PBS -q middleq
#PBS -l mem=30gb|mb|kb,walltime=720:00:00
#HSCHED -s hschedd

cd /database and method code/MDGM database/Species dataset/Sequence/NCBI_Species
tar -xvf all.fna.tar
cat * > all.fna
2bwt-builder all.fna
