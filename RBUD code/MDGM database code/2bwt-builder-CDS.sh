#!/bin/bash
#PBS -q middleq
#PBS -l mem=30gb|mb|kb,walltime=720:00:00
#HSCHED -s hschedd

cd /database and method code/MDGM database/Functional dataset/Sequence/NCBI_CDS
tar -xvf all.ffn.tar
cat * > all.ffn
2bwt-builder all.ffn
