#!/bin/bash
#PBS -q middleq
#PBS -l mem=5gb|mb|kb,walltime=720:00:00
#HSCHED -s hschedd+soap+human
soap -a query.fastq -D genome.fasta.index -o output.soap -r 2 -m 200 -x 1000
 
