#!/bin/bash
#PBS -N Summary
#PBS -l nodes=1:ppn=48
#PBS -q regular
#PBS -l walltime=24:00:00
#PBS -o op.log
#PBS -e op.err
 
set -e

# source ~/.bashrc
export LD_LIBRARY_PATH=/home/daanish/projects/DiploidGenomeAssembly/libdivsufsort/build/lib:$LD_LIBRARY_PATH

LOG=/home/daanish/projects/DiploidGenomeAssembly/summary.txt

cd /home/daanish/projects/DiploidGenomeAssembly/src
make -f makefile_summary clean
make -f makefile_summary

# testing
# MATFASTA=/home/daanish/projects/DiploidGenomeAssembly/test/test0.fasta
# PATFASTA=/home/daanish/projects/DiploidGenomeAssembly/test/test1.fasta

MATFASTA=/home/daanish/data/HG002/chr19_mat.fasta
PATFASTA=/home/daanish/data/HG002/chr19_mat_to_pat.txt

./summary $MATFASTA $PATFASTA $LOG
