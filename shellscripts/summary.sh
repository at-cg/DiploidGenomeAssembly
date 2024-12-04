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

ERR=10000
CHR=19
# ERR=1000

# Data preparation
cd /home/daanish/data/HG002/
if [ ! -d "chr$CHR/p$ERR" ]; then
    mkdir -p "chr$CHR/p$ERR"
fi

# testing
# MATFASTA=/home/daanish/projects/DiploidGenomeAssembly/test/test0.fasta
# PATFASTA=/home/daanish/projects/DiploidGenomeAssembly/test/test1.fasta

MATFASTA="/home/daanish/data/HG002/chr$CHR/mat.fasta"
PATFASTA="/home/daanish/data/HG002/chr$CHR/p$ERR/pat.txt"
EXE=/home/daanish/projects/DiploidGenomeAssembly/src/prepare.py

# python3 $EXE -mpath $MATFASTA -ppath $PATFASTA -err $ERR

# Data analysis
cd /home/daanish/projects/DiploidGenomeAssembly/results
if [ ! -d "chr$CHR/p$ERR" ]; then
    mkdir -p "chr$CHR/p$ERR"
fi

OUT="/home/daanish/projects/DiploidGenomeAssembly/results/chr$CHR/p$ERR/"

cd ../src
# make -f makefile_summary clean
# make -f makefile_summary

# ./summary $MATFASTA $PATFASTA $OUT

# Data Plotting
EXE=/home/daanish/projects/DiploidGenomeAssembly/src/plot.py
GPATH=$OUT/gapstats.txt
MPATH=$OUT/interleaved_mat_stats.txt
PPATH=$OUT/interleaved_pat_stats.txt
SPATH=$OUT/specialinterleaved_stats.txt
OPATH="/home/daanish/projects/DiploidGenomeAssembly/plot/chr$CHR-p$ERR.png"
EP=0.01

python3 $EXE -gpath $GPATH -mpath $MPATH -ppath $PPATH -spath $SPATH -ep $EP -opath $OPATH