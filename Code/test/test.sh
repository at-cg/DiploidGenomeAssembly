#!/bin/bash
#PBS -N RepeatStats-Test
#PBS -l nodes=1:ppn=48
#PBS -q regular
#PBS -l walltime=24:00:00
#PBS -o op.log
#PBS -e op.err
 
set -e

# source ~/.bashrc
export LD_LIBRARY_PATH=/home/daanish/projects/DiploidGenomeAssembly/libdivsufsort/build/lib:$LD_LIBRARY_PATH

cd /home/daanish/projects/DiploidGenomeAssembly/test
make clean
make

./chk