#!/bin/bash
#PBS -N dipCall
#PBS -l nodes=1:ppn=48
#PBS -q regular
#PBS -l walltime=24:00:00
#PBS -o op.log
#PBS -e op.err
 
source ~/.bashrc

EXE=/home/daanish/projects/repeatstatistics/code/mat_triple

cd /home/daanish/projects/repeatstatistics/code
rm op.log op.err
/usr/bin/time -v $EXE


