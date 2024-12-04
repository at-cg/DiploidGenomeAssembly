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

LOG=/home/daanish/projects/DiploidGenomeAssembly/log.txt

cd /home/daanish/projects/DiploidGenomeAssembly/src
# make -f makefile_intradouble clean
# make -f makefile_intradouble
make clean
make

# echo "Starting test" >> $LOG

# echo "Starting intra-double" >> $LOG

# EXE=/home/daanish/projects/DiploidGenomeAssembly/test/test.py

# for i in {1..100}; do
#     INPUT=/home/daanish/projects/DiploidGenomeAssembly/test/test.fasta
#     python3 $EXE -len 500 -opath1 $INPUT

#     OUTPUT=/home/daanish/projects/DiploidGenomeAssembly/test/test_intradoublerepeat.txt
#     ./intradouble $INPUT $OUTPUT

#     echo "Iteration $i finished" >> $LOG
# done

echo "Starting intra-triple" >> $LOG
INPUT=/home/daanish/data/HG002/chr19_mat_to_pat.txt
OUTPUT=/home/daanish/projects/DiploidGenomeAssembly/test/test_intratriplerepeat.txt
./intratriple $INPUT $OUTPUT

# EXE=/home/daanish/projects/DiploidGenomeAssembly/test/test.py

# for i in {1..100}; do
    # INPUT=/home/daanish/projects/DiploidGenomeAssembly/test/test.fasta
    # python3 $EXE -len 100 -opath1 $INPUT

#     OUTPUT=/home/daanish/projects/DiploidGenomeAssembly/test/test_intratriplerepeat.txt
#     ./intratriple $INPUT $OUTPUT

#     echo "Iteration $i finished" >> $LOG
# done

# echo "Starting inter-double" >> $LOG

# for i in {1..100}; do
#     INPUT0=/home/daanish/projects/DiploidGenomeAssembly/test/test0.fasta
#     INPUT1=/home/daanish/projects/DiploidGenomeAssembly/test/test1.fasta
#     python3 $EXE -num 2 -len 100 -snp 10 -opath1 $INPUT0 -opath2 $INPUT1

#     OUTPUT=/home/daanish/projects/DiploidGenomeAssembly/test/test_interdoublerepeat.txt
#     ./interdouble $INPUT0 $INPUT1 $OUTPUT

#     echo "Iteration $i finished" >> $LOG
# done

# echo "Ending test" >> $LOG

