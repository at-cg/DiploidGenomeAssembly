#!/bin/bash
#PBS -N RepeatStats
#PBS -l nodes=1:ppn=48
#PBS -q regular
#PBS -l walltime=24:00:00
#PBS -o op.log
#PBS -e op.err
 
source ~/.bashrc

LOG=/home/daanish/projects/DiploidGenomeAssembly/src/log.txt

cd /home/daanish/projects/DiploidGenomeAssembly/src
make -f makefile_intradouble clean
make -f makefile_intradouble

echo "Starting test" >> $LOG

for i in {1..100}; do
    EXE=/home/daanish/projects/DiploidGenomeAssembly/src/test.py
    python3 $EXE
    INPUT=/home/daanish/projects/DiploidGenomeAssembly/src/test.fasta
    OUTPUT=/home/daanish/projects/DiploidGenomeAssembly/src/test_intradoublerepeat.txt

    export LD_LIBRARY_PATH=/home/daanish/projects/DiploidGenomeAssembly/libdivsufsort/build/lib:$LD_LIBRARY_PATH
    ./intradouble $INPUT $OUTPUT

    echo "Iteration $i finished" >> $LOG
done

echo "Ending test" >> $LOG

