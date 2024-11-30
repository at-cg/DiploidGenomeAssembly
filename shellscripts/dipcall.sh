#!/bin/bash
#PBS -N dipCall
#PBS -l nodes=1:ppn=48
#PBS -q regular
#PBS -l walltime=4:00:00
#PBS -o op.log
#PBS -e op.err
 
EXE=/home/daanish/apps/dipcall.kit/run-dipcall
DIPCALL=/home/daanish/apps/dipcall.kit

CHR19MAT=/home/daanish/data/HG002/chr19_mat.fasta
CHR19PAT=/home/daanish/data/HG002/chr19_pat.fasta
source ~/.bashrc
cd /home/daanish/projects/repeatstatistics
samtools faidx $CHR19MAT
$EXE prefix $CHR19MAT $CHR19PAT $CHR19MAT > prefix.mak
make -j1 -f prefix.mak
#~/miniforge3/bin/tabix -p vcf prefix.dip.vcf.gz
