#!/bin/sh
#PBS -N nano_FM165_all
#PBS -l select=1:ncpus=16:scratch_local=100gb:mem=25gb
#PBS -l walltime=2:00:00

module add python36-modules-gcc
cp /auto/plzen1/home/nikilangelo/Prague/P26_flaminia_FM165/bp* $SCRATCHDIR

cd $SCRATCHDIR
NanoPlot --fastq bp01.fastq -o nano_BP01 -t16
NanoPlot --fastq bp02.fastq -o nano_BP02 -t16
NanoPlot --fastq bp03.fastq -o nano_BP03 -t16
NanoPlot --fastq bp04.fastq -o nano_BP04 -t16
NanoPlot --fastq bp05.fastq -o nano_BP05 -t16
#NanoPlot --fastq bp06.fastq -o nano_BP06 -t16
#NanoPlot --fastq bp07.fastq -o nano_BP07 -t16
#NanoPlot --fastq bp08.fastq -o nano_BP08 -t16
#NanoPlot --fastq bp09.fastq -o nano_BP09 -t16
#NanoPlot --fastq bp10.fastq -o nano_BP10 -t16
#NanoPlot --fastq bp11.fastq -o nano_BP11 -t16
#NanoPlot --fastq bp12.fastq -o nano_BP12 -t16
cp -r nano* /auto/plzen1/home/nikilangelo/Prague/P26_flaminia_FM165/nanoplot_results