#!/bin/bash -x 
#run script
#usage: ./run.sh # (of simulations)

BNGDIR=$HOME/work/BioNetGen-2.2.6-stable
DATADIR=$HOME/work/alex/data
RUNS=$1

rm -rf $DATADIR/*

for i in `seq 1 $RUNS` ; do
	D=$DATADIR/$i
	mkdir -p $D
	perl $BNGDIR/BNG2.pl --outdir $D base.bngl 2>&1 > $D/out.log &
done