#!/bin/bash

DATADIR=$HOME/ork/src/github.com/user/fcer1-gamma-phosphorylation/base+degrad/data/7e-3
cat $DATADIR/*/*.gdat | grep -v '#' | awk '{cnt[$1]++; p0[$1]+=$2; p1[$1]+=$3; p2[$1]+=$4; p3[$1]+=$5; p4[$1]+=$6; p5[$1]+=$7; p6[$1]+=$8; p7[$1]+=$9; p8[$1]+=$10} END { for(i in cnt) print i,p0[i]/cnt[i],p1[i]/cnt[i],p2[i]/cnt[i],p3[i]/cnt[i],p4[i]/cnt[i],p5[i]/cnt[i],p6[i]/cnt[i],p7[i]/cnt[i],p8[i]/cnt[i] }' | sort -g > $DATADIR/avg.gdat
