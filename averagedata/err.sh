#!/bin/bash -x 
#error script
#usage: ./err.sh

EXPDIR=$HOME/Desktop/experimental.gdat
AVGDIR=$HOME/Desktop/7e-5/0.01.gdat
FINDIR=$HOME/Desktop/errors/0.01.gdat

paste $EXPDIR $AVGDIR | head -6 | awk ' { err += (($2 - $8)^2) } END { print err / NR }'> $FINDIR 


