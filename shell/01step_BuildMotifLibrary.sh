#!/usr/bin/env bash
#BATCH -p talon
#SBATCH --time=50:00:00


kin_subtable=$1
MEME_output=$2

for file in $kin_subtable/*.fa
do

filename=$(basename $file .fa)
meme $file -nmotifs 1000 -minw 7 -maxw 15 -mod anr -o $MEME_output/$filename

done









