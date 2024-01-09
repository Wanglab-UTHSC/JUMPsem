#!/usr/bin/env bash
#BATCH -p talon
#SBATCH --time=50:00:00


MAST_input=$1
MAST_output=$2
MEME_output=$3
ev_num=$(grep ">" $MAST_input | wc -l)
ev=$(echo "0.001 * $ev_num" | bc)

file=$(ls $MEME_output)

for filename in $file
do

meme_file=$MEME_output/$filename/meme.txt
mast $meme_file $MAST_input -ev $ev -o $MAST_output/${filename}_MASTout

done


