#!/usr/bin/env bash
#BATCH -p talon
#SBATCH --time=50:00:00


MAST_output=$1
Extract_result_path=$2
file=$(ls $MAST_output)

for filename in $file
do

mast_file=$MAST_output/$filename/mast.html
grep_result=$(grep '"name": "psite.' $mast_file | cut -d '"' -f 4 | sed "s/"psite"/$filename/g")


if (grep -q '"name": "psite' $mast_file)
then
   echo -e "${grep_result}\n" >> $2/extract_result.txt
fi

done


