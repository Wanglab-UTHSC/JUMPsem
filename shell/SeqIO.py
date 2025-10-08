#pip install biopython

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from random import randint
import sys

# python3 SeqIO.py output input_file

output = sys.argv[1] # output file with path
input_file = sys.argv[2] # raw input file with path, csv file

f1 = open(output, "w")

fasta = dict()
with open("/home/groupdirs/wanglab/dehui_kong/Database/mouse/proteomeSeq/mouse.fasta", "r") as squence:
    for record in SeqIO.parse(squence, "fasta"):
        fasta[record.id.split('|')[1]] = record

lengthdict = dict()
with open(input_file) as seqlengths:
    next(seqlengths)
    for line in seqlengths:
        split_IDlength  = line.strip().split(",")
        sub_prot = split_IDlength[0]
        psite = str(split_IDlength[2])
        moid = int(split_IDlength[2][1:])

        lengthdict[(sub_prot+"_"+psite)] = int(moid)
        #print(lengthdict)

#f1.write("SUB_psite" + "\t" + "Peptide" + "\n")
for k, v in lengthdict.items():
    if fasta.get(k.split("_")[0]) is None:
        continue
    #print(k)
        #print(fasta[k].seq[v-5:v+5])
    #f1.write(">psite." + format(k) + "\t")
    f1.write(format(k) + "\t")
    f1.write(str(fasta[k.split("_")[0]].seq[v-9:v+8]) + "\n")
    #print(f1)

f1.close()

