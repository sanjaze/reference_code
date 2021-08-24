#!/usr/bin/python3

from Bio import SeqIO
import pandas
import os
import sys
from Bio.Seq import Seq
import argparse

parser = argparse.ArgumentParser(description = "Generate readthrough proteins.")
parser.add_argument("-min", help="Minimal length of AA until next stop for seq to be included. Default 5",type=int,action="store",default=5)
parser.add_argument("-max", help="Maximal length of AA after stop to check for Stop. Default 200", type=int,action="store",default=200)
parser.add_argument("-prot", help="Number of AA's to include from protein end. Default 50, w for whole Protein",action="store",default="w")
parser.add_argument("--d", help="Add description to fasta output header.", action="store_true")
parser.add_argument("-o", help="Directory for output",default=os.getcwd(),action="store")
args = parser.parse_args()
print(args.min, args.max, args.prot, args.d, args.o)

#protein file
proteins = {}
for record in SeqIO.parse("/prj/transread/data/human_genome_38.p13/GCF_000001405.39_GRCh38.p13_protein.faa", "fasta"):
    proteins[record.description] = record.seq

genome = {}  
for record in SeqIO.parse("/prj/transread/data/human_genome_38.p13/GCF_000001405.39_GRCh38.p13_genomic.fna", "fasta"):
    genome[record.id.split(" ")[0].strip()] = record.seq

c = 0
m = len(proteins)
outfile = args.o+"extendedProts_min"+str(args.min)+"_max"+str(args.max)+"_prot"+str(args.prot)+"_A_test.fasta"
open(outfile, "a").close()

#test= {}
#test["NP_005908.1 malate dehydrogenase, cytoplasmic isoform MDH1 [Homo sapiens]"] = proteins["NP_005908.1 malate dehydrogenase, cytoplasmic isoform MDH1 [Homo sapiens]"]

for k, v in proteins.items():
  c += 1
  if c%1000 == 0:
    print((c/m)*100, "%")
  matches = []
  id_only = k.strip().split(" ")[0]
  test_ = "ID=cds-"+id_only+";"
  
  with open("/prj/transread/data/human_genome_38.p13/GRCh38.p13_genomic_UTRs.gff", "r") as gff:
    for line in gff:
      if test_ in line: #filters matching cds from gff
        matches.append(line.strip().split("\t"))
   
  matches_df = pandas.DataFrame(matches,columns=["chromo","source","feature","start","stop","score","strand","phase","attributes"])
  matches_df = matches_df[matches_df["feature"]== "CDS"]
  print(matches_df)
  if matches_df.strand[0] == "+":
    start = int(max(matches_df.stop)) #-9
    stop = int(max(matches_df.stop)) + (args.max*3) #+ 990
    dna = genome[matches_df.chromo[0]][start:stop]
    
  elif matches_df.strand[0] == "-":
    stop = int(min(matches_df.start))-1 #+ 8
    start = int(min(matches_df.start)) - ((args.max*3)+1)#- 991
    dna = genome[matches_df.chromo[0]][start:stop].reverse_complement()
    
  else:
    print("No GFF entry found for", k)
    sys.exit()
  
  alternate_as = ["R", "C", "W", "S"]
  prot_seq = str(dna.translate())
  #print(prot_seq)
  if "*" not in str(prot_seq): #exclude seqs with no stop after -max AA
    #print("no stop")
    continue
  
  try:
    out_seq = prot_seq.split("*")[0]
    if len(out_seq) < args.min: #exclude seqs smaller than -min
      #print("too short extension")
      continue
  except IndexError as err:
    print(err, k, prot_seq)
    
  description = " ".join(k.strip().split(" ")[1:])
  continue
  with open(outfile,"a") as f:
    f.write(">"+id_only+"_"+"A"+"|"+description+"\n"+str(v)+"\n") #writes A version
    
  for as_ in alternate_as:
    if args.d:
      header = ">"+id_only+"_"+as_+"|"+description+"\n"
    else:
      header = ">"+id_only+"_"+as_+"\n"
    with open(outfile, "a") as f:
      if args.prot == "w":
        f.write(header+str(v)+as_+str(out_seq)+"\n")
      else:
        f.write(header+str(v[(len(v)-(int(args.prot)-1)):len(v)])+as_+str(out_seq)+"\n")
        
print("Done!")
