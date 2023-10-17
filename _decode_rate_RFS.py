#Calculate the codon decoding rate of RiboSeq on different open reading frames
from collections import defaultdict
import re
import sys  
#python _decode_rate_RFS.py iuputfile reference 2 or 0 or 1
arg1 = sys.argv[1]  #input_file_bedgraph
arg2 = sys.argv[2]  #longest_reference
arg3 = sys.argv[3]  #2=plus 0=minus 1=inframe
out=int(arg3)
if out==2:
    u="plus"
if out==0:
    u="minus"
if out==1:
    u="inframe"
file1 = open(arg1+".bedgraph", 'r')
file2 = open(arg2+".fa", 'r')
fo = open("decoding_rate_"+arg1+"_"+u+".txt", 'w')
codontable = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
        'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W'
    }

def cut_text(text,lenth): 
    textArr = re.findall('.{'+str(lenth)+'}', text) 
    return textArr

dic_cds = {}
dic_len = {}
for i in file2:
    a=i.strip()
    if ">" in a:
        name =a[1:]
    else:
        seq =a[30:-30]#Remove 30 nt each upstream and downstream
        dic_cds[name] = seq
        dic_len[name] =len(seq)

#all ribosome footprint reads
list_total=[]
for line in file1:
    a = line.strip().split("\t")
    name = a[0].split("_")[0]
    x=int(a[2])-int(a[1])
    if name in dic_cds.keys() and a[3] != "0" and int(a[1]) >= 130 and x>=2:
        for count in range(x):
            pos=str(int(a[1])+count)
            rfs=str(a[3])
            list_total.append(name+"\t"+pos+"\t"+rfs)
    if name in dic_cds.keys() and a[3] != "0" and int(a[1]) >= 130 and x==1:
        list_total.append(name+"\t"+str(a[1])+"\t"+str(a[3]))

#Select Open Reading Frame of ribosome footprint reads
list_out_frame=[]
for line2 in list_total:
    a = line2.strip().split("\t")
    if int(a[1])%3 ==out and len(str(dic_cds[a[0]])[int(a[1])-130:int(a[1])-127]) >= 3:
        list_out_frame.append(a[0]+"\t"+str(a[1])+"\t"+str(a[2]))

gene_counts = {}
gene_rows = {}
for row in list_out_frame: 
    name = row.split("\t")[0]
    num = row.split("\t")[2]
    if name not in gene_counts:
        gene_counts[name] = int(num)
    else:
        gene_counts[name]  += int(num)
    if name not in gene_rows:
        gene_rows[name] =1
    else:
        gene_rows[name] +=1

list_standard=[]
    
for line3 in list_out_frame:
    name = line3.split("\t")[0]
    pos = line3.split("\t")[1]
    num = line3.split("\t")[2]
    norm =float(int(num)/(int(gene_counts[name])/int(gene_rows[name])))
    if norm <=50 and int(gene_rows[name])>=10:
        list_standard.append(name+"\t"+str(dic_cds[name])[int(pos)-130:int(pos)-127]+"\t"+str(norm))
    
codon_counts={}
codon_num={} 
for line4 in list_standard:
    gene = line4.strip().split("\t")[0]
    aa_name = line4.strip().split("\t")[1]
    aaa_num =line4.strip().split("\t")[2]
    codon=aa_name
    if codon not in codon_counts:
        codon_counts[codon] = 1
    else:
        codon_counts[codon]  += 1
    if codon not in codon_num:
        codon_num[codon] =float(aaa_num)
    else:
        codon_num[codon] +=float(aaa_num)
for k,v in codon_counts.items():
    ratio=codon_num[k]/v
    fo.write(k+"("+codontable[k]+")"+"\t"+str(ratio)+"\n")

        


