#The random decoding rate ranking of HSC is generated by processing bg files, 
#Judge the reliability of the stop codon ranking of the frameshift boxes
#Statistics that rank multiple times by looping
import random  
import sys  
#python _hsc_rank_random.py iuputfile ref.fa 2or0
arg1 = sys.argv[1]  #input_file_bedgraph
arg2 = sys.argv[2]  #longest_reference
arg3 = sys.argv[3]  #2=plus 0=minus
out=int(arg3)
if out==2:
    u="plus"
else:
    u="minus"
file1 = open(arg1+".bg", 'r')
file2 = open(arg2, 'r')
fw = open("rank_random_"+arg1+"_"+u+".txt", 'a')
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

dic_cds = {}
dic_len = {}
for i in file2:
    a=i.strip()
    if ">" in a:
        name =a[1:]
    else:
        seq =a[30:-30] #each end 30nt not include
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

#out_frame ribosome footprint reads
list_out_frame=[]
for line2 in list_total:
    a = line2.strip().split("\t")
    x = int(dic_len[a[0]])
    if int(a[1])%3==2 and out == 2:
        num = random.randint(1, x-5)
        if int(num)%3 ==1 :
            list_out_frame.append(a[0]+"\t"+str(num)+"\t"+str(a[2]))
        elif int(num)%3 ==0 :
            list_out_frame.append(a[0]+"\t"+str(num+1)+"\t"+str(a[2]))
        elif int(num)%3 ==2 :
            list_out_frame.append(a[0]+"\t"+str(num-1)+"\t"+str(a[2]))
    if int(a[1])%3==0 and out ==0:
        num = random.randint(1, x-5)
        if int(num)%3 ==2 :
            list_out_frame.append(a[0]+"\t"+str(num)+"\t"+str(a[2]))
        elif int(num)%3 ==1 :
            list_out_frame.append(a[0]+"\t"+str(num+1)+"\t"+str(a[2]))
        elif int(num)%3 ==0 :
            list_out_frame.append(a[0]+"\t"+str(num-1)+"\t"+str(a[2]))

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
        list_standard.append(name+"\t"+str(dic_cds[name])[int(pos):int(pos)+3]+"\t"+str(norm))
    
codon_counts={}
codon_num={} 
for line4 in list_standard:
    gene = line4.strip().split("\t")[0]
    aa_name = line4.strip().split("\t")[1]#codontable[line4.strip().split("\t")[1]]
    aaa_num =line4.strip().split("\t")[2]
    if aa_name not in codon_counts :
        codon_counts[aa_name] = float(aaa_num)
    elif aa_name in codon_counts :
        codon_counts[aa_name]  += float(aaa_num)

    if aa_name not in codon_num :
        codon_num[aa_name]=1
    elif aa_name in codon_num :
        codon_num[aa_name]+=1

dict_codon={}
for k1,v1 in codon_num.items():
    dict_codon[k1+"("+codontable[k1]+")"]=float(codon_counts[k1]/v1)
# Calculate the decode rate ranking 
keys = sorted(dict_codon, key=dict_codon.get, reverse=True)  
ranked_d = {key: i+1 for i, key in enumerate(keys)}   
fw.write(str(ranked_d['TAG(*)'])+"\t"+str(ranked_d['TAA(*)'])+"\t"+str(ranked_d['TGA(*)'])+"\t"+str(ranked_d['TAG(*)']+ranked_d['TAA(*)']+ranked_d['TGA(*)'])+"\n")
