#Calculate the codon decoding rate of RiboSeq on different open reading frames
from collections import defaultdict
import sys

# Check the number of command line arguments
if len(sys.argv) < 4:
    sys.exit("Usage: python _decode_rate_RFS.py inputfile reference 2 or 0 or 1")

input_file_bedgraph = sys.argv[1] + ".bedgraph"
longest_reference = sys.argv[2] + ".fa"
out = int(sys.argv[3])

if out == 2:
    u = "plus"
elif out == 0:
    u = "minus"
elif out == 1:
    u = "inframe"
else:
    sys.exit("Invalid value for the third argument (2=plus 0=minus 1=inframe)")

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
with open(longest_reference, 'r') as file2:
    name = None
    seq = ''
    for line in file2:
        if line.startswith(">"):
            if name is not None:
                dic_cds[name] = seq
            name = line[1:].strip()
            seq = ''
        else:
            seq += line[30:-30]

list_total = []
with open(input_file_bedgraph, 'r') as file1:
    for line in file1:
        a = line.strip().split("\t")
        name = a[0].split("_")[0]
        x = int(a[2]) - int(a[1])
        if name in dic_cds and a[3] != "0" and int(a[1]) >= 130 and x >= 2:
            for count in range(x):
                pos = str(int(a[1]) + count)
                rfs = str(a[3])
                list_total.append(f"{name}\t{pos}\t{rfs}")
        if name in dic_cds and a[3] != "0" and int(a[1]) >= 130 and x == 1:
            list_total.append(f"{name}\t{a[1]}\t{a[3]}")

list_out_frame = []
for line2 in list_total:
    a = line2.strip().split("\t")
    if int(a[1]) % 3 == out and len(dic_cds[a[0]][int(a[1]) - 130:int(a[1]) - 127]) >= 3:
        list_out_frame.append(f"{a[0]}\t{a[1]}\t{a[2]}")

gene_counts = defaultdict(int)
gene_rows = defaultdict(int)
for row in list_out_frame:
    name, _, num = row.split("\t")
    gene_counts[name] += int(num)
    gene_rows[name] += 1

list_standard = []
for line3 in list_out_frame:
    name, pos, num = line3.split("\t")
    norm = float(int(num) / (int(gene_counts[name]) / int(gene_rows[name])))
    if norm <= 50 and int(gene_rows[name]) >= 10:
        list_standard.append(f"{name}\t{dic_cds[name][int(pos) - 130:int(pos) - 127]}\t{norm}")

codon_counts = defaultdict(int)
codon_num = defaultdict(float)
for line4 in list_standard:
    gene, aa_name, aaa_num = line4.strip().split("\t")
    codon = aa_name
    codon_counts[codon] += 1
    codon_num[codon] += float(aaa_num)

with open(f"decoding_rate_{input_file_bedgraph}_{u}.txt", 'w') as fo:
    for k, v in codon_counts.items():
        ratio = codon_num[k] / v
        fo.write(f"{k}({codontable[k]})\t{ratio}\n")
