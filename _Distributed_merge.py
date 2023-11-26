###Combining multiple sets of results from riboseq data reads distributions into a single plot###
###GongShimin of Yunnan University by2023.11.24###
from collections import defaultdict
import matplotlib.pyplot as plt
import os
import subprocess
dict_list=[]
dic_name={}
o=0
files = os.listdir('.')
for file in files:
    if file.endswith('.txt'):
        #RUN "_Statistical_reads_distribution.py" to get data in txt format
        o+=1
        fname=file[:-18]
        dic_name[o]=fname
        file1 = open(file, 'r')
        dic={}
        for i in file1:
            a=i.strip().split("\t")
            pos=int(a[0])
            va=float(a[1])
            dic[pos]=va
        dict_list.append(dic)

plt.figure(figsize=(10,5))
for i, data in enumerate(dict_list):
    sorted_data = dict(sorted(data.items()))
    plt.plot(list(sorted_data.keys()), list(sorted_data.values()), label=dic_name[i+1])
plt.title('ALL--RPF')
plt.xlabel('position of gene')
plt.ylabel('reads density')
plt.legend()
plt.savefig('ALL_meta_plot.pdf', format='pdf')
