###Reads distribution statistics of Ribo-seq data on genes###
###Before and after each 100bp UTR region was recorded in 10 equal parts###
###CDS region recorded in 100 equal parts###

###Run the script will count all the bg files under the current path###
###generate pdf graphics and the corresponding txt text data###

###GongShimin of Yunnan University by2023.11.24###
from collections import defaultdict
import matplotlib.pyplot as plt
import os
import subprocess

files = os.listdir('.')
for file in files:
    if file.endswith('.bg'):
        fname=file[:-3]
        file1 = open(file, 'r')
        fw = open(fname+'_mateplot_data.txt', 'w')
        cut=10#Setting the outlier threshold#
        least=50#Set the gene to at least read#

        dic_len=defaultdict(list)
        list_total=[]
        for line in file1:
            a = line.strip().split("\t")
            name = a[0]
            p = int(a[2])
            dic_len[name].append(p)
            x=int(a[2])-int(a[1])
            if a[3] != "0" and x>=2:
                for count in range(x):
                    pos=str(int(a[1])+count)
                    rfs=str(a[3])
                    list_total.append(name+"\t"+pos+"\t"+rfs)
            if  a[3] != "0" and x==1:
                list_total.append(name+"\t"+str(a[1])+"\t"+str(a[3]))
        dic_max={}
        for k,v in dic_len.items():
            m=int(max(v))-200
            dic_max[k]=m

        dic_len=defaultdict(list)
        dic_reads={}
        for line2 in list_total:
            a = line2.strip().split("\t")
            na=a[0]
            pos=int(a[1])
            allreads=int(a[2])
            if na not in dic_reads and pos>=100 and pos<=dic_max[na]+100:
                dic_reads[na]=allreads
            elif na in dic_reads and pos>=100 and pos<=dic_max[na]+100:
                dic_reads[na]+=allreads

        dic_rpkm={}
        for k1,v1 in dic_reads.items():
            rpk=float(int(v1)/int(dic_max[k1]))
            dic_rpkm[k1]=rpk
        dic_distrbution={}
        list_all_ratio=[]
        for line3 in list_total:
            a = line3.strip().split("\t")
            na=a[0]
            if na in dic_rpkm.keys():
                ratio=float(int(a[2])/int(dic_reads[na]))*dic_rpkm[na]
                if 100<=int(a[1])<=dic_max[na]+100 and dic_reads[na]>least and ratio<cut:
                    pos=int((int(a[1])-100)*100/int(dic_max[na]))
                    if pos not in dic_distrbution:
                        dic_distrbution[pos]=ratio
                    else:
                        dic_distrbution[pos]+=ratio
                if 100>int(a[1]) and dic_reads[na]>least and ratio<cut:
                    pos=int(int(a[1])/10)-10
                    if pos not in dic_distrbution:
                        dic_distrbution[pos]=ratio
                    else:
                        dic_distrbution[pos]+=ratio
                if int(a[1])>dic_max[na]+100 and dic_reads[na]>least and ratio<cut:
                    pos=int((dic_max[na]+200-int(a[1]))/10)+100
                    if pos not in dic_distrbution:
                        dic_distrbution[pos]=ratio
                    else:
                        dic_distrbution[pos]+=ratio
                        
        ma=0
        for k2,v2 in dic_distrbution.items(): 
            ma+=v2  
        dic_meta_data={}
        for k3,v3 in dic_distrbution.items():   
            dic_meta_data[k3]=float(v3/ma)
            fw.write(str(k3)+"\t"+str(float(v3/ma))+"\n")

        sorted_data = dict(sorted(dic_meta_data.items()))
        keys = sorted_data.keys()
        values = sorted_data.values()
        
        plt.figure(figsize=(10,5))
        plt.plot(keys, values)
        plt.title(str(fname)+'--RPF')
        plt.xlabel('position of gene')
        plt.ylabel('reads density')
        plt.savefig(str(fname)+'_meta_plot.pdf', format='pdf')
        print("===finish " +fname+"===")
 

