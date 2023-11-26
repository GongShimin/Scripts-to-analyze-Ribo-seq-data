###Process bedgraph file of Ribo-seq data to create a sliding box to filter/ 
###for segments with 10 consecutive boxes in which the code/ 
###shift rate is 2 standard deviations above the overall level
###GongShimin of Yunnan University by2023.11.13###
import numpy as np
import argparse
from collections import defaultdict

def process_bedgraph(bedgraph_file):
    with open(bedgraph_file, "r") as infile:
        bedgraph_data = infile.readlines()
    signal_values = {}
    for line in bedgraph_data:
        parts = line.strip().split("\t")
        name = parts[0]
        start = int(parts[1])
        end = int(parts[2])
        signal = int(parts[3])
        for i in range(start, end):
            signal_values[(name, i)] = signal
    return signal_values

def calculate_ratios(signal_data, window_size, skip_size, min_nonzero_count):
    gene_ratios = []
    for i in range(0, len(signal_data) - window_size, skip_size):
        window = signal_data[i:i + window_size]
        # frame shift ratio 
        nonzero_count = len([x for x in window if x != 0])
        if nonzero_count > min_nonzero_count:
            total_sum = sum(window)
            del window[::3]
            nonzero_sum = sum(window)
            ratio = nonzero_sum / total_sum
            gene_ratios.append(ratio)
    return gene_ratios

def main():
    parser = argparse.ArgumentParser(description="Process BedGraph file and calculate ratios.")
    parser.add_argument("bedgraph_file", help="Input BedGraph file")
    parser.add_argument("--w", type=int, default=90, help="Window size for sliding")
    parser.add_argument("--s", type=int, default=3, help="Skip size for sliding")
    parser.add_argument("--m", type=int, default=20, help="Minimum count of nonzero values in a window")
    parser.add_argument("--o", default="FSR_SWIPE_POS.txt", help="Output file for results")
    args = parser.parse_args()
    signal_values = process_bedgraph(args.bedgraph_file)
    
    all_ratios = []
    gene_pos_ratios = {}
    for name in set([entry[0] for entry in signal_values.keys()]):
        signal_data = [signal_values.get((name, i), 0) for i in range(1, max([entry[1] for entry in signal_values.keys()]) + 1)]
        signal_data = signal_data[105:-100]
        gene_ratios = calculate_ratios(signal_data, args.window_size, args.skip_size, args.min_nonzero_count)
        all_ratios.extend(gene_ratios)
        for i, ratio in enumerate(gene_ratios):
            position = i * args.skip_size
            gene_pos_ratios[name + ":" + str(position)] = ratio

    std = np.std(all_ratios)
    average = np.mean(all_ratios) 
    DThreshold_value = average + (2 * std)
    UThreshold_value = average + (3 * std)
    gene_enrich = defaultdict(list)
    for k, v in gene_pos_ratios.items():
        gene = k.split(":")[0]
        position = k.split(":")[1]
        if DThreshold_value < v < UThreshold_value:
            gene_enrich[gene].append(position)
        #Setting the range for filtering

    with open(args.o, "w") as fo:
        for k1, v1 in gene_enrich.items():
            res = []
            for i in range(len(v1)):
                if not res:
                    res.append([v1[i]])
                elif int(v1[i - 1]) + 6 >= int(v1[i]):
                    res[-1].append(v1[i])
                else:
                    res.append([v1[i]])
            for i in res:
                if len(i) > 10:
                    fo.write(k1 + ":" + i[0] + "-" + i[-1] + "\n")
                    #Record information on screening segments
if __name__ == "__main__":
    main()
