#!usr/bin/env python

import argparse

parser = argparse.ArgumentParser(description="Quality Index Swapping")
parser.add_argument("-f1", "--filepath1", help="R1 file", required=True, type=str)
parser.add_argument("-f2", "--filepath2", help="R2 file", required=True, type=str)
parser.add_argument("-f3", "--filepath3", help="R3 file", required=True, type=str)
parser.add_argument("-f4", "--filepath4", help="R4 file", required=True, type=str)

args = parser.parse_args() #sets arguments as call-able 

f1 = args.filepath1
f2 = args.filepath2
f3 = args.filepath3
f4 = args.filepath4

# f1 = "./1294_S1_L008_R1_001_sample.fastq"
# f2 = "./1294_S1_L008_R2_001_sample.fastq"
# f3 = "./1294_S1_L008_R3_001_sample.fastq"
# f4 = "./1294_S1_L008_R4_001_sample.fastq"

def convert_phred(char):
    n = ord(char) - 33
    return n

def read_freq(file, index, n):
    phredscore_dict = {}

    if index == True:
        bp_len = 8
    else:
        bp_len = 101

    # add up quality scores
    with open(file) as fh:
        i = 0  # line counter
        for line in fh:
            i += 1
            line = line.strip("\n")
            if i % 4 == 0:  # picks out only quality score lines
                sum1 = 0
                for char in line:
                    sum1 += ord(char) - 33
                mean = int(sum1 / bp_len)
                if mean in phredscore_dict:
                    phredscore_dict[mean] += 1
                else:
                    phredscore_dict[mean] = 1
    with open("R" + str(n) + "_part2Plot.txt", "w") as output:
        for key in phredscore_dict:
            output.write(str(key) + "\t" + str(phredscore_dict[key])+"\n")

read_freq(f1, False,1)
read_freq(f2, True,2)
read_freq(f3, True,3)
read_freq(f4, False,4)
