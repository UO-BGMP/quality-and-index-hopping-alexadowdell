#!usr/bin/env python

######################################## Format argparse ##########################################################
# 1. Write script parameters to take in the four files as well as the index.tsv file and a min quality score cut off.

import argparse
import numpy as np

parser = argparse.ArgumentParser(description="Quality Index Swapping")
parser.add_argument("-f1", "--filepath1", help="R1 file", required=True, type=str)
parser.add_argument("-f2", "--filepath2", help="R2 file", required=True, type=str)
parser.add_argument("-f3", "--filepath3", help="R3 file", required=True, type=str)
parser.add_argument("-f4", "--filepath4", help="R4 file", required=True, type=str)
parser.add_argument("-f5", "--filepath5", help="Indexes.txt file", required=True, type=str)
parser.add_argument("-c", "--cutoff", help="Min quality score cutoff", required=True, type=int)


args = parser.parse_args() #sets arguments as call-able 
f1 = args.filepath1
f2 = args.filepath2
f3 = args.filepath3
f4 = args.filepath4
f5 = args.filepath5
c = args.cutoff

# f1 = "./1294_S1_L008_R1_001_sample.fastq"
# f2 = "./1294_S1_L008_R2_001_sample.fastq"
# f3 = "./1294_S1_L008_R3_001_sample.fastq"
# f4 = "./1294_S1_L008_R4_001_sample.fastq"
# f5 = "./indexes2.txt"
# c = 0

def convert_phred(letter):
    n = ord(letter) - 33
    return n


######################## Build dictionary containing every 4 nucleotide seq combination ############################

# Every possible combo of 2, 8 base long indexes as key and value = counter(increases by 1 every time seq is 
# uncovered in index reads R2, R3)
nIndexDict = {}
indexComboDict = {}
unmatchedIndexDict = {}
list_index_reads = [] #list holding index reads from said index file

with open(f5) as fh: #open tsv file, set to fh
    for line in fh:
        col_values = line.strip().split("\t") # strip and split by tab for each value in both columns
        index_reads = (col_values[4]) #grab index reads 
        list_index_reads.append(index_reads) # put into list_index_reads array 
    #print(list_index_reads)
    
# Generate all possible combinations of 2 index values and store in dict with value of 0
for currentLine in list_index_reads:
    for index in list_index_reads: #call all indexes in array
        if currentLine == index:
            indexComboDict[currentLine + "_" + index]=0 # set indexes to 0 [index1, index2:0]
        else:
            unmatchedIndexDict[currentLine + "_" + index]=0


######################## Read in four files and assign corresponding lines to variables ##########################
# Initialize variables

headerR1 = ""
headerR2 = ""
headerR3 = ""
headerR4 = ""

seqR1 = ""
seqR2 = ""
seqR3 = ""
seqR4 = ""

plusR1 = ""
plusR2 = ""
plusR3 = ""
plusR4 = ""

phredR1 = ""
phredR2 = ""
phredR3 = ""
phredR4 = ""

N_cntr = 0
matched_cntr = 0
indexHopped_cntr = 0
sequencing_error_cntr = 0
retained_cntr = 0
discarded_cntr = 0
total_cntr = 0

def rc(dna):
    seq_dict = {'A':'T','T':'A','G':'C','C':'G', 'N':'N'}
    return "".join([seq_dict[base] for base in reversed(dna)])
 
with open(f1, "rt") as fh_R1, open(f2, "rt") as fh_R2, open(f3) as fh_R3, open(f4) as fh_R4:
    for line in fh_R1:
        line = line.strip("\n")
        headerR1 = line.strip("\n")
        seqR1 = fh_R1.readline().strip("\n")
        plusR1 = fh_R1.readline().strip("\n")
        phredR1 = fh_R1.readline().strip("\n")

        headerR2 = fh_R2.readline().strip("\n")
        seqR2 = fh_R2.readline().strip("\n")
        plusR2 = fh_R2.readline().strip("\n")
        phredR2 = fh_R2.readline().strip("\n")

        headerR3 = fh_R3.readline().strip("\n")
        seqR3 = rc(fh_R3.readline().strip("\n"))
        plusR3 = fh_R3.readline().strip("\n")
        phredR3 = fh_R3.readline().strip("\n")

        headerR4 = fh_R4.readline().strip("\n")
        seqR4 = fh_R4.readline().strip("\n")
        plusR4 = fh_R4.readline().strip("\n")
        phredR4 = fh_R4.readline().strip("\n")

        val = seqR2 + "_" + seqR3

        total_cntr+=1
        
        if total_cntr%500000000:
            print("On line: ", total_cntr, flush = True)

        score = []
        for base in range(len(phredR2)):        
            score.append(convert_phred(phredR2[base]))
            score.append(convert_phred(phredR3[base]))
        
        if all(items >= c for items in score):  
            retained_cntr+=1 
            if val in indexComboDict:           
                indexComboDict[val] += 1
                matched_cntr += 1
            
            if val in unmatchedIndexDict:    
                unmatchedIndexDict[val] += 1
                indexHopped_cntr += 1
            
            if "N" in val:
                if val in nIndexDict:
                    nIndexDict[val]+=1
                    N_cntr += 1
                else:
                    nIndexDict[val]=1
                    N_cntr += 1
            
            if val not in indexComboDict and val not in unmatchedIndexDict and val not in nIndexDict:
                sequencing_error_cntr += 1      
    
        else:
            discarded_cntr += 1 
            
                
with open("1294_S1_L008_R"+ "_Output_table_Fixed_Raw-2", 'w') as output1:
    j = 0
    k = 0   
    m = 0
    output1.write("\n" + "Number of Undetermined N Reads: " + str(N_cntr))
    output1.write("\n" + "Number of Expected Paired Index Reads : " + str(matched_cntr))
    output1.write("\n" + "Number of Index Swapped Reads: " + str(indexHopped_cntr))
    output1.write("\n" + "Number of Sequencing Error Reads: " + str(sequencing_error_cntr))
    output1.write("\n" + "Number of Retained Reads: " + str(retained_cntr))
    output1.write("\n" + "Number of Discarded Reads: " + str(discarded_cntr))
    output1.write("\n" + "Total Number of Reads: " + str(total_cntr))
    output1.write("\n" + "Percentage of Reads Index Swapped: " + str((indexHopped_cntr/total_cntr)*100))
    output1.write("\n" + "Percentage of Expected Paired Index Reads: " + str((matched_cntr/total_cntr)*100))
    output1.write("\n" + "Cutoff: " + str(c) + "\n")
    output1.write("\n")
    output1.write("Dual Indexed Pairs" + "\t" + "Counts" + "\n")
    for i in indexComboDict:
        output1.write(str(i) + "\t" + str(indexComboDict[i]) + "\n")
        j+=1       
    output1.write("\n")
    output1.write(("Index Hopped Pairs" + "\t" + "Counts" + "\n"))
    for i in unmatchedIndexDict:
        output1.write(str(i) + "\t" + str(unmatchedIndexDict[i]) + "\n")
        k+=1
    output1.write("\n")
    output1.write(("N Indexed Pairs" + "\t" + "Counts" + "\n"))
    for i in nIndexDict:
        output1.write(str(i) + "\t" + str(nIndexDict[i]) + "\n")
        m+=1       
   

# N_cntr = 0
# matched_cntr = 0
# indexHopped_cntr = 0
# sequencing_error_cntr = 0
# retained_cntr = 0
# discarded_cntr = 0
# total_cntr = 0


