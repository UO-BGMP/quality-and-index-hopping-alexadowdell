# quality-and-index-hopping-alexadowdell

1. Contains:

-  `Final_IndexHoppingProject.Rmd`: Final Rmarkdown compilation of Index Hopping Project
-  `Final_IndexHoppingProject_html.pdf`: Final rendered HTML compilation of Index Hopping Project exported as a PDF
-  `Final_IndexHoppingScript.html`: Final rendered HTML compilation of Index Hopping Project
-  `Final_IndexHoppingScript.pdf`: Final rendered PDF compilation of Index Hopping Project

2. `outputfiles/` contains: 
  - Part 1-1: Distributions for Avg Quality Score for each Base Position
    - R1 file: `R1_Part1_dist_output.tsv`
    - R2 file: `R2_Part1_dist_output.tsv`
    - R3 file: `R3_Part1_dist_output.tsv`
    - R4 file: `R4_Part1_dist_output.tsv`
  
  - Part 1-2: Frequency Distributions based on Avg Quality Score
    - R1 file: `R1_part2Plot.txt`
    - R2 file: `R2_part2Plot.txt`
    - R3 file: `R3_part2Plot.txt`
    - R4 file: `R4_part2Plot.txt`
    
  - Part 2: Demultiplexing
    - Demultiplexing script output (contains cutoff 30 statistics, expected index pair dictionary, index swapped pair dictionary, undetermined (N) index pair dictionary): `1294_S1_L008_R_Output_table_Fixed_30-2`
    - Frequency Distribution of Expected Index Pairs at Cutoff 30: `expectedOutputTableFixed30.tsv`
    - Frequency Distribution of Swapped Index Pairs at Cutoff 30: `swappedOutputTableFixed.tsv`
    - Frequency Distribution of Swapped Index Pairs at Cutoff 30 Reformatted for Heatmap plotting:`swappedOutputTableFixedHeatmapReformatted.tsv`
    - Frequency Distribution of Expected Index Pairs at Cutoff 0: `expectedOutputTableFixed0-2.tsv`
    - Frequency Distribution of Swapped Index Pairs at Cutoff 0: `swappedIndexesRaw2.tsv`


3. `scripts/` contains: 
  - Part 1-1:
    - Mean Distribution for R1 file: `part1Hist1.py`
    - Mean Distribution for R2 file: `part1Hist2.py`
    - Mean Distribution for R3 file: `part1Hist3.py`
    - Mean Distribution for R4 file: `part1Hist4.py`
  
  - Part 1-2: 
    - Frequency Distribution for R1, R2, R3, R4 files: `part2FreqDistFixed.py`
    
  - Part 2:
    - Demultiplexing Cutoff 30: `outputTableFixed30-2.py`
    - Demultiplexing Cutoff 0: `outputTableFixed0-2.py`
  
  
  
  
  
  
