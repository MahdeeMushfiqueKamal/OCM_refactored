# Offline (conservative) count-min sketch

A hash-based data structure for approximating k-mer count from a fasta file in (fixed) sublinear space. 
Performs better than count-min-sketch and conservative-count-min-sketch

## How to Run 
**Requirements:** 
- `build-essential` on linux. 
- `MinGW` on windows. 


**Compilation:**

Unix/Linux (using GCC):
```
#!/bin/bash
g++ -std=c++17 -o ocm ocm.cpp
```

Windows (using MinGW):
```
@echo off
g++ -std=c++17 -o ocm.exe ocm.cpp
```

macOS (using clang++, which comes pre-installed):
```
#!/bin/bash
clang++ -std=c++17 -o ocm ocm.cpp
```

**Count:**
Run `./ocm count`

Command Line Arguments: 

- `-k [value]`: Set the length of k-mers.
- `-h [value]`: Set the number of hash function / height for the Count Min Sketch. Default value = 7
- `-w [value]`: Set the width parameter for the Count Min Sketch. Default Value = 1048576
- `-n [value]`: Set the total number of rounds. Default value = 4
- `-o [file]`: Specify the output file for the sketch.
- `-fa [file]`: Specify the input fasta file containing DNA sequences.
- `-r:` Disable canonicalization of k-mers. Default value: Automatically canonicalize
- `-c:` Enable conservative update strategy. Default = False

Example: `./ocm count -k 22 -fa input/lhg22L20MC5x.fa -o sketch_file.sketch`

**Query:**
Run `./ocm query`

Command Line Arguments: 
- `-f [value]` : Input Sketch File. 
- `-q [value]` : Query File Name. Contains k-mers for which we want to report estimated counts. 
- `-o [value]` : Output File Name. Returns the k-mers and their estimated count in a csv format. 

Example: `./ocm query -f sketch_file.sketch -q input/test_exact_count_lhg22L20MC5x_20000.txt -o output/query_result.csv`

**Test:**
Run `./ocm test `

Command line arguments are similar to Query command. But, the query file contains k-mer and their true count. The output file reports k-mer, true_count, estimated_count.

Example: `./ocm query -f sketch_file.sketch -q input/test_exact_count_lhg22L20MC5x_20000.txt -o output/query_result.csv`

