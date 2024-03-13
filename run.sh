#!bin/bash
rm ocm
g++ -std=c++17 -o ocm ocm.cpp
./ocm count -c -k 22 -h 7 -w 1048576 -n 4 -fa input/lhg22L20MC5x.fa -o sketch_file.sketch --bitsize 2
./ocm test -f sketch_file.sketch -q input/test_exact_count_lhg22L20MC5x_20000.txt -o output/query_result.csv --bitsize 2