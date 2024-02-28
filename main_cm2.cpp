#include<iostream>
#include"include/cm2.h"
#include <chrono>


using namespace std;

int main(){
    CountMinSketch2 <uint64_t, WangHash> countMinSketchObject;
    countMinSketchObject.createCountMinSketch(7,1048576*2, true, true);
    
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    countMinSketchObject.updateCountFromFile("input/lhg22L20MC5x.fa",22);
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

    std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << "[sec]" << std::endl;
    countMinSketchObject.save_sketch("sketch_file.sketch");


    countMinSketchObject.createCountMinSketch("sketch_file.sketch");

    string query_file_name="input/test_exact_count_lhg22L20MC5x_20000.txt", query_result_file_name="output/query_result.csv";
    ifstream infile(query_file_name);
    if(!infile.good()){cout<<"Couldn't open input query file\n";}
    ofstream query_result_file;
    query_result_file.open(query_result_file_name, ios::out);
    if(!query_result_file.good()){cout<<"Couldn't open output query file file\n";}
    query_result_file<<"kmer,true_count,estimated_count\n";

    string kmer; int true_count;
    while (infile >> kmer >> true_count){
        // cout<<kmer<<","<<true_count<<","<<countMinSketchObject.estimate_count(kmer)<<endl;
        query_result_file <<kmer<<","<<true_count<<","<<countMinSketchObject.estimate_count(kmer)<<endl;
    }

    
}