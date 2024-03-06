#include<iostream>
#include"include/cm.h"
#include <chrono>

using namespace std;

int main(int argc, char *argv[]){
    string mode(argv[1]); //"count" or "query"
    int kmer_len = 0, height = 7, width = 1048576, total_round = 4;
    bool canonicalize = true, conservative = false;
    string output_sketch_file = "", input_fasta_file = "";
    if(mode == "count"){
        for(int i= 2; i<argc-1; i++){
            string arg(argv[i]);
            //string param(argv[++i]);
            if(arg=="-k") kmer_len=stoi(argv[++i]);
            else if(arg=="-h") height =stoi(argv[++i]);
            else if(arg=="-w") width = stoi(argv[++i]);
            else if(arg=="-o") output_sketch_file = argv[++i];
            else if(arg=="-fa") input_fasta_file=argv[++i];
            else if(arg=="-r") canonicalize = false;
            else if(arg=="-c") conservative = true;
        }

        if(kmer_len <= 0 || input_fasta_file=="" || output_sketch_file==""){
            cout<<"INVALID INPUT"; return 0;
        }

        CountMinSketch <uint64_t, WangHash> countMinSketchObject;
        // countMinSketchObject.createCountMinSketch(height, width, conservative, canonicalize);
        countMinSketchObject.createCountMinSketch(height, width, conservative, canonicalize);

        std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
        countMinSketchObject.updateCountFromFile(input_fasta_file,kmer_len);
        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

        std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << "[sec]" << std::endl;

        countMinSketchObject.save_sketch("sketch_file.sketch");
    }

    else if(mode == "query"){
        string input_sketch_name;
        string query_file_name;
        string query_result_file_name;

        for(int i=2; i<argc-1; i++)
        {
            string arg(argv[i]);
            string param(argv[++i]);
            if(arg=="-f") input_sketch_name=param;
            else if(arg=="-q") query_file_name=param;
            else if(arg=="-o") query_result_file_name=param;
        }
        cout<<"query"<<endl;

        CountMinSketch <uint64_t, WangHash> CountMinSketchObject;
        CountMinSketchObject.createCountMinSketch(input_sketch_name);

        ifstream infile(query_file_name);
        if(!infile.good()){cout<<"Couldn't open input query file\n";}

        ofstream query_result_file;
        query_result_file.open(query_result_file_name, ios::out);
        if(!query_result_file.good()){cout<<"Couldn't open output query file file\n";}
        query_result_file<<"kmer,estimated_count\n";

        string kmer;

        while (infile >> kmer){
            query_result_file <<kmer<<","<<CountMinSketchObject.estimate_count(kmer)<<endl;
        }
    }

    else if(mode == "test"){
        string input_sketch_name;
        string query_file_name;
        string query_result_file_name;

        for(int i=2; i<argc-1; i++)
        {
            string arg(argv[i]);
            string param(argv[++i]);
            if(arg=="-f") input_sketch_name=param;
            else if(arg=="-q") query_file_name=param;
            else if(arg=="-o") query_result_file_name=param;
        }
        cout<<"test"<<endl;

        CountMinSketch <uint64_t, WangHash> CountMinSketchObject;
        CountMinSketchObject.createCountMinSketch(input_sketch_name);

        ifstream infile(query_file_name);
        if(!infile.good()){cout<<"Couldn't open input query file\n";}

        ofstream query_result_file;
        query_result_file.open(query_result_file_name, ios::out);
        if(!query_result_file.good()){cout<<"Couldn't open output query file file\n";}
        query_result_file<<"kmer,true_count,estimated_count\n";

        string kmer; int true_count;

        while (infile >> kmer >> true_count){
            query_result_file <<kmer<<","<<true_count<<","<<CountMinSketchObject.estimate_count(kmer)<<endl;
        }
    }
}