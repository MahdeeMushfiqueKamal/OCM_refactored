#include<bits/stdc++.h>
#include <assert.h>
#include "hash.h"

//support function
uint64_t addChar(uint64_t k_mer, char ch){
    switch(ch)
    {
        case 'A':
                k_mer = k_mer<<2; // A=00
                break;
        case 'T':
                k_mer = k_mer<<2;
                k_mer = k_mer | 1;   //T=01
                break;
        case 'G':
                k_mer = k_mer<<2;  //G=10
                k_mer = k_mer | 2;
                break;
        case 'C':
                k_mer = k_mer<<2; //C=11
                k_mer = k_mer | 3;
                break;
    }
    return k_mer;
}

int64_t reverse_compliment(uint64_t cal_kmer, int kmer_length)
    {
        uint64_t k_mer = 0;
        uint64_t mask = 3;

        for(int i=0; i<kmer_length; i++)
        {
            switch(cal_kmer & mask)
            {

            case 0:
                k_mer = k_mer<<2;
                k_mer = k_mer | 1;   //A=00->T=01
                break;
            case 1:
                k_mer = k_mer<<2; //T=01->A=00
                break;
            case 2:
                k_mer = k_mer<<2; //G=10->C=11
                k_mer = k_mer | 3;
                break;
            case 3:
                k_mer = k_mer<<2;  //C=11->G=10
                k_mer = k_mer | 2;
                break;

            }
            cal_kmer=cal_kmer>>2;
        }
        return k_mer;
    }

using std::allocator;
template<typename CounterType=int32_t , typename HashStruct = WangHash >
class ccmbase{
    std::vector<CounterType, allocator<CounterType>> core_;     //resisters of the hash table
    uint32_t np_;                                               // no of column (W) is 2^np_
    uint32_t nh_;                                               // number of hash functions
    uint64_t mask_;                                             // and-ing a number with mask will give X mod W
    uint64_t seedseed_;
    const bool conservative_;
    HashStruct hf_;
    std::vector<uint64_t, allocator<uint64_t>> seeds_;

public:
    ccmbase(unsigned np, unsigned nh=10, unsigned seedseed=137, bool conservative = false):
        np_(np),nh_(nh),mask_((1ull << np_) - 1),seedseed_(seedseed),conservative_(conservative)
        {
            core_.resize(nh_ << np_);
            assert(core_.size() == (nh_ << np_) );
            std::mt19937_64 mt(seedseed_ + 4);
            while(seeds_.size() < static_cast<unsigned>(nh_)) seeds_.emplace_back(mt());
        }

    void update_count(uint64_t val) {
        std::vector<uint64_t> pos(nh_);
        for(unsigned added = 0; added < nh_; added++){
            CounterType hv = hf_(val ^ seeds_[added]);
            pos[added] = (hv & mask_) + (added << np_);   // exact positions where we will increase the counter by one.
        }

        if(conservative_ == false){
            for(unsigned added = 0; added < nh_; added++) core_[pos[added]]++;
        }
        else if(conservative_ == true){

            CounterType min_count = std::numeric_limits<CounterType>::max();
            for(unsigned added = 0; added < nh_; added++){
                min_count = (std::min<CounterType>)(core_[pos[added]], min_count);
            }

            for(unsigned added = 0; added < nh_; added++){
                if(core_[pos[added]] == min_count) core_[pos[added]]++;
            }
        }
    }

    CounterType estimate_count(uint64_t val) const {
        std::vector<uint64_t> pos(nh_);
        for(unsigned added = 0; added < nh_; added++){
            CounterType hv = hf_(val ^ seeds_[added]);
            pos[added] = (hv & mask_) + (added << np_);   // exact positions where we will increase the counter by one.
        }

        CounterType min_count = std::numeric_limits<CounterType>::max();
        for(int i=0; i< nh_; i++){
            min_count = (std::min<CounterType>)(core_[pos[i]], min_count);
        }
        return min_count;
    }

    void save_sketch(std::string output_file_name){
        std::cout<<"saving sketch------\n";
        std::ofstream outputfile;
        outputfile.open(output_file_name, std::ios::out | std::ios::binary);
        // Write Binary File //
        if(outputfile.is_open()){
            outputfile.write(reinterpret_cast<char*>(&np_), sizeof(np_));
            outputfile.write(reinterpret_cast<char*>(&nh_), sizeof(nh_));
            outputfile.write(reinterpret_cast<char*>(&seedseed_), sizeof(seedseed_));

            for(unsigned i = 0; i < core_.size(); i++){
                outputfile.write(reinterpret_cast<char*>(&core_[i]), sizeof(core_[i]));
            }
            outputfile.close();
        }
        std::cout<<"Sketch is saved------\n";
    }

    void load_from_sketch(std::string input_file_name){
        std::cout<<"loading data from sketch file------\n";
        std::ifstream input_file;
        input_file.open(input_file_name, std::ios::in | std::ios::binary);
        if(input_file.is_open()){
            input_file.seekg(sizeof(uint32_t)*2 + sizeof(uint64_t), std::ios::beg);
            for(uint32_t i=0; i< (nh_<<np_) ; i++){
                input_file.read(reinterpret_cast<char *>(&core_[i]), sizeof(core_[i]));
                //std::cout<<core_[i]<<std::endl;
            }
            std::cout<<"core is read from sketch file------"<<std::endl;
        }

    }

    void printInfo(){
        std::cout<<np_<<" "<<nh_<<" \n";
    }
};


template<typename CounterType=int32_t , typename HashStruct = WangHash >
class CountMinSketch{
    ccmbase<CounterType, HashStruct> *sketch;
    bool conservative; 
    bool canonicalize;
    int seed; 
public:
    void createCountMinSketch(int height, int width, bool conservative, bool canonicalize){
        int np = log2(width);
        sketch = new ccmbase<CounterType, HashStruct>(height, np, seed, conservative);
    }

    void updateCountFromFile(std::string filename, int kmerLen){
        std::ifstream fasta_file(filename);
        int64_t currentKmer = 0; int current_len = 0;
        int chunk_size = 1000;
        char arr_chunk[chunk_size];
        bool isInHeader = false;
        uint64_t MASK = (1ull << (2*kmerLen)) -1;
        // std::cout<<"value of mask: "<<MASK<<std::endl;
        // int chunk_count = 0; int line_no =1;
        while(!fasta_file.eof())
        {
            // chunk_count++;
            fasta_file.read(arr_chunk, chunk_size);
            int i = 0;

            while(i<chunk_size)
            {
                char ch = arr_chunk[i];
                if(ch==EOF)break;            
                // if (ch == '\n')line_no++;
                // std::cout<<"char found "<<ch<<std::endl;
                if(ch == '>'){
                    isInHeader = true;
                    currentKmer = 0; current_len = 0;
                    i++; continue;
                }
                else if (isInHeader==true && ch=='\n'){
                    isInHeader = false;
                    i++; continue;
                }
                if(isInHeader){i+=1;continue;}
                if(ch=='\n' || ch=='\r' || ch==' '){i+=1;continue;}
                if(ch=='N'){
                    currentKmer = 0;current_len = 0;
                    i+=1;continue;
                }
                else
                {
                    // kmer_count++;
                    if(current_len < kmerLen){
                        currentKmer = addChar(currentKmer, ch);
                        current_len++;
                    }
                    else{
                        currentKmer = addChar(currentKmer, ch) & MASK;
                    }
                    // std::cout<<currentKmer<<" "<<current_len<<std::endl;

                    if(current_len==kmerLen)
                    {
                        //GOT KMER -- do the necessary things
                        sketch->update_count(currentKmer);
                        if(canonicalize) sketch->update_count(reverse_compliment(currentKmer, kmerLen));
                    }
                }
                i+=1;
            }
        }
    }

    void printInfo(){
        sketch->printInfo();
    }
};


