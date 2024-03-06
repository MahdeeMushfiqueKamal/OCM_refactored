#include <vector>
#include <random>
#include <limits>
#include <fstream>
#include <cmath>
#include <assert.h>
#include "hash.h"
#include "support.h"

using std::allocator;

template<typename CounterType=int32_t , typename HashStruct = WangHash >
class ocmbase{
    std::vector<CounterType, allocator<CounterType>> core_; //resisters of the hash table
    std::vector<unsigned int> collision_;  // will keep track of collision after each round
    uint32_t np_;  // no of column (W) is 2^np_
    uint32_t nh_;  // number of hash functions
    uint64_t mask_;   // and-ing a number with mask will give X mod W
    uint64_t seedseed_;
    const bool conservative_;
    HashStruct hf_;
    std::vector<uint64_t, allocator<uint64_t>> seeds_;

public:
    ocmbase(unsigned np, unsigned nh=10, unsigned seedseed=137, bool conservative = false):
        np_(np),nh_(nh),mask_((1ull << np_) - 1),seedseed_(seedseed),conservative_(conservative)
        {
            core_.resize(nh_ << np_);
            collision_.resize(nh_ << np_);
            for(int i=0;i<collision_.size(); i++) collision_[i] = 0;
            assert(core_.size() == (nh_ << np_) );
            std::mt19937_64 mt(seedseed_ + 4);
            while(seeds_.size() < static_cast<unsigned>(nh_)) seeds_.emplace_back(mt());
        }

    void clear_core(){
        core_.clear();
        core_.shrink_to_fit();
        core_.resize(nh_ << np_);
    }

    void clear_collision(){
        collision_.clear();
        collision_.shrink_to_fit();
        collision_.resize(nh_ << np_);
    }

    void update_count(uint64_t val, int round, int total_round) {  
        int min_collision = std::numeric_limits<int>::max();
        std::vector<uint64_t> pos(nh_);
        for(unsigned added = 0; added < nh_; added++){
            CounterType hv = hf_(val ^ seeds_[added]);
            pos[added] = (hv & mask_) + (added << np_);   // exact positions where we will increase the counter by one.
            min_collision = std::min(min_collision, (int)collision_[pos[added]]);
        }

        for(unsigned added = 0; added < nh_; added++){
            if( collision_[pos[added]] == min_collision) core_[pos[added]]++;
        }
    }

    void update_collision(uint64_t val, int round, int total_round){
        int min_collision = std::numeric_limits<int>::max();
        std::vector<uint64_t> pos(nh_);
        for(unsigned added = 0; added < nh_; added++){
            CounterType hv = hf_(val ^ seeds_[added]);
            pos[added] = (hv & mask_) + (added << np_);   // exact positions where we will increase the counter by one.
            min_collision = std::min(min_collision, (int)collision_[pos[added]]);
        }


        if (min_collision >= round-2){
            // find min-count
            CounterType min_count = std::numeric_limits<CounterType>::max();
            for(unsigned added = 0; added < nh_; added++){
                min_count = (std::min<CounterType>)(core_[pos[added]], min_count);
            }

            for(unsigned added = 0; added < nh_; added++){
                if (core_[pos[added]] > min_count){
                    collision_[pos[added]] = round - 1;
                }
            }
        }
    }

    void update_count_collision(uint64_t val, int round, int total_round){
        int min_collision = std::numeric_limits<int>::max();
        std::vector<uint64_t> pos(nh_);
        for(unsigned added = 0; added < nh_; added++){
            CounterType hv = hf_(val ^ seeds_[added]);
            pos[added] = (hv & mask_) + (added << np_);   // exact positions where we will increase the counter by one.
            min_collision = std::min(min_collision, (int)collision_[pos[added]]);
        }

        if(min_collision < round - 1){
            // #>=1 cell without collision in prev round
            CounterType min_count = std::numeric_limits<int>::max();
            for(unsigned added=0; added< nh_; added++){
                if(collision_[pos[added]] == min_collision) min_count = std::min(min_count, core_[pos[added]]);
            }

            for(unsigned added=0; added< nh_; added++){
                if(collision_[pos[added]] == min_collision){
                    // c[i,j] = min(C[i,j]+1, min_count)        Changed code from paper
                    if(core_[pos[added]] == min_count){
                        core_[pos[added]] = min_count+ 1;
                    }
                }
            }
        }

        else{
            // every cell has a collision in the prev round
            CounterType min_count = std::numeric_limits<int>::max();
            for(unsigned added=0; added< nh_; added++){
                min_count = std::min(min_count, core_[pos[added]]);
            }

            for(unsigned added=0; added< nh_; added++){
                if(round < total_round && core_[pos[added]] > min_count){
                    collision_[pos[added]] = round;
                }
                // Changed Code from paper
                //core_[pos[added]] = std::min( core_[pos[added]] + 1, min_count);
                if(core_[pos[added]] == min_count){
                    core_[pos[added]] = min_count + 1;
                }
            }
        }
    }

    CounterType estimate_count(uint64_t val) const {
        int min_collision = std::numeric_limits<int>::max();
        std::vector<uint64_t> pos(nh_);
        // map kmers and find min collision
        for(unsigned added = 0; added < nh_; added++){
            CounterType hv = hf_(val ^ seeds_[added]);
            pos[added] = (hv & mask_) + (added << np_);   // exact positions where we will increase the counter by one.
            min_collision = std::min(min_collision, (int)collision_[pos[added]]);
        }

        CounterType min_count = std::numeric_limits<CounterType>::max();

        // min count with smallest collision number
        for(unsigned added = 0; added < nh_; added++){
            if( collision_[pos[added]] == min_collision) min_count = std::min(min_count,core_[pos[added]]);
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
            std::cout<<"core is saved------\n";
            for(unsigned i = 0; i < collision_.size(); i++){
                int temp = collision_[i];
                outputfile.write(reinterpret_cast<char*>(&temp), sizeof(temp));
                //outputfile.write(reinterpret_cast<char*>(&collision_[i]), sizeof(collision_[i]));
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
            input_file.seekg(sizeof(uint32_t)*2 + sizeof(uint64_t), std::ios::beg);            //future work
            for(uint32_t i=0; i< (nh_<<np_) ; i++){
                input_file.read(reinterpret_cast<char *>(&core_[i]), sizeof(core_[i]));
            }
            std::cout<<"core is read from sketch file------"<<std::endl;
            for(uint32_t i=0; i< (nh_<<np_); i++){
                int temp;
                input_file.read(reinterpret_cast<char *>(&temp), sizeof(temp));
                collision_[i]=temp;
            }
        }

    }
};


template<typename CounterType=int32_t , typename HashStruct = WangHash >
class OfflineCountMinSketch{
    ocmbase<CounterType, HashStruct> *sketch;
    int total_round;
    bool conservative;
    bool canonicalize;
    int seed = 137;
    // the function pointer for all the update functions
    typedef void (ocmbase<CounterType, HashStruct>::*UpdateFunctionPointer)(uint64_t, int, int);
    UpdateFunctionPointer updateFunctionPointer;
public:
    void createOfflineCountMinSketch(int np, int nh, int total_round, bool conservative, bool canonicalize){
        sketch = new ocmbase<CounterType, HashStruct>(np, nh, seed, conservative);
        this->total_round = total_round;
        this->conservative = conservative;
        this->canonicalize = canonicalize;
        if(sketch != NULL){
            if(conservative) std::cout<<"New Offline Conservative Count Min Sketch Created with height: "<< nh << " width: "<< (1<<np) << "\n";
            else std::cout<<"New Offline Count Min Sketch Created\n";
        }
    }

    void createOfflineCountMinSketch(std::string input_sketch_name){
        if(!sketch) delete sketch;
        ifstream input_sketch_file(input_sketch_name, std::ios::in | std::ios::binary);
        if(!input_sketch_file.is_open()){
            cout<<"Can't open sketch file\n";
            return;
        }

        uint32_t NP, NH;  uint64_t SEED;
        input_sketch_file.read(reinterpret_cast<char *>(&NP), sizeof(NP));
        input_sketch_file.read(reinterpret_cast<char *>(&NH), sizeof(NH));
        input_sketch_file.read(reinterpret_cast<char *>(&SEED), sizeof(SEED));
        input_sketch_file.close();

        sketch = new ocmbase<CounterType, HashStruct>(NP, NH, SEED, false);

        // note while creating a new object from sketch file, we are setting conservative = False. 
        // This is because, we are not storing the conservative flag in the sketch file.
        sketch->load_from_sketch(input_sketch_name);
    }

    void updateFromFile(std::string filename, int kmerLen, int currentRound){
        std::ifstream fasta_file(filename);
        int64_t currentKmer = 0; int current_len = 0;
        int chunk_size = 1000;
        char arr_chunk[chunk_size];
        bool isInHeader = false;
        uint64_t MASK = (1ull << (2*kmerLen)) -1;
        uint64_t checksum = 0;
        // std::cout<<"value of mask: "<<MASK<<std::endl;
        // int chunk_count = 0; int line_no =1;
        while(!fasta_file.eof())
        {
            // chunk_count++;
            fasta_file.read(arr_chunk, chunk_size);
            std::streamsize originalChunkSize = fasta_file.gcount();
            // https://stackoverflow.com/a/20911639/10941480
            // std::cout<<"Chunk found datasize: "<<dataSize<<"\n";
            // int i = 0;

            // while(i<chunk_size)
            for(int i=0; i< originalChunkSize; i++){
                char ch = arr_chunk[i];         
                // if (ch == '\n')line_no++;
                // std::cout<<"char found "<<ch<<std::endl;
                if(ch == '>'){
                    isInHeader = true;
                    currentKmer = 0; current_len = 0;
                    continue;
                }
                else if (isInHeader==true && ch=='\n'){
                    isInHeader = false;
                    continue;
                }
                if(isInHeader){continue;}
                if(ch=='\n' || ch=='\r' || ch==' '){continue;}
                if(ch=='N' || ch == 'n' || ch == 'X' || ch == 'x'){
                    currentKmer = 0;current_len = 0;
                    continue;
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
                        checksum += currentKmer;
                        // call the function pointer
                        (sketch->*updateFunctionPointer)(currentKmer, currentRound, total_round);
                        if(canonicalize) (sketch->*updateFunctionPointer)(reverse_compliment(currentKmer, kmerLen), currentRound, total_round);
                    }
                }
            }
        }
        std::cout<<"Updated Kmers from "<<filename<<", checksum "<<checksum<<"\n";
    }

    void constructOfflineCountMinSketch(std::string input_file, int kmerLen){
        sketch->clear_core();
        sketch->clear_collision();

        for(int current_round = 1; current_round <= total_round; current_round++){
            std::cout << "Performing Round : "<< current_round <<"\n";
            if(current_round > 1){
                std::cout<< "Will now perform update_collision\n";
                updateFunctionPointer = &ocmbase<CounterType, HashStruct>::update_collision;

                updateFromFile(input_file, kmerLen, current_round);
            }
            // clear counter array
            sketch->clear_core();

            if(conservative){
                std::cout<< "Will now perform update_count_collision\n";
                updateFunctionPointer = &ocmbase<CounterType, HashStruct>::update_count_collision;
            }
            else{
                std::cout<< "Will now perform update_count\n";
                updateFunctionPointer = &ocmbase<CounterType, HashStruct>::update_count;
            
            }
            updateFromFile(input_file, kmerLen, current_round);
        }

    }

    CounterType estimate_count(std::string kmer_str){
        uint_fast64_t kmer_val = cal(kmer_str);
        return sketch->estimate_count(kmer_val);
    }

    void save_sketch(std::string output_file_name){
        sketch->save_sketch(output_file_name);
    }
};


