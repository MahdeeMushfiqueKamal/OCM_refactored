#include<bits/stdc++.h>
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
};


