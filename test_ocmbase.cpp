#include<iostream>
#include "include/ocm.h"


// for round r = 1, . . . , R do
// if r > 1 then
// for kmer k in kmers do
// UPDATE-COLLISION(k, r)
// clear counter array C
// for kmer k in kmers do
// UPDATE-COUNT(k)

using namespace std;
int main(){
    ocmbase <uint64_t, WangHash> *sketch;

    vector<int> v = {1,2,3,4,5,6,7,8,9,10};

    sketch = new ocmbase<uint64_t, WangHash>(20, 7, 137, 0);


    for (int round = 1; round <= 4; round++){
        for(int val : v){
            if (round > 1){
                sketch->update_collision(val, round, 4);
            }
            sketch->clear_core();
            sketch->update_count_collision(val, round, 4);
        }
    }

    for( int val : v){
        cout << "val: " << val << " estimated count: " << sketch->estimate_count(val) << endl;
    }
}

