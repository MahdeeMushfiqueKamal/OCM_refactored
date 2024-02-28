#include<iostream>
#include "include/ocm.h"

using namespace std;
int main(){
    ocmbase <uint64_t, WangHash> *sketch;
    sketch = new ocmbase<uint64_t, WangHash>(20, 7, 137, 0);


    sketch->update_count(100,1,4);

    cout<< sketch->estimate_count(100);
}

