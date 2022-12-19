#include<iostream>
#include "include/cm.h"
// #include "include/hash.h"

using namespace std;
int main(){
    ccmbase <uint64_t, WangHash> *sketch;
    sketch = new ccmbase<uint64_t, WangHash>(20, 7, 137, 0);


    sketch->update_count(100);
}

