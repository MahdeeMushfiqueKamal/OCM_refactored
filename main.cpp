#include<iostream>
#include"include/cm.h"


using namespace std;

int main(){
    CountMinSketch <uint64_t, WangHash> countMinSketchObject;
    countMinSketchObject.createCountMinSketch(7,1048576, false, true);
    countMinSketchObject.printInfo();
    countMinSketchObject.updateCountFromFile("input/lhg22L20MC5x.fa",22);
    // countMinSketchObject.updateCountFromFile("input/last.txt",22);
}