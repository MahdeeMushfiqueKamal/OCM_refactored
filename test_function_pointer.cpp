#include <iostream>
#include "include/hash.h"

template<typename CounterType=int32_t , typename HashStruct = WangHash >
class DummyClass {
public:
    void dummyFunction1() {
        std::cout << "Dummy function 1 called" << std::endl;
    }

    void dummyFunction2() {
        std::cout << "Dummy function 2 called" << std::endl;
    }

    void dummyFunction3() {
        std::cout << "Dummy function 3 called" << std::endl;
    }
};

template<typename CounterType=int32_t , typename HashStruct = WangHash >
class DummyCaller {
private:
    DummyClass<CounterType, HashStruct> dummyObj;
    void (DummyClass<CounterType, HashStruct>::*functionPointer)();

public:
    void callDummyFunctions() {
        functionPointer = &DummyClass<CounterType, HashStruct>::dummyFunction1;
        (dummyObj.*functionPointer)();
        
        functionPointer = &DummyClass<CounterType, HashStruct>::dummyFunction2;
        (dummyObj.*functionPointer)();

        functionPointer = &DummyClass<CounterType, HashStruct>::dummyFunction3;
        (dummyObj.*functionPointer)();
    }
};

int main() {
    DummyCaller<> caller; // You can specify template parameters here if needed
    caller.callDummyFunctions();
    return 0;
}
