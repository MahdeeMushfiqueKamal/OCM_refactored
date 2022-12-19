#include<iostream>
using namespace std;
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

int64_t reverse_compliment(uint64_t cal_kmer, int kmer_length){
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


uint64_t cal(string str_k_mer)
{
    const char* char_array = str_k_mer.c_str();
    uint64_t k_mer = 0;
    for(int j=0; char_array[j]!='\0'; j++)
    {
        switch(char_array[j])
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
    }
    //cout<<"cal value: "<<k_mer;
    //my_binary(k_mer);
    return k_mer;

}