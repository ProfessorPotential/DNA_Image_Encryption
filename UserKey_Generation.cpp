#include "UserKey_Generation.h"
#include <iostream>
#include <bitset>
#include <vector>
#include "include/ap_int.h"



using namespace std;

vector<ap_int<21>> GenerateKey() {
    vector<ap_int<21>> keyVector;
    //User Generated 168 Bit Key
    bitset<168> key("10101001100111111010101011011010011110111011001111111010101011011001110001101011011011001001001010111011100110000110110101010111111101010110100100110101111100011011011101000101011110");

    for (int i = 0; i < 168; i += 21) {
        std::bitset<21> bitsSubset;
        for (int j = 0; j < 21; ++j) {
            bitsSubset[j] = key[i + j];
        }
        keyVector.push_back(static_cast<ap_int<21>>(bitsSubset.to_ulong()));
    }

    return keyVector;
}
