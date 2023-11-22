#include <stdint.h>
#include <iostream>
#include <iomanip>

/*
 REFERENCES:  
 [1] R. A. Elmanfaloty and E. Abou-Bakr, “Random property enhancement of
    a 1D chaotic PRNG with finite precision implementation,” Chaos Solitons
    Fractals, vol. 118, pp. 134–144, 2019.
 [2] 
*/
using namespace std;
#define BYTE_TO_BINARY_PATTERN "%c%c%c%c%c%c%c%c"
#define BYTE_TO_BINARY(byte)  \
  ((byte) & 0x80 ? '1' : '0'), \
  ((byte) & 0x40 ? '1' : '0'), \
  ((byte) & 0x20 ? '1' : '0'), \
  ((byte) & 0x10 ? '1' : '0'), \
  ((byte) & 0x08 ? '1' : '0'), \
  ((byte) & 0x04 ? '1' : '0'), \
  ((byte) & 0x02 ? '1' : '0'), \
  ((byte) & 0x01 ? '1' : '0') 

#define MASK_40BIT 0xFFFFFFFFFF
#define NUM_BIT_KEY 8
#define NUM_BIT_XOR 32
#define XOR_GROUP_SIZE 4

uint64_t skew_tent(uint64_t xn, uint64_t p, uint64_t a1, uint64_t a2);
uint8_t cross_couple(uint64_t xn, uint64_t yn);
uint8_t shift_xor(uint64_t x);
uint64_t multiply_fixed_pt(uint64_t a, uint64_t b);
long double conv_double(uint64_t val); // Helper function for troubleshooting


int main(){
    /* Values loaded per publication [1] */
    uint64_t xo = (0x0000b00000);
    uint64_t p1 = (0x0040000000);
    uint64_t a1 = (0x0400000000);
    uint64_t a2 = (0x0155555555);

    uint64_t yo = (0x0011000000);
    uint64_t p2 = (0x0059999999);
    uint64_t b1 = (0x02db6db6db);
    uint64_t b2 = (0x0189d89d89);

    /* 
     * Quick sanity checks.
     * This part is to be disgarded for HLS trnasfer. 
     */
    cout << "********************************************************" << endl;
    cout << "Initial parameters for cross coupled skew-tent chaos map" << endl;
    cout << setprecision(10);  // Note 2^-32 will have precisions of 10 decimal points
    cout << "xo = " << conv_double(xo) << endl;
    cout << "yo = " << conv_double(yo) << endl;
    cout << "a1 = " << conv_double(a1) << endl << "a2 = " << conv_double(a2) << endl; 
    cout << "b1 = " << conv_double(b1) << endl << "b2 = " << conv_double(b2) << endl; 
    cout << "p1 = " << conv_double(p1) << endl << "p2 = " << conv_double(p2) << endl;
    cout << "********************************************************" << endl << endl;
    /*
     * PRNG (pseudo-random number generation) process  
     */
    uint64_t xprev = xo;
    uint64_t xn, yn; 
    uint64_t yprev = yo;
    uint8_t key;
    for(int i=0; i<5; i++){
        yn = skew_tent(xprev, p2, a1, a2);
        xn = skew_tent(yprev, p1, b1, b2);
        key = cross_couple(xn, yn); 
        cout << "Key for index " << i << ": ";
        printf(BYTE_TO_BINARY_PATTERN "\n", BYTE_TO_BINARY(key));
        yprev = yn; 
        xprev = xn;
    }
    return 0;
}

uint64_t skew_tent(uint64_t xn, uint64_t p, uint64_t a1, uint64_t a2){
    /*
     * Single skew tent map function.
     * Same variable names chosen as the publication [1]
     */
    if(xn <= p) {
        return multiply_fixed_pt(xn, a1);
    }
    return multiply_fixed_pt((uint64_t) (~xn + 1ULL + (1ULL << 32)) & MASK_40BIT, a2);
}

uint8_t shift_xor(uint64_t x){
    uint8_t xor_result = 0;
    for(int i=0; i<NUM_BIT_XOR; i++){
        bool bit = 0;
        for(int j=0; j<XOR_GROUP_SIZE; j++){
            bit ^= (x >> (i + (j * NUM_BIT_KEY))) % 2;
        }
        xor_result |= (bit << i);
    }
    return xor_result;
}
uint8_t cross_couple(uint64_t xn, uint64_t yn){
    return shift_xor(xn) ^ shift_xor(yn);
}

uint64_t multiply_fixed_pt(uint64_t a, uint64_t b){
    return (uint64_t)(((__uint128_t) a *(__uint128_t)  b)>>32);

}

long double conv_double(uint64_t val){
    return (long double)((val) / (long double)(1ULL << 32));
}