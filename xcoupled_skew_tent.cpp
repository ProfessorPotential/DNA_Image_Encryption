#include "xcoupled_skew_tent.h"

// int main(){
//     /* Values loaded per publication [1] */
//     uint64_t xo = (0x0000b00000);
//     uint64_t p1 = (0x0040000000);
//     uint64_t a1 = (0x0400000000);
//     uint64_t a2 = (0x0155555555);

//     uint64_t yo = (0x0011000000);
//     uint64_t p2 = (0x0059999999);
//     uint64_t b1 = (0x02db6db6db);
//     uint64_t b2 = (0x0189d89d89);

//     /* 
//      * Quick sanity checks.
//      * This part is to be disgarded for HLS trnasfer. 
//      */
//     cout << "********************************************************" << endl;
//     cout << "Initial parameters for cross coupled skew-tent chaos map" << endl;
//     cout << setprecision(10);  // Note 2^-32 will have precisions of 10 decimal points
//     cout << "xo = " << conv_double(xo) << endl;
//     cout << "yo = " << conv_double(yo) << endl;
//     cout << "a1 = " << conv_double(a1) << endl << "a2 = " << conv_double(a2) << endl; 
//     cout << "b1 = " << conv_double(b1) << endl << "b2 = " << conv_double(b2) << endl; 
//     cout << "p1 = " << conv_double(p1) << endl << "p2 = " << conv_double(p2) << endl;
//     cout << "********************************************************" << endl << endl;
//     /*
//      * PRNG (pseudo-random number generation) process  
//      */
//     uint64_t xprev = xo;
//     uint64_t xn, yn; 
//     uint64_t yprev = yo;
//     uint8_t key;
//     for(int i=0; i<5; i++){
//         yn = skew_tent(xprev, p2, a1, a2);
//         xn = skew_tent(yprev, p1, b1, b2);
//         key = cross_couple(xn, yn); 
//         cout << "Key for index " << i << ": ";
//         printf(BYTE_TO_BINARY_PATTERN "\n", BYTE_TO_BINARY(key));
//         yprev = yn; 
//         xprev = xn;
//     }
//     return 0;
// }

void generate_keys_xcoupled_skew_tent(uint64_t x0, uint64_t y0, uint64_t p1, uint64_t p2, uint32_t num_keys, uint8_t* keys){
    uint64_t xprev = x0;
    uint64_t yprev = y0;
    uint8_t key;
    
    /*
     * Multiplication is faster than division. 
     * Pre-calculate divisions in terms of inverse 
     * multiplication form.
     * p1 -> a1 = 1/p1 & a2 = 1/(1-p1)
     * p2 -> b1 = 1/p2 & b2 = 1/(1-p2)
     */
    uint64_t xn, yn, a1, a2, b1, b2;
    a1 = (1ULL << NUM_FRAC_BIT * 2) / p1; // 1÷p1
    a2 = (1ULL << NUM_FRAC_BIT * 2) / subtract_fixed_pt(1ULL << NUM_FRAC_BIT, p1); // 1÷(1-p1)
    b1 = (1ULL << NUM_FRAC_BIT * 2) / p2; // 1÷p2
    b2 = (1ULL << NUM_FRAC_BIT * 2) / subtract_fixed_pt(1ULL << NUM_FRAC_BIT, p2); // 1÷(1-p2)

    /*
     * PRNG (pseudo-random number generation) process  
     */
    for(uint32_t i=0; i<num_keys; i++){
        xn = skew_tent(yprev, p2, b1, b2);
        yn = skew_tent(xprev, p1, a1, a2);
        key = cross_couple(xn, yn);
        xprev = xn;
        yprev = yn; 
        keys[i] = key;
    }
}

uint64_t skew_tent(uint64_t xn, uint64_t p, uint64_t a1, uint64_t a2){
    /*
     * Single skew tent map function.
     * Same variable names chosen as the publication [1]
     */
    if(xn <= p) {
        return multiply_fixed_pt(xn, a1);
    }
    return multiply_fixed_pt((uint64_t) (~xn + 1ULL + (1ULL << 32)) & MASK_DATA, a2);
}

uint8_t shift_xor(uint64_t x){
    uint8_t xor_result = 0;
    for(int i=0; i<NUM_KEY_BIT; i++){
        bool bit = 0;
        for(int j=0; j<XOR_GROUP_SIZE; j++){
            bit ^= (x >> (i + (j * NUM_KEY_BIT))) % 2;
        }
        xor_result |= (bit << i);
    }
    return xor_result;
}
uint8_t cross_couple(uint64_t xn, uint64_t yn){
    return shift_xor(xn) ^ shift_xor(yn);
}

uint64_t subtract_fixed_pt(uint64_t a, uint64_t b){
    return (a + (~b + 1ULL)) & MASK_DATA;// a - b = a + twos(b)
}

uint64_t multiply_fixed_pt(uint64_t a, uint64_t b){
    return (uint64_t)(((__uint128_t) a *(__uint128_t)  b)>>32);

}

long double conv_double(uint64_t val){
    return (long double)((val) / (long double)(1ULL << 32));
}