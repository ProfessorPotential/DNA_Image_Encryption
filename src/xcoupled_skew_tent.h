/*
 REFERENCES:  
 [1] R. A. Elmanfaloty and E. Abou-Bakr, “Random property enhancement of
    a 1D chaotic PRNG with finite precision implementation,” Chaos Solitons
    Fractals, vol. 118, pp. 134–144, 2019.
*/

#ifndef XCOUPLED_SKEW_TENT_H
#define XCOUPLED_SKEW_TENT_H

#include <stdint.h>
#include <iostream>
#include <iomanip>

using namespace std;

/* Helper MACRO to print 8-bit key generate*/
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

/* Definitions for cross-coupled skew tent map PRNG */
#define NUM_DATA_BIT 56
#define NUM_FRAC_BIT 48
#define MASK_DATA 0xFFFFFFFFFFFFFF
#define MASK_FRAC 0xFFFFFFFFFFFF
#define NUM_KEY_BIT 8
#define XOR_GROUP_SIZE NUM_FRAC_BIT / NUM_KEY_BIT

#define ONE_PT_ZERO (1ULL << NUM_FRAC_BIT)

/* Main function */
void generate_sequence_xcoupled_skew_tent(uint64_t x0, uint64_t y0, uint64_t p1, uint64_t p2, uint32_t num_seq, uint8_t* sequence);
void generate_permutation_xcoupled_skew_tent(uint64_t x0, uint64_t y0, uint64_t p1, uint64_t p2, uint32_t num_seq, uint32_t* permutation_map);



/* All other functions defined */
void pre_calculate_multipliers(uint64_t p, uint64_t* a1, uint64_t* a2); 
uint64_t skew_tent(uint64_t xn, uint64_t p, uint64_t a1, uint64_t a2);
uint8_t cross_couple(uint64_t xn, uint64_t yn);
uint8_t shift_xor(uint64_t x);
uint64_t subtract_fixed_pt(uint64_t a, uint64_t b);
uint64_t multiply_fixed_pt(uint64_t a, uint64_t b);
uint64_t divide_fixed_pt(uint64_t a, uint64_t b);
long double conv_double(uint64_t val); // Helper function to print the value in long double

#endif
