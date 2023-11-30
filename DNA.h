#include <stdint.h>
#include <iostream>

#ifndef DNA_H
#define DNA_H

#define NUM_RULES 8

#define NUM_NUCLEOTIDES 4
#define NUM_NUCLEOTIDES_IN_DATA 4
#define NUCLEOTIDE_MASK 0b11
#define NUCLEOTIDE_LEN 2

#define A_0 0b00
#define G_0 0b01
#define C_0 0b10
#define T_0 0b11

#define A_1 0b00
#define G_1 0b10
#define C_1 0b01
#define T_1 0b11

#define A_2 0b01
#define G_2 0b00
#define C_2 0b11
#define T_2 0b10

#define A_3 0b10
#define G_3 0b00
#define C_3 0b11
#define T_3 0b01

#define A_4 0b01
#define G_4 0b11
#define C_4 0b00
#define T_4 0b10

#define A_5 0b10
#define G_5 0b11
#define C_5 0b00
#define T_5 0b01

#define A_6 0b11
#define G_6 0b01
#define C_6 0b10
#define T_6 0b00

#define A_7 0b11
#define G_7 0b10
#define C_7 0b01
#define T_7 0b00


const uint8_t DNA_LUT[NUM_RULES][NUM_NUCLEOTIDES] = {
    {A_0, G_0, C_0, T_0},
    {A_1, G_1, C_1, T_1},
    {A_2, G_2, C_2, T_2},
    {A_3, G_3, C_3, T_3},
    {A_4, G_4, C_4, T_4},
    {A_5, G_5, C_5, T_5},
    {A_6, G_6, C_6, T_6},
    {A_7, G_7, C_7, T_7}
};

void DNA_encode(unsigned char* data, unsigned char* output_data, int size, int rule);
void DNA_add(unsigned char* data1, unsigned char* data2, unsigned char* output_data, int size);
void DNA_subtract(unsigned char* data1, unsigned char* data2, unsigned char* output_data, int size);
uint8_t DNA_add_nucleotides(uint8_t a, uint8_t b);
void DNA_Encryption(unsigned char* img, unsigned char* eimage, int size, int r);



#endif