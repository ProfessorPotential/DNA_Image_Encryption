#include "DNA.h"


// following three functions maybe more resource-saving for HLS
void DNA_code(unsigned char* data, unsigned char* output_data, int size, int rule){
    for(int i=0; i<size; i++){
        output_data[i] = 0;
        for(int j=0; j < NUM_NUCLEOTIDES_IN_DATA; j++){
            uint8_t nucleotide = (data[i] >> j * NUCLEOTIDE_LEN) & NUCLEOTIDE_MASK;
            output_data[i] |= (DNA_LUT[rule][nucleotide]) << j * NUCLEOTIDE_LEN;
        }
    }
}

void DNA_subtract(unsigned char* data1, unsigned char* data2, unsigned char* output_data, int size){
    for(int i=0; i<size; i++){
        output_data[i] = 0;
        for(int j=0; j < NUM_NUCLEOTIDES_IN_DATA; j++){
            uint8_t nucleotide_x = (data1[i] >> j * NUCLEOTIDE_LEN) & NUCLEOTIDE_MASK;
            uint8_t nucleotide_y = (data2[i] >> j * NUCLEOTIDE_LEN) & NUCLEOTIDE_MASK;
            output_data[i] |= (DNA_add_nucleotides(nucleotide_x, (~nucleotide_y)+1) & NUCLEOTIDE_MASK )<< (j * NUCLEOTIDE_LEN);
        }
    }
}


void DNA_add(unsigned char* data1, unsigned char* data2, unsigned char* output_data, int size){
    for(int i=0; i<size; i++){
        output_data[i] = 0;
        for(int j=0; j < NUM_NUCLEOTIDES_IN_DATA; j++){
            uint8_t nucleotide_x = (data1[i] >> j * NUCLEOTIDE_LEN) & NUCLEOTIDE_MASK;
            uint8_t nucleotide_y = (data2[i] >> j * NUCLEOTIDE_LEN) & NUCLEOTIDE_MASK;
            output_data[i] |= DNA_add_nucleotides(nucleotide_x, nucleotide_y) << (j * NUCLEOTIDE_LEN);
        }
    }
}

uint8_t DNA_add_nucleotides(uint8_t a, uint8_t b){
    return (a + b) & NUCLEOTIDE_MASK;
}
