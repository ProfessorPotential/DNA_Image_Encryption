#include "DNA.h"


// following three functions maybe more resource-saving for HLS
void DNA_encode(unsigned char* data, unsigned char* output_data, int size, int rule){
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

void DNA_Encryption(unsigned char* img, unsigned char* eimage, int size, int r){
    //A = 00 | C = 01 | G = 10 | T = 11
    for (int i = 0; i < size; ++i) {
        for (int j = 6; j >= 0; j-=2) {    
            switch (r) {
                case 1: //Rule One  Encoding
                    eimage[i] = img[i];
                    eimage[i+1] = img[i+1];
                    break;
                case 2: //Rule Two Encoding
                    //Modifying bits for ('00') Case
                    if(((img[i] >> (j+1)) & 0x01) == 0 && ((img[i] >> j) & 0x01) == 0){
                        eimage[i] &= ~(1 << (j+1)); //setting to 0
                        eimage[i] &= ~(1 << j); //setting to 0
                    }
                    //Modifying bits for ('01') Case
                    if(((img[i] >> (j+1)) & 0x01) == 0 && ((img[i] >> j) & 0x01) == 1){
                        eimage[i] |= (1 << (j+1)); //setting to 1
                        eimage[i] &= ~(1 << j); //setting to 0
                    }
                    //Modifying bits for ('10') Case
                    if(((img[i] >> (j+1)) & 0x01) == 1 && ((img[i] >> j) & 0x01) == 0){
                        eimage[i] &= ~(1 << (j+1)); //setting to 0
                        eimage[i] |= (1 << j); //setting to 1
                    } 
                    //Modifying bits for ('11') Case
                    if(((img[i] >> (j+1)) & 0x01) == 1 && ((img[i] >> j) & 0x01) == 1){
                        eimage[i] |= (1 << (j+1)); //setting to 1
                        eimage[i] |= (1 << j); //setting to 1
                    }
                    break;
                case 3: //Rule Three Encoding
                    //Modifying bits for ('00') Case
                    if(((img[i] >> (j+1)) & 0x01) == 0 && ((img[i] >> j) & 0x01) == 0){
                        eimage[i] |= (1 << (j+1)); //setting to 1
                        eimage[i] |= (1 << j); //setting to 1
                    }
                    //Modifying bits for ('01') Case
                    if(((img[i] >> (j+1)) & 0x01) == 0 && ((img[i] >> j) & 0x01) == 1){
                        eimage[i] &= ~(1 << (j+1)); //setting to 0
                        eimage[i] |= (1 << j); //setting to 1
                    }
                    //Modifying bits for ('10') Case
                    if(((img[i] >> (j+1)) & 0x01) == 1 && ((img[i] >> j) & 0x01) == 0){
                        eimage[i] |= (1 << (j+1)); //setting to 1
                        eimage[i] &= ~(1 << j); //setting to 0
                    } 
                    //Modifying bits for ('11') Case
                    if(((img[i] >> (j+1)) & 0x01) == 1 && ((img[i] >> j) & 0x01) == 1){
                        eimage[i] &= ~(1 << (j+1)); //setting to 0
                        eimage[i] &= ~(1 << j); //setting to 0
                    }
                    break;
                case 4: //Rule Four Encoding
                    //Modifying bits for ('00') Case
                    if(((img[i] >> (j+1)) & 0x01) == 0 && ((img[i] >> j) & 0x01) == 0){
                        eimage[i] |= (1 << (j+1)); //setting to 1
                        eimage[i] |= (1 << j); //setting to 1
                    }
                    //Modifying bits for ('01') Case
                    if(((img[i] >> (j+1)) & 0x01) == 0 && ((img[i] >> j) & 0x01) == 1){
                        eimage[i] |= (1 << (j+1)); //setting to 1
                        eimage[i] &= ~(1 << j); //setting to 0
                    }
                    //Modifying bits for ('10') Case
                    if(((img[i] >> (j+1)) & 0x01) == 1 && ((img[i] >> j) & 0x01) == 0){
                        eimage[i] &= ~(1 << (j+1)); //setting to 0
                        eimage[i] |= (1 << j); //setting to 1
                    } 
                    //Modifying bits for ('11') Case
                    if(((img[i] >> (j+1)) & 0x01) == 1 && ((img[i] >> j) & 0x01) == 1){
                        eimage[i] &= ~(1 << (j+1)); //setting to 0
                        eimage[i] &= ~(1 << j); //setting to 0
                    }
                    break;
                case 5: //Rule five encoding
                    //Modifying bits for ('00') Case
                    if(((img[i] >> (j+1)) & 0x01) == 0 && ((img[i] >> j) & 0x01) == 0){
                        eimage[i] &= ~(1 << (j+1)); //setting to 0
                        eimage[i] |= (1 << j); //setting to 1
                    }
                    //Modifying bits for ('01') Case
                    if(((img[i] >> (j+1)) & 0x01) == 0 && ((img[i] >> j) & 0x01) == 1){
                        eimage[i] &= ~(1 << (j+1)); //setting to 0
                        eimage[i] &= ~(1 << j); //setting to 0
                    }
                    //Modifying bits for ('10') Case
                    if(((img[i] >> (j+1)) & 0x01) == 1 && ((img[i] >> j) & 0x01) == 0){
                        eimage[i] |= (1 << (j+1)); //setting to 1
                        eimage[i] |= (1 << j); //setting to 1
                    } 
                    //Modifying bits for ('11') Case
                    if(((img[i] >> (j+1)) & 0x01) == 1 && ((img[i] >> j) & 0x01) == 1){
                        eimage[i] |= (1 << (j+1)); //setting to 1
                        eimage[i] &= ~(1 << j); //setting to 0
                    }
                    break;
                case 6: //Rule six encoding
                    //Modifying bits for ('00') Case
                    if(((img[i] >> (j+1)) & 0x01) == 0 && ((img[i] >> j) & 0x01) == 0){
                        eimage[i] &= ~(1 << (j+1)); //setting to 0
                        eimage[i] |= (1 << j); //setting to 1
                    }
                    //Modifying bits for ('01') Case
                    if(((img[i] >> (j+1)) & 0x01) == 0 && ((img[i] >> j) & 0x01) == 1){
                        eimage[i] |= (1 << (j+1)); //setting to 1
                        eimage[i] |= (1 << j); //setting to 1
                    }
                    //Modifying bits for ('10') Case
                    if(((img[i] >> (j+1)) & 0x01) == 1 && ((img[i] >> j) & 0x01) == 0){
                        eimage[i] &= ~(1 << (j+1)); //setting to 0
                        eimage[i] &= ~(1 << j); //setting to 0
                    } 
                    //Modifying bits for ('11') Case
                    if(((img[i] >> (j+1)) & 0x01) == 1 && ((img[i] >> j) & 0x01) == 1){
                        eimage[i] |= (1 << (j+1)); //setting to 1
                        eimage[i] &= ~(1 << j); //setting to 0
                    }
                    break;
                case 7:
                    //Modifying bits for ('00') Case
                    if(((img[i] >> (j+1)) & 0x01) == 0 && ((img[i] >> j) & 0x01) == 0){
                        eimage[i] |= (1 << (j+1)); //setting to 1
                        eimage[i] &= ~(1 << j); //setting to 0
                    }
                    //Modifying bits for ('01') Case
                    if(((img[i] >> (j+1)) & 0x01) == 0 && ((img[i] >> j) & 0x01) == 1){
                        eimage[i] &= ~(1 << (j+1)); //setting to 0
                        eimage[i] &= ~(1 << j); //setting to 0
                    }
                    //Modifying bits for ('10') Case
                    if(((img[i] >> (j+1)) & 0x01) == 1 && ((img[i] >> j) & 0x01) == 0){
                        eimage[i] |= (1 << (j+1)); //setting to 1
                        eimage[i] |= (1 << j); //setting to 1
                    } 
                    //Modifying bits for ('11') Case
                    if(((img[i] >> (j+1)) & 0x01) == 1 && ((img[i] >> j) & 0x01) == 1){
                        eimage[i] &= ~(1 << (j+1)); //setting to 0
                        eimage[i] |= (1 << j); //setting to 1
                    }
                    break;
                case 8:
                    //Modifying bits for ('00') Case
                    if(((img[i] >> (j+1)) & 0x01) == 0 && ((img[i] >> j) & 0x01) == 0){
                        eimage[i] |= (1 << (j+1)); //setting to 1
                        eimage[i] &= ~(1 << j); //setting to 0
                    }
                    //Modifying bits for ('01') Case
                    if(((img[i] >> (j+1)) & 0x01) == 0 && ((img[i] >> j) & 0x01) == 1){
                        eimage[i] |= (1 << (j+1)); //setting to 1
                        eimage[i] |= (1 << j); //setting to 1
                    }
                    //Modifying bits for ('10') Case
                    if(((img[i] >> (j+1)) & 0x01) == 1 && ((img[i] >> j) & 0x01) == 0){
                        eimage[i] &= ~(1 << (j+1)); //setting to 0
                        eimage[i] &= ~(1 << j); //setting to 0
                    } 
                    //Modifying bits for ('11') Case
                    if(((img[i] >> (j+1)) & 0x01) == 1 && ((img[i] >> j) & 0x01) == 1){
                        eimage[i] &= ~(1 << (j+1)); //setting to 0
                        eimage[i] |= (1 << j); //setting to 1
                    }
                    break;
                default:
                    eimage[i] = img[i];
                    eimage[i+1] = img[i+1];
                    break;
            } 
        }
    }  
    return; 
}

