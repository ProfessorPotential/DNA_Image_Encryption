#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#include "include/ap_int.h"
#include <iostream>
#include <bitset>
#include <vector>
using namespace std;

// Function to DNA Encode an image
unsigned char* DNA_Encryption(unsigned char* img, int size, int r){
    //A = 00 | C = 01 | G = 10 | T = 11
    unsigned char* eimage = new unsigned char[size];
    memset(eimage, 0, size); //initialize with 0s

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
    return eimage;
}

int main(){
    // Read an image from file
    const char* imagePath = "test_image.jpeg";
    int width, height, channels;

    unsigned char* image = stbi_load(imagePath, &width, &height, &channels, 0);
    int size = width*height*channels;
    unsigned char* r_image = new unsigned char[size];
    memset(r_image, 0, size); //initialize with 0s

    r_image = DNA_Encryption(image, size, 2);
    stbi_write_jpg("output.jpg", width, height, channels, r_image, 100);

    return 0;
}