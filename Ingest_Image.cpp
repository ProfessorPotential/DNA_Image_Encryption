#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#include "include/ap_int.h"
#include <iostream>
#include <bitset>
#include <vector>
//#include "UserKey_Generation.h"

using namespace std;

/*--------------------------------GLOBAL VARIABLES-----------------------*/

/*-------------------------------FUNCTIONS------------------------------*/

// Function to convert an integer to binary representation
string intToBinary(int value, int numBits) {
    return bitset<32>(value).to_string().substr(32 - numBits);
}

// Function to calculate the hamming distance. Takes in a matrixed bitplane (8 bits) and outputs a 168 bit hamming distance
vector<ap_int<21>> HammingDistance(vector<uint8_t> bplane, uint8_t channels){
     //Collector for Hamming Distance
    uint64_t dist = 0ULL;
    uint8_t xor_bits;
    vector<ap_int<21>> result_array(8); 

    //Length of the Bitplanes
    int length = bplane.size();

    for (int i = 0; i < length; i+=2){
        for (int c = 0; c < channels; c++){
            xor_bits = bplane[i+c] & bplane[i+c+channels];
            while(xor_bits){
                dist += xor_bits & 1;
                xor_bits >>= 1;
            }
        }
    }
    //Max-Value = 256*256*8 = 524288 = 10000000000000000000 (20 bit)
    for (int f = 0; f < 8; f++){
        result_array[f] = ap_int<21>(dist);
    }
    return result_array;
}

// Function to check the size of a bitset to verify it is 168 bits for the key
template <size_t N>
bool check_size_168(bitset<N>& bitsetNumber){
    if (bitsetNumber.size() == 168){
        return true;
    }
    else return false;
}

// Function to DNA Encode an image
vector<uint8_t> DNA_Encode_Bitplane(vector<uint8_t> img, int r){
    //A = 00 | G = 01 | C = 10 | T = 11
    vector<uint8_t> eimage = img;
    for (int i = 0; i < img.size(); ++i) {
        for (int j = 0; j < 8; j+=2) {    
            switch (r) {
                case 1:
                    eimage[i] = img[i];
                    eimage[i+1] = img[i+1];
                    break;
                case 2:
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
                case 3:
                    //Modifying bits for A ('00') Case
                    if((img[i] >> j) == 0 && (img[i] >> (j+1)) == 0){
                        eimage[i] &= ~(1 << j); //setting to 0
                        eimage[i] |= (1 << (j+1)); //setting to 1
                    }
                    //Modifying bits for G ('01') Case
                    if((img[i] >> j) == 0 && (img[i] >> (j+1)) == 1){
                        eimage[i] &= ~(1 << j); //setting to 0
                        eimage[i] &= ~(1 << (j+1)); //setting to 0
                    }
                    //Modifying bits for G ('10') Case
                    if((img[i] >> j) == 1 && (img[i] >> (j+1)) == 0){
                        eimage[i] |= (1 << j); //setting to 1
                        eimage[i] |= (1 << (j+1)); //setting to 1
                    } 
                    //Modifying bits for T ('11') Case
                    if((img[i] >> j) == 1 && (img[i] >> (j+1)) == 1){
                        eimage[i] |= (1 << j); //setting to 1
                        eimage[i] &= ~(1 << (j+1)); //setting to 0
                    }
                    break;
                case 4:
                    break;
                case 5:
                    break;
                case 6:
                    break;
                case 7:
                    break;
                case 8:
                    break;
                default:
                    break;
            } 
        }
    }   
    return eimage;
}

unsigned char* encoded_image(vector<uint8_t> img){
    //Initializing Variables
    vector<unsigned char> result(img.size(),0);

    //Translating 8-bit values to unsigned char
    for (int i = 0; i < img.size(); ++i) {
        result[i] = static_cast<unsigned char>(img[i]); //converting to integer
    }

    //Translating vector into unsigned char*
    unsigned char* mydata = result.data();
    return mydata;
}

// Function to calculate the bitplane from a vector of intetgers

int main() {
    /*------------------------------IMAGE INGESTATION--------------------------------*/
    // Read an image from file
    const char* imagePath = "test_image.jpeg";
    int width, height, channels;

    unsigned char* image = stbi_load(imagePath, &width, &height, &channels, 0);

    // Check if the image was successfully loaded
    if (!image) {
        cout << "Error: Could not read the image." << endl;
        return -1;
    }

    // Print image information
    cout << "Width: " << width << ", Height: " << height << ", Channels: " << channels << endl;
    int size = width*height*channels;
    /*-------------------------------------------------------------------------------*/


    //Modify the unsigned char variables into 8-bit (uint_8) variables
    vector<uint8_t> bplane(size);
    for (int i = 0; i < width * height * channels; ++i) {
        bplane[i] = image[i];
    }

    /*---------------------------TEST FUNCTION TO VISUALIZE----------------------------*/
    // Access pixel values and perform operations
    // For example, print the RGB values of the pixel at (x,y)
    int x = 25;
    int y = 25;

    // Set the 3rd bit (make it 1)
    //myVariable |= (1 << 2);

    // Clear the 3rd bit (make it 0)
    // myVariable &= ~(1 << 2);

    int index = (y * width + x) * channels;
    cout << "Pixel at (" << x << ", " << y << "): " << endl;
    for (int c = 0; c < channels; ++c) {
        cout << "Integer: " << static_cast<int>(image[index + c]) << " ";
        cout << "Binary: " << intToBinary(static_cast<int>(bplane[index + c]),8) << " ";
        cout << "Bit 0 is: " << ((bplane[index + c] >> 1) & 1) << endl;;
    }
    cout << endl;

    /*
    // Print the size of each bitplane
    for (int i = 0; i < 8; ++i) {
        cout << "Size of Bitplane " << i << ": " << bitplanes[i].size() << endl;
    }
    */

    /*-------------------------------------------------------------------------------*/
   
    // Calculate the Hamming Distance
    vector<ap_int<21>> HD = HammingDistance(bplane, channels);

    //DNA Encode the Original Matrix
    vector<uint8_t> encoded_bp = DNA_Encode_Bitplane(bplane,2);

    //Translating the encoded bit-plane back to unsigned char form
    unsigned char* e_img = encoded_image(encoded_bp);

    /*------------------------------KEY GENERATION-----------------------------------*/
    //Defining the user generated key
    vector<ap_int<21>> keyVector;
    bitset<168> key("10101001100111111010101011011010011110111011001111111010101011011001110001101011011011001001001010111011100110000110110101010111111101010110100100110101111100011011011101000101011110");
    int z = 0;
    for (int i = 0; i < 168; i += 21) {
        bitset<21> bitsSubset;
        for (int j = 0; j < 21; ++j) {
            if (j == 20){
                bitsSubset[j] = 0;
            }
            else
                bitsSubset[j] = key[i + j];
        }
        keyVector.push_back(ap_int<21>(bitsSubset.to_ulong()) ^ HD[z]);
        z += 1;
    }
    /*-------------------------------------------------------------------------------*/

    /*
    //cout << "168 Bit Hamming Distance:" << endl << HD << endl;
    //cout << "User Key: " << endl << key << endl;
    //cout << "Special Key:" << endl << skey << endl;
    */

    // Write Image for Debugging
    stbi_write_jpg("output.jpg", width, height, channels, e_img, 100);

    // Free the allocated memory for the image
    stbi_image_free(image);

    return 0;
}
