#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#include "ap_int.h"
#include <iostream>
#include <bitset>
#include <vector>

using namespace std;

/*--------------------------------GLOBAL VARIABLES-----------------------*/
//User Generated 168 Bit Key
bitset<168> key("10101001100111111010101011011010011110111011001111111010101011011001110001101011011011001001001010111011100110000110110101010111111101010110100100110101111100011011011101000101011110");

//Chaotic Parameters (using float since C++ treats these as 32-bit numbers instead of 64-bit like double data type)
bitset<40> x_0("0000000000000000101100000000000000000000"); //0000b00000 hexadecimal
bitset<40> y_0("0000000000010001000000000000000000000000"); //0011000000 hexadecimal
bitset<40> p_1("0000000001000000000000000000000000000000"); //0040000000 hexadecimal
bitset<40> p_2("0000000001011001100110011001100110011001"); //0059999999 hexadecimal

ap_uint<168> var;

/*-------------------------------FUNCTIONS------------------------------*/

// Function to convert an integer to binary representation
string intToBinary(int value, int numBits) {
    return bitset<32>(value).to_string().substr(32 - numBits);
}

// Function to calculate the hamming distance. Takes in a matrixed bitplane (8 bits) and outputs a 168 bit hamming distance
uint64_t HammingDistance(vector<uint8_t> bitplane) {
    //Collector for Hamming Distance
    uint64_t dist = 0;

    //Length of the Bitplanes
    int length = bitplane.size();
    uint64_t result;

    for (int i = 0; i < length; i++){
        for (int j = 0; j < 4; j++){
            int result = ((bitplane[i] >> j) & 1) & ((bitplane[i] >> (j+4)) & 1);
            dist = dist + result;
        }
    }

    //Concatenate 4-times
    return dist;
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
                    //Modifying bits for A ('00') Case
                    if((img[i] >> j) == 0 & (img[i] >> (j+1)) == 0){
                        eimage[i] &= ~(1 << j); //setting to 0
                        eimage[i] &= ~(1 << (j+1)); //setting to 0
                    }
                    //Modifying bits for G ('01') Case
                    if((img[i] >> j) == 0 & (img[i] >> (j+1)) == 1){
                        eimage[i] |= (1 << j); //setting to 1
                        eimage[i] &= ~(1 << (j+1)); //setting to 0
                    }
                    //Modifying bits for G ('10') Case
                    if((img[i] >> j) == 1 & (img[i] >> (j+1)) == 0){
                        eimage[i] &= ~(1 << j); //setting to 0
                        eimage[i] |= (1 << (j+1)); //setting to 1
                    } 
                    //Modifying bits for T ('11') Case
                    if((img[i] >> j) == 1 & (img[i] >> (j+1)) == 1){
                        eimage[i] |= (1 << j); //setting to 1
                        eimage[i] |= (1 << (j+1)); //setting to 1
                    }
                    break;
                case 3:
                    //Modifying bits for A ('00') Case
                    if((img[i] >> j) == 0 & (img[i] >> (j+1)) == 0){
                        eimage[i] &= ~(1 << j); //setting to 0
                        eimage[i] |= (1 << (j+1)); //setting to 1
                    }
                    //Modifying bits for G ('01') Case
                    if((img[i] >> j) == 0 & (img[i] >> (j+1)) == 1){
                        eimage[i] &= ~(1 << j); //setting to 0
                        eimage[i] &= ~(1 << (j+1)); //setting to 0
                    }
                    //Modifying bits for G ('10') Case
                    if((img[i] >> j) == 1 & (img[i] >> (j+1)) == 0){
                        eimage[i] |= (1 << j); //setting to 1
                        eimage[i] |= (1 << (j+1)); //setting to 1
                    } 
                    //Modifying bits for T ('11') Case
                    if((img[i] >> j) == 1 & (img[i] >> (j+1)) == 1){
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
    /*----------------------------------------------------VARIABLES-----------------------------------*/
    //Defining the user generated key

    //Check that Key is 168 Bits
    if(!check_size_168(key)){
        return -1;
    }

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

    //Modify the unsigned char variables into 8-bit (uint_8) variables
    vector<uint8_t> bplane(size);
    for (int i = 0; i < width * height * channels; ++i) {
        bplane[i] = image[i];
    }

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
    
    // Calculate the Hamming Distance
    uint64_t(HD) = HammingDistance(bplane);
    cout << "HammingDistance is: " << HD << " or " << intToBinary(static_cast<int>(HD),32);

    /*
    // Check HammingDistance is 168 bit
    if(!check_size_168(HD)){
        return -1;
    }
    */

    //DNA Encode the Original Matrix
    vector<uint8_t> encoded_bp = DNA_Encode_Bitplane(bplane,2);

    
    //Translating the encoded bit-plane back to unsigned char form
    unsigned char* e_img = encoded_image(encoded_bp);

    /*
    // Calculating the key for chaotic cores using XOR
    bitset<168>(skey) = HD^key;

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
