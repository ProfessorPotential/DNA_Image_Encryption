#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
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

/*-------------------------------FUNCTIONS------------------------------*/

// Function to convert an integer to binary representation
string intToBinary(int value, int numBits) {
    return bitset<32>(value).to_string().substr(32 - numBits);
}

// Function to calculate the hamming distance. Takes in a matrixed bitplane (8 bits) and outputs a 168 bit hamming distance
bitset<168> HammingDistance(vector<vector<string> > bitplane) {
    //Collector for Hamming Distance
    int dist = 0;

    //Length of the Bitplanes
    int length = bitplane[0].size();

    for (int i = 0; i < length; i++){
        for (int j = 0; j < 4; j++){
            int result = int(bitplane[j][i]==bitplane[j+4][i]);
            dist = dist + result;
        }
    }
    bitset<42>(short_key) = bitset<42>(dist);
    bitset<168>(long_key);

    //Concatenate 4-times
    for (int i = 0; i < 4; ++i) {
        long_key <<= 42;  // Left shift by 42 bits
        long_key |= bitset<168>(short_key.to_string()); //OR operator
    }
    return bitset<168>(long_key);
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
vector<vector<string> > DNA_Encode_Bitplane(vector<vector<string> > img, int r){
    //A = 00 | G = 01 | C = 10 | T = 11
    vector<vector<string> > eimage = img;
    for (int i = 0; i < img[0].size(); ++i) {
        for (int j = 0; j < img.size(); j+=2) {    
            switch (r) {
                case 1:
                    eimage[j][i] = img[j][i];
                    eimage[j+1][i] = img[j+1][i];
                    break;
                case 2:
                    if(eimage[j][i] == "0" & eimage[j+1][i] == "0"){
                        eimage[j][i] = "0";
                        eimage[j+1][i] = "0";
                    }
                    if(eimage[j][i] == "0" & eimage[j+1][i] == "1"){
                        eimage[j][i] = "1";
                        eimage[j+1][i] = "0";
                    }
                    if(eimage[j][i] == "1" & eimage[j+1][i] == "0"){
                        eimage[j][i] = "0";
                        eimage[j+1][i] = "1";
                    } 
                    if(eimage[j][i] == "1" & eimage[j+1][i] == "1"){
                        eimage[j][i] = "1";
                        eimage[j+1][i] = "1";
                    }
                    break;
                case 3:
                    if(eimage[j][i] == "0" & eimage[j+1][i] == "0"){
                        eimage[j][i] = "0";
                        eimage[j+1][i] = "1";
                    }
                    if(eimage[j][i] == "0" & eimage[j+1][i] == "1"){
                        eimage[j][i] = "0";
                        eimage[j+1][i] = "0";
                    }
                    if(eimage[j][i] == "1" & eimage[j+1][i] == "0"){
                        eimage[j][i] = "1";
                        eimage[j+1][i] = "1";
                    } 
                    if(eimage[j][i] == "1" & eimage[j+1][i] == "1"){
                        eimage[j][i] = "1";
                        eimage[j+1][i] = "0";
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

unsigned char* encoded_image(vector<vector<string> > img){
    //Initializing Variables
    vector<unsigned char> result(img[0].size(),0);

    //Building 8-Bit values
    for (int i = 0; i < img[0].size(); ++i) {
        std::bitset<8> bits; 
        for (int j = 0; j < img.size(); ++j){
            bits[j] = stoi(img[j][i]);
        }
    result[i] = static_cast<unsigned char>(bits.to_ulong()); //converting to integer
    }

    //Translating vector into unsigned char
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
    
    // Access pixel values and perform operations
    // For example, print the RGB values of the pixel at (x,y)
    int x = 25;
    int y = 25;

    int index = (y * width + x) * channels;
    cout << "Pixel at (" << x << ", " << y << "): " << endl;
    for (int c = 0; c < channels; ++c) {
        cout << "Integer: " << static_cast<int>(image[index + c]) << " ";
        cout << "Binary: " << intToBinary(static_cast<int>(image[index + c]),8) << " ";
        cout << endl;
    }
    cout << endl;
    
    //Converting the image to a 2D binary bit plane. 
    //Dimension 0 = Bit place (1 through 8)
    //Dimension 1 = 1D image pixel index 
    vector<vector<string> > bitplanes(8, vector<string>(width * height * channels));
    for (int i = 0; i < width * height * channels; ++i) {
        for (int j = 0; j < 8; ++j) {
            // Extracting the jth from the image
            int bit = (image[i] >> j) & 1;
            // Store the bit in the corresponding bitplane
            bitplanes[j][i] = to_string(bit);
        }
    }
    /*
    // Print the size of each bitplane
    for (int i = 0; i < 8; ++i) {
        cout << "Size of Bitplane " << i << ": " << bitplanes[i].size() << endl;
    }
    */

    // Calculate the Hamming Distance
    bitset<168>(HD) = HammingDistance(bitplanes);
    // Check HammingDistance is 168 bit
    if(!check_size_168(HD)){
        return -1;
    }

    //DNA Encode the Original Matrix
    vector<vector<string> > encoded_bp = DNA_Encode_Bitplane(bitplanes,4);

    //Translating the encoded bit-plane back to unsigned char form
    unsigned char* e_img = encoded_image(encoded_bp);

    // Calculating the key for chaotic cores using XOR
    bitset<168>(skey) = HD^key;

    //cout << "168 Bit Hamming Distance:" << endl << HD << endl;
    //cout << "User Key: " << endl << key << endl;
    //cout << "Special Key:" << endl << skey << endl;

    // Write Image for Debugging
    stbi_write_jpg("output.jpg", width, height, channels, e_img, 100);

    // Free the allocated memory for the image
    stbi_image_free(image);

    return 0;
}
