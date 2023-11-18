#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#include <iostream>
#include <bitset>
#include <vector>

using namespace std;

/*
This will only work with a static size image. The input function is agnostic to sizing, but does not resize
*/

// Function to convert an integer to binary representation
string intToBinary(int value, int numBits) {
    return bitset<32>(value).to_string().substr(32 - numBits);
}

// Function to calculate the hamming distance. Takes in a matrixed bitplane (8 bits) and outputs a 168 bit hamming distance
string HammingDistance(vector<vector<string> > bitplane) {
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
    return bitset<168>(dist).to_string();
}

int main() {
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

    cout << "1D Size: " << sizeof(image) << endl;

    
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
    vector<vector<string> > bitplanes(8, vector<string>(width * height));

    for (int i = 0; i < width * height; ++i) {
        for (int j = 0; j < 8; ++j) {
            // Extracting the jth (8-bit) from the image
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
    string HD = HammingDistance(bitplanes);
    cout << "168 Bit Hamming Distance:" << HD << endl;

    // Free the allocated memory for the image
    stbi_image_free(image);

    return 0;
}
