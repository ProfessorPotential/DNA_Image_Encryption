#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#include <iostream>

using namespace std;

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

    // Access pixel values and perform operations
    // For example, print the RGB values of the pixel at (x=10, y=10)
    int x = 25;
    int y = 25;
    int index = (y * width + x) * channels;
    cout << "Pixel at (" << x << ", " << y << "): ";
    for (int c = 0; c < channels; ++c) {
        cout << static_cast<int>(image[index + c]) << " ";
    }
    cout << endl;

    // Free the allocated memory for the image
    stbi_image_free(image);

    return 0;
}
