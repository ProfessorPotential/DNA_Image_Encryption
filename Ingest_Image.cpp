#include <iostream>
#include <opencv2/opencv.hpp>

using namespace std;
using namespace cv;

int main() {
    // Read an image from file
    Mat inputImage = imread("path/to/your/image.jpg");

    // Check if the image was successfully loaded
    if (inputImage.empty()) {
        cout << "Error: Could not read the image." << endl;
        return -1;
    }

    // Display the original image
    imshow("Original Image", inputImage);
    waitKey(0);

    // Convert the image to a matrix (if not already in matrix format)
    Mat imageMatrix = inputImage;

    // Now 'imageMatrix' is your image represented as a matrix
    // You can access pixel values and perform operations on the matrix

    return 0;
}
