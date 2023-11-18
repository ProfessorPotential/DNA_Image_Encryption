#include <iostream>
#include <opencv2/opencv.hpp>

using namespace std;
using namespace cv;

void encryptRGBImage(double x, int m, int n) {
    Mat r(m, n, CV_8UC1);
    Mat g(m, n, CV_8UC1);
    Mat b(m, n, CV_8UC1);

    r.at<uchar>(0) = static_cast<uchar>(x + 0.000000001);
    g.at<uchar>(0) = static_cast<uchar>(x);
    b.at<uchar>(0) = static_cast<uchar>(x - 0.000000001);

    for (int i = 1; i < m * n; ++i) {
        r.at<uchar>(i) = static_cast<uchar>(1 - 2 * r.at<uchar>(i - 1) * r.at<uchar>(i - 1));
        g.at<uchar>(i) = static_cast<uchar>(1 - 2 * g.at<uchar>(i - 1) * g.at<uchar>(i - 1));
        b.at<uchar>(i) = static_cast<uchar>(1 - 2 * b.at<uchar>(i - 1) * b.at<uchar>(i - 1));
    }

    for (int i = 0; i < m * n; ++i) {
        r.at<uchar>(i) = static_cast<uchar>(round(r.at<uchar>(i) * 98273917129) % 256);
        g.at<uchar>(i) = static_cast<uchar>(round(g.at<uchar>(i) * 98273917129) % 256);
        b.at<uchar>(i) = static_cast<uchar>(round(b.at<uchar>(i) * 98273917129) % 256);
    }

    Mat encryptedImage(m, n, CV_8UC3);
    for (int i = 0; i < m * n; ++i) {
        encryptedImage.at<Vec3b>(i)[0] = r.at<uchar>(i);
        encryptedImage.at<Vec3b>(i)[1] = g.at<uchar>(i);
        encryptedImage.at<Vec3b>(i)[2] = b.at<uchar>(i);
    }

    imshow("Encrypted Image", encryptedImage);
    waitKey(0);
}

int main() {
    double x = 0.5; // You can change the initial value as needed
    int m = 256;   // You can change the dimensions as needed
    int n = 256;

    encryptRGBImage(x, m, n);

    return 0;
}
