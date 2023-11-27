#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#include "include/ap_int.h"
#include <iostream>
#include <bitset>
#include <vector>
//#include "UserKey_Generation.h"

/*
 REFERENCES:  
 [1] 	Patel, Sakshi & Veeramalai, Thanikaiselvan. (2022). 
 		Image Encryption Using a Spectrally Efficient Halton Logistics Tent (HaLT) Map 
 		and DNA Encoding for Secured Image Communication. Entropy. 
 		24. 803. 10.3390/e24060803. 
 [2] 
*/

using namespace std;

string intToBinary(int value, int numBits) {
    return bitset<32>(value).to_string().substr(32 - numBits);
}

int main() {
	// Read an image from file
    const char* imagePath1 = "test_image.jpeg";
    const char* imagePath2 = "test_image2.jpeg";
    int width1, height1, channels1, width2, height2, channels2;

    unsigned char* image1 = stbi_load(imagePath1, &width1, &height1, &channels1, 0);
    unsigned char* image2 = stbi_load(imagePath2, &width2, &height2, &channels2, 0);
    int size = width1*height1*channels1;
    unsigned char* r_image = new unsigned char[size];
    memset(r_image, 0, size); //initialize with 0s

    int pixel = 320;
    /*
    cout << "Binary Value: " << intToBinary(static_cast<int>(image2[pixel]),8) << endl;
    cout << "Bit 1 is: " << ((image2[pixel] >> 7) & 0x01) << endl;
    cout << "Bit 2 is: " << ((image2[pixel] >> 6) & 0x01) << endl;
    cout << "Bit 3 is: " << ((image2[pixel] >> 5) & 0x01) << endl;
    cout << "Bit 4 is: " << ((image2[pixel] >> 4) & 0x01) << endl;
    cout << "Bit 5 is: " << ((image2[pixel] >> 3) & 0x01) << endl;
    cout << "Bit 6 is: " << ((image2[pixel] >> 2) & 0x01) << endl;
    cout << "Bit 7 is: " << ((image2[pixel] >> 1) & 0x01) << endl;
    cout << "Bit 8 is: " << ((image2[pixel] >> 0) & 0x01) << endl;
	cout << "Binary Value: " << intToBinary(static_cast<int>(image2[pixel]),8) << endl;
	cout << "Masking bit 2 and 3 to 0" << endl;
	image2[pixel] &= ~(1 << 6);
	image2[pixel] &= ~(1 << 5);
    cout << "Bit 2 is: " << ((image2[pixel] >> 6) & 0x01) << endl;
    cout << "Bit 3 is: " << ((image2[pixel] >> 5) & 0x01) << endl;
	cout << "Binary Value: " << intToBinary(static_cast<int>(image2[pixel]),8) << endl;

    cout << static_cast<int>(image1[2] >> 1) << endl;
	*/

    //A = 00 | C = 01 | G = 10 | T = 11
  	for(int i = 0; i <= size; i++){
  		 for (int j = 6; j >= 0; j-=2) {
  		 	if((((image1[i] >> (j+1)) & 0x01) == 0) && (((image1[i] >> j) & 0x01) == 0)){
  		 		//A('00') + A('00')
  		 		if(((image2[i] >> (j+1)) & 0x01) == 0 && ((image2[i] >> j) & 0x01) == 0){
  		 			r_image[i] &= ~(1 << (j+1)); //setting to 0
                    r_image[i] |= (1 << j); //setting to 1
            	}
            	//A('00') + C('01)
  		 		if(((image2[i] >> (j+1)) & 0x01) == 0 && ((image2[i] >> j) & 0x01) == 1){
  		 			r_image[i] &= ~(1 << (j+1)); //setting to 0
                    r_image[i] &= ~(1 << j); //setting to 0
                }
                //A('00') + G('10')
  		 		if(((image2[i] >> (j+1)) & 0x01) == 1 && ((image2[i] >> j) & 0x01) == 0){
  		 			r_image[i] |= (1 << (j+1)); //setting to 1
                    r_image[i] |= (1 << j); //setting to 1
                }
                //A('00') + T('11')
  		 		if(((image2[i] >> (j+1)) & 0x01) == 1 && ((image2[i] >> j) & 0x01) == 1){
  		 			r_image[i] |= (1 << (j+1)); //setting to 1
                    r_image[i] &= ~(1 << j); //setting to 0
                }
  			}
  			if((((image1[i] >> (j+1)) & 0x01) == 0) && (((image1[i] >> j) & 0x01) == 1)){
  				//C('01') + A('00')
  		 		if(((image2[i] >> (j+1)) & 0x01) == 0 && ((image2[i] >> j) & 0x01) == 0){
  		 			r_image[i] &= ~(1 << (j+1)); //setting to 0
                    r_image[i] &= ~(1 << j); //setting to 0
            	}
  				//C('01') + C('01')
  		 		if(((image2[i] >> (j+1)) & 0x01) == 0 && ((image2[i] >> j) & 0x01) == 1){
  		 			r_image[i] &= ~(1 << (j+1)); //setting to 0
                    r_image[i] |= (1 << j); //setting to 1
                }
                //C('01') + G('10')
  		 		if(((image2[i] >> (j+1)) & 0x01) == 1 && ((image2[i] >> j) & 0x01) == 0){
  		 			r_image[i] |= (1 << (j+1)); //setting to 1
                    r_image[i] &= ~(1 << j); //setting to 0
                }
                //C('01') + T('11')
  		 		if(((image2[i] >> (j+1)) & 0x01) == 1 && ((image2[i] >> j) & 0x01) == 1){
  		 			r_image[i] |= (1 << (j+1)); //setting to 1
                    r_image[i] |= (1 << j); //setting to 1
                }
  			}
  			if((((image1[i] >> (j+1)) & 0x01) == 1) && (((image1[i] >> j) & 0x01) == 0)){
  			  	//G('10') + A('00')
  		 		if(((image2[i] >> (j+1)) & 0x01) == 0 && ((image2[i] >> j) & 0x01) == 0){
  		 			r_image[i] |= (1 << (j+1)); //setting to 1
                    r_image[i] |= (1 << j); //setting to 1
            	}
  				//G('10') + C('01')
  		 		if(((image2[i] >> (j+1)) & 0x01) == 0 && ((image2[i] >> j) & 0x01) == 1){
  		 			r_image[i] |= (1 << (j+1)); //setting to 1
                    r_image[i] &= ~(1 << j); //setting to 0
                }
                 //G('10')+ G('10')
  		 		if(((image2[i] >> (j+1)) & 0x01) == 1 && ((image2[i] >> j) & 0x01) == 0){
  		 			r_image[i] &= ~(1 << (j+1)); //setting to 0
                    r_image[i] |= (1 << j); //setting to 1
                }
                //G('10') + T('11')
  		 		if(((image2[i] >> (j+1)) & 0x01) == 1 && ((image2[i] >> j) & 0x01) == 1){
  		 			r_image[i] &= ~(1 << (j+1)); //setting to 0
                    r_image[i] &= ~(1 << j); //setting to 0
                } 				
  			}
  			if((((image1[i] >> (j+1)) & 0x01) == 1) && (((image1[i] >> j) & 0x01) == 1)){
  				//T('11') + A('00')
  		 		if(((image2[i] >> (j+1)) & 0x01) == 0 && ((image2[i] >> j) & 0x01) == 0){
  		 			r_image[i] |= (1 << (j+1)); //setting to 1
                    r_image[i] &= ~(1 << j); //setting to 0
            	}
  				//T('11') + C('01')
  		 		if(((image2[i] >> (j+1)) & 0x01) == 0 && ((image2[i] >> j) & 0x01) == 1){
  		 			r_image[i] |= (1 << (j+1)); //setting to 1
                    r_image[i] |= (1 << j); //setting to 1
                }
                //T('11')+ G('10')
  		 		if(((image2[i] >> (j+1)) & 0x01) == 1 && ((image2[i] >> j) & 0x01) == 0){
  		 			r_image[i] &= ~(1 << (j+1)); //setting to 0
                    r_image[i] &= ~(1 << j); //setting to 0
                }
  				//T('11') + T('11')
  				if(((image2[i] >> (j+1)) & 0x01) == 1 && ((image2[i] >> j) & 0x01) == 1){
  		 			r_image[i] &= ~(1 << (j+1)); //setting to 0
                    r_image[i] |= (1 << j); //setting to 1
                } 	
  			}

  		}
  	}

    cout << "Image 1: Pixel " << pixel << " " << intToBinary(static_cast<int>(image1[pixel]),8) << endl;
    cout << "Image 2: Pixel " << pixel << " " << intToBinary(static_cast<int>(image2[pixel]),8) << endl;
	cout << "Resulting Addition: Pixel " << pixel << " " << intToBinary(static_cast<int>(r_image[pixel]),8) << endl;
   

    stbi_write_jpg("output2.jpg", width1, height1, channels1, r_image, 100);

	return 0;
}