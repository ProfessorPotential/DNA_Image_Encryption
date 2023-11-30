#define STB_IMAGE_IMPLEMENTATION
#include "./src/stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "./src/stb_image_write.h"

#include "./src/xcoupled_skew_tent.h"
#include "./src/DNA.h"

using namespace std;
/*
 * References:
 * [1] R. Elmanfaloty, A. Alnajim, and E. Abou-Bakr, “A finite precision
 *     implementation of an image encryption scheme based on dna encoding
 *     and binarized chaotic cores.” IEEE Access, Access, IEEE, vol. 9, pp.
 *     136 905 – 136 916, 2021.
 * [2] R. A. Elmanfaloty and E. Abou-Bakr, “Random property enhancement of
 *     a 1D chaotic PRNG with finite precision implementation,” Chaos Solitons
 *     Fractals, vol. 118, pp. 134–144, 2019.
 */

#define SECRET_KEY_LEN 24 // 26 * 8 = 192
#define GEN_KEY_LEN 4
#define SECRET_KEY_GROUP_LEN (SECRET_KEY_LEN / GEN_KEY_LEN)
#define KEY_CHAR_SIZE 8

uint64_t HammingDistance2(uint8_t* image, uint32_t height, uint32_t width){
    uint8_t xor_bits;
    uint64_t h_dist = 0ULL;

    for(uint32_t y=0; y<height; y++){
        for(uint32_t x=0, dx=1; x < width && dx < width; x +=2, dx += 2){
            xor_bits = image[x  + (y * width)] ^ image[dx + (y * width)];
            while(xor_bits){
                h_dist += xor_bits % 2;
                xor_bits >>= 1; 
            }
        }
    }
    return h_dist;
}

int main() {
    // Read an image from file
    const char* imagePath = "./input_image/test_image.jpeg";
    int width, height, channels;

    unsigned char* image = stbi_load(imagePath, &width, &height, &channels, 1);
    unsigned char output_image[width*height];

    /* Pseudo random sequence generation */
    cout <<" *** Generating random sequences *** " << endl;
    uint64_t dist = HammingDistance2((uint8_t*) image, height, width);
    
    uint8_t secret_key1[SECRET_KEY_LEN+1] = "tH15_iS_S3cReT_k3Y_nO.1*";
    uint8_t secret_key2[SECRET_KEY_LEN+1] = "This-1Z-sECr3t-KEV-N0_2!";
    
    uint64_t K1[GEN_KEY_LEN], K2[GEN_KEY_LEN]; // K[4] = {x, y, p1, p2};
    for(int i=0; i<GEN_KEY_LEN; i++){
        K1[i] = dist;
        K2[i] = dist;
        uint64_t sk1=0;
        uint64_t sk2=0;
        for(int j=0; j<SECRET_KEY_GROUP_LEN; j++){
            sk1 |= (uint64_t) secret_key1[i * SECRET_KEY_GROUP_LEN + j] << (KEY_CHAR_SIZE * j);
            sk2 |= (uint64_t) secret_key2[i * SECRET_KEY_GROUP_LEN + j] << (KEY_CHAR_SIZE * j);
        }
        K1[i] ^= sk1;
        K2[i] ^= sk2;
    }
    uint32_t permute[width * height];
    uint8_t random_sequence[width * height], random_sequence_enc[width * height], image_enc[width*height], encrypted_img[width*height];// output_img[width*height];

    generate_sequence_xcoupled_skew_tent(K1[0], K1[1], K1[2], K1[3], width*height, random_sequence);
    generate_permutation_xcoupled_skew_tent(K2[0], K2[1], K2[2], K2[3], width*height, permute);
    
    int img_dim = width*height;
    /* Testing Permutation (test only not part of the algorithm)*/
    cout <<" *** Testing random sequence generation by performing pixel scrambling *** " << endl;
    for(int i=0; i<img_dim; i++){
        unsigned char tmp = image[i];
        image[i] = image[permute[i]];
        image[permute[i]] = tmp;
    }
    stbi_write_jpg("./output_image/permuted.jpg", width, height, 1, image, 100);

    for(int i=img_dim - 1; i>=0; i--){
        unsigned char tmp = image[i];
        image[i] = image[permute[i]];
        image[permute[i]] = tmp;
    }
    stbi_write_jpg("./output_image/unpermuted.jpg", width, height, 1, image, 100);


    /* Encryption */
    cout <<" *** Encrypting image *** " << endl;
    DNA_Encryption(random_sequence, random_sequence_enc, img_dim, 5);
    DNA_Encryption(image, image_enc, img_dim, 5);
    DNA_add(image_enc, random_sequence_enc, encrypted_img, img_dim);

    for(int i=0; i<width * height; i++){
        unsigned char tmp = encrypted_img[i];
        encrypted_img[i] = encrypted_img[permute[i]];
        encrypted_img[permute[i]] = tmp;
    }

    stbi_write_jpg("./output_image/encrypted_img.jpg", width, height, 1, encrypted_img, 100);
    
    /* Decryption */
    cout <<" *** Decrypting image *** " << endl;
    uint8_t decrpyted_img_enc[img_dim], decrpyted_img_dec[img_dim], decrpyted_img[img_dim];
    for(int i=0; i<height*width; i++){
        decrpyted_img_enc[i] = encrypted_img[i];
    }
    for(int i=img_dim - 1; i>=0 ;i--){
        unsigned char tmp = decrpyted_img_enc[i];
        decrpyted_img_enc[i] = decrpyted_img_enc[permute[i]];
        decrpyted_img_enc[permute[i]] = tmp;
    }

    DNA_subtract(decrpyted_img_enc, random_sequence_enc, decrpyted_img_dec, img_dim);
    DNA_Encryption(decrpyted_img_dec, decrpyted_img, img_dim, 5);
    stbi_write_jpg("./output_image/decrypted_img.jpg", width, height, 1, decrpyted_img, 100);
    return 0;
}
