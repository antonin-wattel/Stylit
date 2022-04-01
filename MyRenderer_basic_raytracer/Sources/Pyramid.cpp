#include "Pyramid.h"



//const Image_multichannel Pyramid::upsample(const Image_multichannel& A) { //these should go in the Image_multichannel class
//
////https://www.geeksforgeeks.org/spatial-resolution-down-sampling-and-up-sampling-in-image-processing/
//
//}


const Image_multichannel Pyramid::gaussian_filter(const Image_multichannel& A) { //these should go in the Image_multichannel class

//https://stackoverflow.com/questions/42186498/gaussian-blur-image-processing-c
//(best and easy would be to use opencv)

    Image_multichannel res(A.width(), A.height(), A.num_channels());
    return A;


    //TO DO: USE A BIGGER KERNEL
    /*int kernel[3][3] = { 1, 2, 1,
                         2, 4, 2,
                         1, 2, 1 };

    for (int channel = 0; channel < A.num_channels(); channel++) {
        for (int x = 0; x < A.width(); x++) {
            for (int y = 0; y < A.height(); y++) {
                for (int k = 0; k < 3; k++) {
                    
                    res.operator()(channel, x, A.height() - 1 - y)[k] = accessPixel(arr, col, row, k, width, height);

                }
            }
        }

    }*/
}


 