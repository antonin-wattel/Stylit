#include "stylit.h"


//input: multi - channel images A = { A, A0} and Bk = { B, B0k }


using namespace std;

void run_stylit(std::pair<Image_multichannel, Image> A, std::pair<Image_multichannel, Image> B) {

    int patch_size = 5;
    float mu = 2;//guidance influence
    int downsampling_ratio = 2;
    int num_em_iterations = 6;//optimization iterations for each level
    int pyramid_levels = 2;

    

    //compute pyramid with downsampling ratio 2
    //...

    
    /*std::pair<Image_multichannel, Image_multichannel> A;
    std::pair<Image_multichannel, Image> B;*/

    //coarse to fine
    for (int pyramid_level = 0; pyramid_level < pyramid_levels; pyramid_level++) {

        std::cout << "pyramid depth" << pyramid_level << std::endl;
        //downsample using the pyramid


        //b = rescale_images(b_original, 0.5 * *(num_pyramid_levels - 1 - pyramid_level))
 
        //EM-like iteration scheme to minimize energy
        for (int iteration_step = 0; iteration_step < num_em_iterations; iteration_step++) {
            std::cout << "iteration " << iteration_step << std::endl;

            //initialize the NNF
            std::vector<glm::vec2> NNF;
            NNF.resize(B.second.width() * B.second.height(), glm::vec2(0, 0));

            std::cout << "width" << B.second.width() << std::endl;
#pragma omp parallel for
            //go over patches in B and store the NNF in A
            for (int x = 1; x < B.first.width() - 1; x++) {
                std::cout <<x << std::endl;
#pragma omp parallel for
                for (int y = 1; y < B.first.height() - 1; y++) {
                    glm::vec2 q(x, y);
                    NNF[y * B.second.width() + x] = get_NNF(q, A, B, mu); // stores indices of best matching patches for each pixel in 
                }
            }

#pragma omp parallel for
            //synthesize the image b'
            for (int x = 1; x < B.first.width() - 1; x++) {
#pragma omp parallel for
                for (int y = 1; y < B.first.height() - 1; y++) {
                    glm::vec2 q(x, y);
                    B.second.operator()(x, y) = average(A, NNF, q);
                }
            }
        }
    }

    B.second.save("synthesised_image.png");

}


//USEFUL: TO FIX !
float error(pair<Image_multichannel, Image> A, pair<Image_multichannel, Image> B, glm::vec2 p, glm::vec2 q, float mu){

    //to do: implement feature vector
    /*a_p = A.second.get_patch(p);
    b_p = A.second.get_patch(q);
    ap = A.first.get_patch(p);
    bp = B.first.get_patch(q);
    res = (a_p - b_p).squaredNorm() + mu*(ap-bp).squaredNorm();
    return res;*/

    return 0.f;
    //we probably need an eigen matrix here
}

glm::vec2 get_NNF(glm::vec2 q, pair<Image_multichannel, Image> A, pair<Image_multichannel, Image> B, float mu) {

    //to do: speed this up using patch match

    float min_E = std::numeric_limits<float>::infinity();;
    glm::vec2 argmin(-1, -1);
    for (int x = 1; x < A.first.width(); x++) {
        for (int y = 1; y < A.first.height(); y++) {

            glm::vec2 p(x, y);
            float E = error(A, B, p, q, mu);
            
            if (E < min_E) {
                min_E = E;
                argmin = glm::vec2(x, y);
            }
        }
    }

    return argmin;

}

glm::vec3 average(pair<Image_multichannel, Image> A, std::vector<glm::vec2> NNF, glm::vec2 q) {
    //compute the average color of colocated pixels in neighbour patches
    return glm::vec3(0, 0, 0);
}









//main loop

// this can be of great help !
//https://github.com/Vottivott/3d-style-transfer/blob/master/stylit.py
//
//
//    // Our pyramid uses 2 for the downsampling ratio
//    
//source, target, channel_weights = load_stylit_images()
//
//num_pyramid_levels = 7  # 6
//a_w, a_h = 400, 400
//a = target
//a_guides_original = np.copy(target[:, : , 3 : ])
//b_original = source
//a = rescale_images(a, 0.5 * *(num_pyramid_levels - 1))
//patch_size = 5
//smallest_b_size = (np.array(b_original.shape[:2]) * 0.5 * *(num_pyramid_levels - 1)).astype(int)
//
//b = rescale_images(b_original, smallest_b_size) # b = rescale_images(b_original, 0.5 * *(num_pyramid_levels - 1))
//offsets = init_random_offsets(a, b, patch_size, channel_weights)
//a[:, : , : 3] = reconstruct_image(offsets, a[:, : , : 3], b[:, : , : 3], patch_size//2) # Initial reconstruction from the randomly initialized nearest-neighbor field
//
////we run the synthesis on num_pyramid_levels levels -> what does this mean ?
//for pyramid_level in range(num_pyramid_levels) :
//    
//    b = rescale_images(b_original, 0.5 * *(num_pyramid_levels - 1 - pyramid_level))
//    //6 optimization steps on each level
//    for iteration in range(6) :
//
//        if patch_size >= min(a_w, a_h, b.shape[0], b.shape[1]) ://just a check, not sure we need this
//            continue
//
//        //To accelerate the retrieval of nearest neighbours, we use PatchMatch
//        //what are offsets ????
//        offsets = patchmatch(a, b, offsets, patchmatch_iterations = 4, patch_size = patch_size, iteration_callback = iteration_callback, channel_weights = channel_weights)
//        a[:, : , : 3] = reconstruct_image(offsets, a[:, : , : 3], b[:, : , : 3], patch_size//2)
//        print("Finished texture syntesis iteration " + str(iteration))
//        show_image(a[:, : , : 3])
//
//
//    rescaled_output = rescale_images(a[:, : , : 3], 2.0)
//    save_image(a[:, : , : 3], "results2/output_level_" + str(pyramid_level) + ".png")
//    if pyramid_level < num_pyramid_levels - 1:
//        a = np.concatenate((rescaled_output, rescale_images(a_guides_original, rescaled_output.shape[:2])), 2)
//        b = rescale_images(b_original, smallest_b_size * 2 * *(pyramid_level + 1))
//        offsets = upscaled_offsets(offsets, a, b, patch_size, channel_weights)