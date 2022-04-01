#include "stylit.h"
using namespace std;


void run_stylit(std::pair<Image_multichannel, Image_multichannel> &A, std::pair<Image_multichannel, Image_multichannel> &B) {

    std::chrono::high_resolution_clock clock;

    float mu = 2;
    int downsampling_ratio = 2;
    int num_em_iterations = 4;
    size_t pyramid_levels = 6;

    ////compute image pyramids
    Pyramid ap(pyramid_levels, A.first);
    Pyramid a_p(pyramid_levels, A.second);
    Pyramid bp(pyramid_levels, B.first);
    Pyramid b_p(pyramid_levels, B.second);
    std::cout << "created pyramids" << std::endl;
   
    //save images pyramid for report
    // for (int level=0; level<pyramid_levels; level ++){
    //     ap[level].save("pyramid"+std::to_string(level)+".png", 0);
    // }

    //Coarse to fine resolution
    for (int pyramid_level = pyramid_levels-1; pyramid_level >= 0; pyramid_level--) {

        std::cout << "------------------------------- pyramid level = " << pyramid_level << std::endl;
        std::chrono::time_point<std::chrono::high_resolution_clock> before_level = clock.now();//clock

        std::pair<Image_multichannel, Image_multichannel> Ap(ap[pyramid_level], a_p[pyramid_level]);
        std::pair<Image_multichannel, Image_multichannel> Bp(bp[pyramid_level], b_p[pyramid_level]);

        //EM-like iteration scheme to minimize energy
        for (int iteration_step = 0; iteration_step < num_em_iterations; iteration_step++) {
            std::cout << "iteration " << iteration_step << std::endl;
            std::chrono::time_point<std::chrono::high_resolution_clock> before_iteration = clock.now();//clock

            //Initialize the NNF
            std::vector<glm::vec2> NNF;
            NNF.resize(Bp.second.width() * Bp.second.height(), glm::vec2(0, 0));

            //----------------------------------------------------------------------
            //Go over patches in B and store the NNF in A
#pragma omp parallel for collapse(2)
            for (int x = 1; x < Bp.first.width()-1; x++) {
                for (int y = 1; y < Bp.first.height()-1; y++) {
                    glm::vec2 q(x, y);
                    NNF[y * Bp.second.width() + x] = get_NNF(q, Ap, Bp, mu); 
                }
            }
            std::chrono::time_point<std::chrono::high_resolution_clock> after_iteration = clock.now();
            double elapsed = (double)std::chrono::duration_cast<std::chrono::milliseconds>(after_iteration - before_iteration).count();
            std::cout<<"iteration "+std::to_string(iteration_step)+" :" + std::to_string(elapsed) + "ms"<<std::endl; 
            //----------------------------------------------------------------------


            //----------------------------------------------------------------------
            //synthesize the image b'
            //std::cout <<"synthesizing"<< std::endl;
#pragma omp parallel for collapse(2)
            for (int x = 0; x < Bp.first.width(); x++) {
                for (int y = 0; y < Bp.first.height(); y++) {
                    glm::vec2 q(x, y);
                    Bp.second.operator()(0, x, y) = average(Ap, NNF, q);
                }
            }
            //----------------------------------------------------------------------
           
            //upsample b' and put it in the next pyramid 
            if (pyramid_level >0)
                b_p[pyramid_level - 1] = Bp.second.upsample(2);

            Bp.second.save("progress/synthesised_"+std::to_string(pyramid_level)+"_" + std::to_string(iteration_step) + ".png", 0);
        }  

        std::chrono::time_point<std::chrono::high_resolution_clock> after_level = clock.now();
        double elapsedTime = (double)std::chrono::duration_cast<std::chrono::milliseconds>(after_level - before_level).count();
        std::cout<<"level "+std::to_string(pyramid_level)+" :" + std::to_string(elapsedTime) + "ms"<<std::endl; 
    }
}


float error(const pair<Image_multichannel, Image_multichannel> & A, const pair<Image_multichannel, Image_multichannel>& B, const glm::vec2 &p, const glm::vec2 &q, const float & mu){
   
    MatrixXd a_p ;
    MatrixXd b_p;
    MatrixXd ap ;
    MatrixXd bp;
    A.second.get_patch(p, a_p);
    B.second.get_patch(q, b_p);
    A.first.get_patch(p, ap);
    B.first.get_patch(q, bp);

    return (a_p - b_p).squaredNorm() + mu*(ap-bp).squaredNorm();
}


glm::vec2 get_NNF(const glm::vec2 & q, const pair<Image_multichannel, Image_multichannel>&A, const pair<Image_multichannel, Image_multichannel>&B, const float& mu) {
    //brute force searching -> way too slow
    //TO DO: speed this up using patch match
    float min_E = std::numeric_limits<float>::infinity();;
    glm::vec2 min(0, 0);
    for (int x = 1; x < A.first.width()-1; x++) {
        for (int y = 1; y < A.first.height()-1; y++) {

            glm::vec2 p(x, y);
            float E = error(A, B, p, q, mu);
            
            if (E < min_E) {
                min_E = E;
                min = p;
            }
        }
    }
    return min;
}


glm::vec3 average(const pair<Image_multichannel, Image_multichannel>& A, const std::vector<glm::vec2>& NNF, const glm::vec2& q) {
    //compute the average color of colocated pixels in neighbour patches
    glm::vec3 avg(0., 0., 0.);
    MatrixXd a_p;
    glm::vec2 q_ = q;

    //deal with borders
    if (q.x == 0) q_.x = 1;
    if (q.x == A.first.width()-1) q_.x -= 1;
    if (q.y == 0) q_.y = 1;
    if (q.y == A.first.height()-1) q_.y -= 1;

    //get average
    A.second.get_patch(NNF[q_.y * A.second.width() + q_.x], a_p);
    for (int k = 0; k < 3; k++)
        avg[k] = (a_p.block<3, 3>(0, k * 3)).mean();
    return avg;
   
    //glm::vec3 avg(0., 0., 0.);
    //int num_neigh = 0;
    //MatrixXd a_p;
    ////get neighborhood patches
    //for (int i = -1; i <= 1; i++) {
    //    for (int j = -1; j <= 1; j++) {
    //        int px = q.x + i;
    //        int py = q.y + j;
    //        if ((2 <= px < (A.second.width() - 2)) && (2 <= py < (A.second.height() - 2))) {
    //            //std::cout << "ij"<< i<<", "<< std::endl;     
    //            //A.second.get_patch(NNF[q.y * A.second.width() + q.x], a_p); // make sure about this !!!
    //            A.second.get_patch(NNF[py * A.second.width() + px], a_p); // make sure about this !!!
    //            for (int k = 0; k < 3; k++)
    //                avg[k] += (a_p.block<3, 3>(0, k * 3)).mean();
    //            num_neigh++;
    //        }
    //    }
    //}
    //return avg/float(num_neigh);
}



