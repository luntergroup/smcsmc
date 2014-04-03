#include"particle.hpp"

double unifRand();

double exponSample(double lambda);



//////////////////////////////////////////////////////////////////////////////////////////////////////
//                                NOT USED, FOR REMOVAL
//////////////////////////////////////////////////////////////////////////////////////////////////////


double cal_L_hat(Vparticle * init_pointer_to_current_particles, int pool_size);

void update_cum_sum_array(Vparticle * pointer_to_current_particles,std::valarray<double> & weight_cum_sum,size_t Num_of_particles); // REMOVE???
