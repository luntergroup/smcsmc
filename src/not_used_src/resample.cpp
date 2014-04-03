#include"resample.hpp"

/*! \fn double unifRand()
 * \ingroup group_systemetic
 * \brief Simulate uniform srandom variable between 0 and 1.
 */
double unifRand(){
    //MTRand_closed return_value;
    //return return_value();
    //return rand()/(double(RAND_MAX)+1);

    return rand() / double(RAND_MAX); //generates a psuedo-random float between 0.0 and 0.999...
} 

double exponSample(double lambda){
    //double random_unit=unifRand();
    //double expsample_rtn=-log( 1- random_unit)/ lambda;
    //std::cout << "random unit " << random_unit<< "  exp sample " <<  expsample_rtn<<std::endl;
    //return expsample_rtn;
    return log( 1-unifRand() )/ lambda;
        
}



/*! \todo to be removed, not used at the moment???*/
void update_cum_sum_array(Vparticle * pointer_to_current_particles, std::valarray<double> & weight_cum_sum,size_t Num_of_particles){
    weight_cum_sum=0; //Reinitialize the cum sum array 
    for (size_t i=0; i<Num_of_particles ;i++){
    //for (size_t i=0; i<pointer_to_current_particles->particles.size() ;i++){
        //update the cum sum array
        weight_cum_sum[i+1]=weight_cum_sum[i]+pointer_to_current_particles->particles[i]->weight();
        //cout<<weight_cum_sum[i]<<endl;
    }
    //check for the cum weight
    for (size_t i=0;i<weight_cum_sum.size();i++){dout<<weight_cum_sum[i]<<"  ";}dout<<std::endl;
}

/*! 
 * @ingroup group_pf_init
 * \brief Calculate the total likelihood for the use of rejection sampling */
double cal_L_hat(Vparticle * init_pointer_to_current_particles, int pool_size){
    double L_hat=0;
    size_t initial_num_particles=init_pointer_to_current_particles->particles.size();
    for (size_t i=0; i< initial_num_particles ;i++){
        L_hat+=init_pointer_to_current_particles->particles[i]->weight();
    }
    return L_hat/initial_num_particles*pool_size;
}








