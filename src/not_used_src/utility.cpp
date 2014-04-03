#include<utility.hpp>
//resampling




/*! \fn double unifRand()
 * \brief Simulate random variable between 0 and 1.
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
