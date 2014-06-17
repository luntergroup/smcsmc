/*
 * pf-ARG is short for particle filters for ancestral recombination graphs. 
 * This is a free software for demographic inference from genome data with particle filters. 
 * 
 * Copyright (C) 2013, 2014 Sha (Joe) Zhu and Gerton Lunter
 * 
 * This file is part of pf-ARG.
 * 
 * pf-ARG is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

#include "particleContainer.hpp"

/*! \brief Particle filtering Initialization 
 * Create particle initial states in the simulation
 * 
 * \ingroup group_pf_init
 */ 
ParticleContainer::ParticleContainer(
                    Model* model, 
                    RandomGenerator* rg, 
                    size_t Num_of_states, 
                    //vector <bool> data_for_init_states, 
                    //bool withdata,
                    double initial_position,
                    bool heat_bool ){
    this->heat_bool_ = heat_bool;
    this->random_generator_ = new MersenneTwister(rg->seed());  /*! Initialize random generator for particle filter */
    //this->set_random_generator(rg);
    this->set_ESS(0);
    this->set_current_printing_base(0);    
	dout << " --------------------   Particle Initial States   --------------------" << std::endl;	
	for ( size_t i=0; i < Num_of_states ; i++ ){
		ForestState *  new_state = new ForestState(model,rg);  // create a new state, using scrm; scrm always starts at 0.
//new_state->random_generator_ = new MersenneTwister(rg->seed()+i);  /*! Setting each particle to independent random generator */ //DEBUG
        new_state->setSiteWhereWeightWasUpdated( initial_position );
		new_state->setAncestor ( i );
        this->push(new_state, 1.0/Num_of_states );        
	    }	    
    }


/*! \brief Resampling step
 *  If the effective sample size is less than the ESS threshold, do a resample, currently using systemetic resampling scheme.
 */ 
void ParticleContainer::ESS_resampling(valarray<double> weight_cum_sum, valarray<int> &sample_count, int mutation_at, double ESSthreshold, int num_state){    
    dout << "ESS is " <<  this->ESS() <<", number of particle is " <<  num_state <<endl;
    if (this->ESS() < ESSthreshold){ // resample if the effective sample size is small, to check this step, turn the if statement off
        dout << " ### PROGRESS: ESS resampling" << endl;
        this->systemetic_resampling( weight_cum_sum, sample_count, num_state);
        //this->trivial_resampling ( sample_count, num_state );
        this->resample(sample_count);  
        }    
    }


/*! \brief Particle filtering Resampling step
 * 
 *  Update particles according to the sample counts, and set the particle probabilities to 1
 * 
 * \ingroup group_resample
 */ 
void ParticleContainer::resample(valarray<int> & sample_count){	
    dout << " ****************************** Start making list of new states ****************************** " << std::endl;
	dout << " will make total of " << sample_count.sum()<<" particle states" << endl;
	for (size_t old_state_index = 0; old_state_index < sample_count.size(); old_state_index++){
        if ( (sample_count[old_state_index] > 0)){
            ForestState * current_state = this->particles[old_state_index] ;
            this->push(current_state); // As this pushed step, sets the particle weight to 1, by default value.
            
            // create new copy of the resampled particle 
            for (int ii=2; ii <= sample_count[old_state_index]; ii++){ 
			  	//cout << "       Make a copy of the " << old_state_index << "th particle" << endl;                
				ForestState * new_copy_state = new ForestState( *this->particles[old_state_index] );
				this->push(new_copy_state); // As this pushed step, sets the particle weight to 1, by default value.
                }
            } 
        else {
            delete this->particles[old_state_index]; 
            }
        
        this->particles[old_state_index]=NULL;
        }
    
    this->shifting(sample_count.sum());
    
    dout<<"sample_count.sum() = "<<sample_count.sum()<<endl;
    assert(sample_count.sum() == int(this->particles.size()));
    dout<<this->particles.size()<<endl;
	dout << " ****************************** End of making list of new particles ****************************** " << std::endl;
	assert(this->check_state_orders());
    }


/*!
 * The current implementation append ForestState to the end of the particleContatiner during the resampling step, and the pointers to the ForestState
 * in the previous iteration have been removed. This is house keeping.
 */ 
void ParticleContainer::shifting(int number_of_particles){	
    /*! \todo SHIFTING */
    // SHIFTING AND RESIZING THE PARTICLE ARRAY
    // SHIFTING, this is not pretty, there must be some other ways to do this ...
    for (int i = 0; i < number_of_particles; i++){
        this->particles[i] = this->particles[i+number_of_particles];
        }
    this -> particles.resize(number_of_particles);
    this -> particles.shrink_to_fit();
    // SHIFTING AND RESIZING THE PARTICLE ARRAY FINISHED    
    }


/*! \brief Update particle weight according to the haplotype data
 *	@ingroup group_pf_update
 */
void ParticleContainer::update_state_weights_at_A_single_site(
    double mutation_at,
    double mutation_rate, 
    //bool withdata,
    bool empty_file,
    vector <bool> haplotypes_at_tips
    ){
			
	// now update the weights of all particles, by calculating the likelihood of the data over the previous segment	
	for (size_t particle_i=0; particle_i < this->particles.size(); particle_i++){
		this->particles[particle_i]->include_haplotypes_at_tips(haplotypes_at_tips);

		//double likelihood_of_haplotypes_at_tips = this->particles[particle_i]->calculate_likelihood(withdata); // DEBUG 
        double likelihood_of_haplotypes_at_tips = this->particles[particle_i]->calculate_likelihood( !empty_file ); // DEBUG , if it is not empty_file, calculate the likelihood
        dout << "updated weight =" << this->particles[particle_i]->weight()  << "*" <<  likelihood_of_haplotypes_at_tips <<endl;

        this->particles[particle_i]->setParticleWeight( this->particles[particle_i]->weight() * likelihood_of_haplotypes_at_tips);
		dout << "particle " <<  particle_i<<" done" << endl;
        }
    
    this->normalize_probability(); // It seems to converge slower if it is not normalized ...
	dout << endl;
    }


/*!
 * ParticleContatiner destructor
 */ 
ParticleContainer::~ParticleContainer(){
    delete random_generator_;
    // The following message may show up very often if it is passed by reference ...
    // dout << "ParticleContainer destructor is called" << endl;
    }


/*!
 * Proper particleContatiner destructor, remove pointers of the ForestState
 */ 
void ParticleContainer::clear(){
	// When this is called, this should be the difference between number of forestStates ever built minus ones have already been removed. this should be equal to the size for particles.
    // cout<<"Forest state was created " << new_forest_counter << " times" << endl;  // DEBUG
    // cout<<"Forest state destructor was called " << delete_forest_counter << " times" << endl; // DEBUG
    
    dout << "ParticleContainer clear() is called" << endl;
	for (size_t i = 0; i < this->particles.size(); i++){
		if (this->particles[i]!=NULL){
			delete this->particles[i];
			this->particles[i]=NULL;
            }
        }
	particles.clear();
	dout << "Particles are deleted" << endl;
    }
 

/*!
 * Append new ForestState to the end of the ParticleContainer with weight.
 */ 
void ParticleContainer::push(ForestState* state, double weight){
	state->setParticleWeight(weight);
	this->particles.push_back(state);              
    }


/*! 
 * @ingroup group_pf_resample
 * @ingroup group_pf_update
 * \brief Calculate the effective sample size, and update the cumulative weight of the particles
 */
void ParticleContainer::update_cum_sum_array_find_ESS(std::valarray<double> & weight_cum_sum){
	double wi_sum=0;
	double wi_sq_sum=0;
	double Num_of_states = this->particles.size();
	weight_cum_sum=0; //Reinitialize the cum sum array 
	
	for (size_t i=0; i<Num_of_states ;i++){
		//update the cum sum array
		double w_i=this->particles[i]->weight();
		weight_cum_sum[i+1]=weight_cum_sum[i]+w_i;
		wi_sum = wi_sum + w_i;
		//wi_sq_sum = wi_sq_sum + pow(this->particles[i]->weight,2);
		wi_sq_sum = wi_sq_sum + w_i * w_i;
        }

	//check for the cum weight
	dout << "### particle weights ";
    for (size_t i=0;i<Num_of_states;i++){
        dout << this->particles[i]->weight()<<"  ";
        } dout << std::endl<<std::endl;
	
    dout << "### updated cum sum of particle weight ";
    for (size_t i=0;i<weight_cum_sum.size();i++){
        dout << weight_cum_sum[i]<<"  ";
        } dout << std::endl;

	this->set_ESS(wi_sum * wi_sum / wi_sq_sum);	
    }



/*! 
 * @ingroup group_pf_update
 * \brief Update the current state to the next state, at the given site, update all particles to it's latest genealogy state.  Also include the likelihood for no mutations.
 */
void ParticleContainer::extend_ARGs( double mutation_at, double mutation_rate, bool withdata ){
    dout << endl<<" We are extending particles" << endl<<endl;
 	for (size_t particle_i=0; particle_i < this->particles.size(); particle_i++){
        dout << "We are updating particle " << particle_i << endl;
        /*! 
         * For each particle, extend the current path until the the site such that the next genealogy change is beyond the mutation
         * Invariant: the likelihood is correct up to 'updated_to'
         */
        double updated_to = this->particles[particle_i]->site_where_weight_was_updated();
        dout << "Particle current base is at " << this->particles[particle_i]->current_base() << " weight is updated to " << updated_to <<endl;
        assert (updated_to >= this->particles[particle_i]->current_base());
        while ( updated_to < mutation_at ) {
            
            dout << "  Now at " <<this->particles[particle_i]->current_base()<< " updated_to " << updated_to << " and extending to " << mutation_at << endl;            
            /*!
             * First, update the likelihood up to either mutation_at or the end of this state
             */
            double update_to = min( mutation_at, this->particles[particle_i]->next_base() );
            double length_of_local_tree = this->particles[particle_i]->local_tree_length(); // in generations
            double likelihood_of_segment;
            if (withdata) {
                likelihood_of_segment = exp( -mutation_rate * length_of_local_tree * (update_to - updated_to) ); // assume infinite site model
                //likelihood_of_segment = pow(0.5 - 0.5*exp(-length_of_local_tree*2*mutation_rate ), (update_to - updated_to)); // Finite site model
                dout << " Likelihood of no mutations in segment of length" << (update_to - updated_to) << " is " << likelihood_of_segment << endl;
                } 
            else {
                likelihood_of_segment = 1;
                dout << " no data" << endl;
                }
            this->particles[particle_i]->setParticleWeight( this->particles[particle_i]->weight() * likelihood_of_segment);
            updated_to = update_to;  // rescues the invariant
            
            /*!
             * Next, if we haven't reached mutation_at now, add a new state and iterate
             */
            if (updated_to < mutation_at) {
                //dout<<"calling here"<<endl;
                this->particles[particle_i]->sampleNextGenealogy();
                
                if ( this->heat_bool_ ){
                    TmrcaState tmrca( this->particles[particle_i]->site_where_weight_was_updated(), this->particles[particle_i]->local_root()->height() );
                    this->particles[particle_i]->TmrcaHistory.push_back ( tmrca );
                    }
                
                }
            
            }
        
        assert (updated_to == mutation_at);        
        this->particles[particle_i]->setSiteWhereWeightWasUpdated( mutation_at );
        //cout<<"current_base() = "<<this->particles[particle_i]->current_base()<<" mutation at "<<mutation_at<< " next_base = "<<this->particles[particle_i]->next_base() <<endl;
        }
    
    /*! normalize the probability upon until the mutation */
    //this->normalize_probability(); // This normalization doesn't seem to do much ...
    }


/*!
 * Normalize the particle weight, inorder to prevent underflow problem
 */ 
void ParticleContainer::normalize_probability(){
    double total_probability = 0;
    for ( size_t particle_i = 0;particle_i < this->particles.size(); particle_i++ ){
        total_probability += this->particles[particle_i]->weight();
        }
    for ( size_t particle_i = 0; particle_i < this->particles.size(); particle_i++ ){
        this->particles[particle_i]->setParticleWeight( this->particles[particle_i]->weight() / total_probability);
        }
    }


void ParticleContainer::update_state_to_data(
                        Vcf * VCFfile, 
                        Model * model, 
                        valarray<double> & weight_cum_sum
                        //bool finite_bool
                        ){
    dout <<  " ******************** Update the weight of the particles  ********** " <<endl;
    dout << " ### PROGRESS: update weight at " << VCFfile->site()<<endl;
    double mutation_at = VCFfile->site();
    bool withdata = !VCFfile->missing_data();
    double mutation_rate = model->mutation_rate();
    
    /*!
     * \verbatim
                 
                                  previous VCFfile->site()
                              p0  mut0  p1  p2  p3         p4  mut1
                              .   .     .   .   .          .   . 
                              .   .     .   .   .          .   VCFfile->site()
                              .   .         .   .          .   .
                                  .     1       .              .
                                  .     x---o   .          4   .
                              0   .     |   |   .          x-------o
                              x---------o   |   .          |   .
                              |   .         |              |   .
                              |   .         x---o          |   .
         x--------o           |   .         2   |          |   .
                  |           |   .             x----------o   .
                  |           |   .             3              .
                  x-----------o   .                            .
                                  .                            .
                                  .                            .
      
      \endverbatim
     * At the start of calling this function, the tail ForestState is state 0.
     * Let p_i denotes the position of the start of State i.
     * mut0 was the previous mutation position
     * mut1 is the current mutation position
     * 
     * At the start of this function, likelihood of this particle is updated upto mut0
     * When extending the ARG, we first update the likelihood from mut0 to p1 for 
     * not observing any mutation on this segment of sequence. The new state occurs at p1,
     * until p4. The probability of not observing mutation is updated until mut1.
     * 
     * Call function update_state_weights_at_A_single_site(), update the particle weight for
     * observe mutation at site mut1. The probability of all particles are updated until mut1.
     * 
     * Calculate the cumulated likelihood of all particles and the effective sample size (ESS).
     * 
     * 
     */
    
    //Extend ARGs and update weight for not seeing mutations along the equences
    this->extend_ARGs( mutation_at, mutation_rate, withdata );
    
    //Update weight for seeing mutation at the position 
    this->update_state_weights_at_A_single_site( mutation_at, mutation_rate, VCFfile->empty_file(), VCFfile->vec_of_sample_alt_bool ); 
    ////this->update_state_weights_at_A_single_site( mutation_at, mutation_rate, withdata, VCFfile->vec_of_sample_alt_bool ); // DEBUG
    
    //Update the cumulated probabilities, as well as computing the effective sample size
    this->update_cum_sum_array_find_ESS(weight_cum_sum); 
    }


/*! 
 * @ingroup group_naive 
 * \brief Use simple random sampling to resample
 */
void ParticleContainer::trivial_resampling( std::valarray<int> & sample_count, size_t num_state ){
    sample_count=0;
    for (size_t i=0; i < num_state ;i++){
        size_t index = random_generator()->sampleInt(num_state); 
        sample_count[index]=sample_count[index]+1;
        }
        //cout << sample_count.sum() <<endl;
        assert( sample_count.sum() == num_state );
    }
    

/*! 
 * @ingroup group_systemetic 
 * \brief Use systemetic resampling \cite Doucet2008 to generate sample count for each particle
 */
void ParticleContainer::systemetic_resampling(std::valarray<double> cum_sum, std::valarray<int>& sample_count, int sample_size){
    size_t interval_j = 0;
    size_t sample_i = 0;
    size_t N = sample_size;
    //double u_j = rand() / double(RAND_MAX) / N;
    double u_j = this->random_generator()->sample() / N;
    double cumsum_normalization = cum_sum[cum_sum.size()-1];

    dout << std::endl<<"systematic sampling procedue" << std::endl;
    for (size_t i=0;i<cum_sum.size();i++){dout <<  (cum_sum[i]/cumsum_normalization )<<"  ";}dout << std::endl;    

    sample_count[sample_i] = 0;
    while (sample_i < N) {
        dout << "Is " <<  u_j<<" in the interval of " << std::setw(10)<< (cum_sum[interval_j]/ cumsum_normalization) << " and " << std::setw(10)<< (cum_sum[interval_j+1]/ cumsum_normalization);
        /* invariants: */    
        assert( (cum_sum[interval_j] / cumsum_normalization) < u_j ); 
        assert( sample_i < N );
        /* check whether u_j is in the interval [ cum_sum[interval_j], cum_sum[interval_j+1] ) */
        if ( (sample_i == N) || cum_sum[interval_j+1] / cumsum_normalization > u_j ) {
            sample_count[interval_j] += 1;
            sample_i += 1;
            dout << "  yes, update sample count of particle " << interval_j<<" to " << sample_count[interval_j] <<std::endl;
            u_j += 1.0/double(N);    
            } 
        else {
            dout << "   no, try next interval " << std::endl;
            //assert( sample_i < N-1 );
            interval_j += 1;
            sample_count[ interval_j ] = 0;
            }    
        }
    interval_j=interval_j+1;
    for (;interval_j<N;interval_j++){
        sample_count[ interval_j ] = 0;    
        }
    
    dout << "systematic sampling procedue finished with total sample count " << sample_count.sum()<<std::endl<<std::endl;
    for (size_t i=0;i<sample_count.size();i++){dout << sample_count[i]<<"  ";}  dout << std::endl;
    assert(sample_count.sum()==sample_size);
    }


bool ParticleContainer::appendingStuffToFile( double x_end,  PfParam &pfparam){
    // Record the TMRCA and weight when a heatmap is generated
    if (!pfparam.heat_bool){
        return true;
        }
           /*!
             *  \verbatim
           remove all the state prior to the minimum of
            current_printing_base and 
                previous_backbase     
                .                      backbase                     VCFfile->site()
                .                      .                            .
                .                      .     3                      .
                .                      .     x---o              6   .
                .                  2   .     |   |              x-------o
                .                  x---------o   |              |   .
                .                  |   .         |              |   .
             0  .                  |   .         x---o          |   .
             x---------o           |   .         4   |          |   .
                .      |           |   .             x----------o   .
                .      |           |   .             5              .
                .      x-----------o   .                            .
                .      1               .-------------lag------------.
                .                      .                            .
                Count::update_e_count( .                            ParticleContainer::update_state_to_data(
           x_start = previous_backbase .                              mutation data comes in here
           x_end = backbase            .
                                       .
                                       ParticleContainer::appendingStuffToFile(
                                           x_end = backbase, 
          \endverbatim
         *
         * Likelihood of the particle is updated up until state 6, but because of the lagging we are using
         * report the TMRCA up until state 2
         *           
         */ 
    if (x_end < this->current_printing_base()){
        return true;    
        }
    do {
        //this->set_current_printing_base(x_end);

        if (this->current_printing_base() > 0){

            ofstream TmrcaOfstream   ( pfparam.TMRCA_NAME.c_str()    , ios::out | ios::app | ios::binary);
            ofstream WeightOfstream  ( pfparam.WEIGHT_NAME.c_str()   , ios::out | ios::app | ios::binary); ;
            //ofstream BLOfstream      ( pfparam.BL_NAME.c_str()       , ios::out | ios::app | ios::binary);;
            ofstream SURVIVORstream  ( pfparam.SURVIVOR_NAME.c_str() , ios::out | ios::app | ios::binary);
            
            TmrcaOfstream  << this->current_printing_base();
            WeightOfstream << this->current_printing_base();
            //BLOfstream     << this->current_printing_base();
            SURVIVORstream << this->current_printing_base();
            
            for ( size_t i = 0; i < this->particles.size(); i++){
                ForestState * current_state_ptr = this->particles[i];
                WeightOfstream <<"\t" << current_state_ptr->weight();
                SURVIVORstream <<"\t" << current_state_ptr->ancestor();
                                
                //TmrcaOfstream  << "\t" << current_state_ptr->local_root()->height() / (4 * current_state_ptr->model().default_pop_size); // Normalize by 4N0
                double current_tmrca = current_state_ptr->local_root()->height();
                for (size_t tmrca_i = current_state_ptr->TmrcaHistory.size(); tmrca_i > 0; tmrca_i--){
                    if (current_state_ptr->TmrcaHistory[tmrca_i].base < this->current_printing_base()){
                        break;
                        }
                    current_tmrca = current_state_ptr->TmrcaHistory[tmrca_i].tmrca ;
                    }
                TmrcaOfstream  << "\t" << current_tmrca / (4 * current_state_ptr->model().default_pop_size); // Normalize by 4N0

                //BLOfstream     << "\t" << current_state_ptr->local_tree_length()    / (4 * current_state_ptr->model().default_pop_size); // Normalize by 4N0
                current_state_ptr=NULL;
                }
                
            TmrcaOfstream  << endl;
            WeightOfstream << endl;
            //BLOfstream     << endl;
            SURVIVORstream << endl;
            
            TmrcaOfstream.close();
            WeightOfstream.close();
            //BLOfstream.close();
            SURVIVORstream.close();
            }
        this->set_current_printing_base(this->current_printing_base() + pfparam.heat_seq_window);
        } while ( this->current_printing_base() < x_end);
    return true;
    }
    
    
void ParticleContainer::set_particles_with_random_weight(){
    for (size_t i = 0; i < this->particles.size(); i++){
        this->particles[i]->setParticleWeight( this->random_generator()->sample() );
        }
    }


void ParticleContainer::cumulate_recomb_opportunity_at_seq_end( double seqend ){
    for (size_t i = 0; i < this->particles.size(); i++){
        double opportunity_x = seqend - this->particles[i]->current_base();
        double opportunity_y = this->particles[i]->local_tree_length();
        double recomb_opportunity = opportunity_x * opportunity_y;
        this->particles[i]->record_Recombevent(0, 
                                               //0, 0, 
                                               recomb_opportunity, NOEVENT);        
        }
    }


void ParticleContainer::print_particle_probabilities(){
    for (size_t i = 0; i < this->particles.size(); i++){
        cout<<"weight = "<<this->particles[i]->weight()<<endl;
        }
    }
