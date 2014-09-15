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
                    //size_t random_seed,
                    MersenneTwister *rg,
                    size_t Num_of_states, 
                    //vector <bool> data_for_init_states, 
                    //bool withdata,
                    double initial_position,
                    bool heat_bool ){
    this->heat_bool_ = heat_bool;
    //this->random_generator_ = new MersenneTwister( random_seed );  /*! Initialize random generator for particle filter */
    this->random_generator_ = rg;
    this->set_ESS(0);
    this->set_current_printing_base(0);    
	dout << " --------------------   Particle Initial States   --------------------" << std::endl;	
	for ( size_t i=0; i < Num_of_states ; i++ ){
        size_t new_seed = (size_t) this->random_generator_->sampleInt( INT_MAX );
		RandomGenerator* new_rg = new MersenneTwister( new_seed , this->random_generator_->ff() ); 
		//RandomGenerator* new_rg = new MersenneTwister( new_seed ); 

        Model * new_model =  new Model(*model);
		ForestState* new_state = new ForestState( new_model, new_rg );  // create a new state, using scrm; scrm always starts at 0.  Use a random generator, and model per particle for multithreading
        //ForestState* new_state = new ForestState( model, new_rg );  // create a new state, using scrm; scrm always starts at 0.  Use a random generator, and model per particle for multithreading
        
        // Initialize members of FroestState (derived class members)
        new_state->init_EventContainers( model );    
        new_state->buildInitialTree();
        new_state->setSiteWhereWeightWasUpdated( initial_position );
		new_state->setAncestor ( i );
        if ( this->heat_bool_ ){
            TmrcaState tmrca( 0, new_state->local_root()->height() );
            new_state->TmrcaHistory.push_back ( tmrca );
            }
        this->push(new_state, 1.0/Num_of_states );        
	    }	    
    }


/*! \brief Resampling step
 *  If the effective sample size is less than the ESS threshold, do a resample, currently using systemetic resampling scheme.
 */ 
void ParticleContainer::ESS_resampling(valarray<double> weight_cum_sum, valarray<int> &sample_count, int mutation_at, double ESSthreshold, int num_state){    
    //dout << "ESS is " <<  this->ESS() <<", number of particle is " <<  num_state <<endl;
    dout << "ESS is " <<  this->ESS() <<", number of particle is " <<  num_state << ", and ESSthreshold is " << ESSthreshold <<endl;
    double ESS_diff = ESSthreshold - this->ESS();
    if ( ESS_diff > 1e-6 ){ // resample if the effective sample size is small, to check this step, turn the if statement off
        dout<<"ESS_diff = " << ESS_diff<<endl;
        dout << " ### PROGRESS: ESS resampling" << endl;
        this->systemetic_resampling( weight_cum_sum, sample_count, num_state);
        //this->trivial_resampling ( sample_count, num_state );
        this->resample(sample_count);  
        }    
    }


void ParticleContainer::resample_for_check(valarray<int> & sample_count){	
    this->duplicate_particles ( sample_count );
    size_t Num_of_states = this->particles.size();
    
    for ( size_t i = 0 ;  i < Num_of_states; i++ ){
        if ( this->particles[i]->ForestState_copies.size() == 0 ){
            delete this->particles[i];
            this->particles[i]=NULL;
            continue;
            }
        for ( size_t j = 0; j < this->particles[i]->ForestState_copies.size(); j++){
            this->push( this->particles[i]->ForestState_copies[j] );
            }
        this->particles[i] = NULL;
        }
            
    this->shifting(sample_count.sum());
    assert(sample_count.sum() == int(this->particles.size()));
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
    dout<<"resampling is called"<<endl;
	for (size_t old_state_index = 0; old_state_index < sample_count.size(); old_state_index++){
        if ( sample_count[old_state_index] > 0 ) {
            ForestState * current_state = this->particles[old_state_index];
			// we need at least one copy of this particle; it keeps its own random generator
            this->push(current_state); // The 'push' implementation sets the particle weight to 1
            // create new copy of the resampled particle 
            for (int ii = 2; ii <= sample_count[old_state_index]; ii++) { 
			  	dout << "       Make a copy of the " << old_state_index << "th particle" << endl;                
				ForestState* new_copy_state = new ForestState( *this->particles[old_state_index] );
				// Give the new particle its own random generator (for multithreading)
				size_t new_seed = (size_t) this->particles[old_state_index]->random_generator_->sampleInt( INT_MAX );
                new_copy_state->random_generator_ = new MersenneTwister( new_seed , this->random_generator_->ff() ); 
                //new_copy_state->random_generator_ = new MersenneTwister( new_seed ); 
                new_copy_state->model_ = new Model(*this->particles[old_state_index]->model_);
                //cout << "new random seed is " << this->particles[old_state_index]->random_generator_->seed() + this->particles.size() << endl;
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



/*!
 * ParticleContatiner destructor
 */ 
ParticleContainer::~ParticleContainer(){
    //delete random_generator_;
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
            //delete this->particles[i]->random_generator_; //MULTITRHREADING
            //this->particles[i]->random_generator_ = NULL;
			delete this->particles[i];
			this->particles[i]=NULL;
            }
        }
	this->particles.clear();
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

	for (size_t i=0; i < Num_of_states ;i++){
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
                        VariantReader * VCFfile, 
                        Model * model, 
                        valarray<double> & weight_cum_sum
                        //bool finite_bool
                        ){
    dout <<  " ******************** Update the weight of the particles  ********** " <<endl;
    dout << " ### PROGRESS: update weight at " << VCFfile->site()<<endl;
    double P1 = VCFfile->site();
    double p1 = VCFfile->seg_end_site();

    bool previous_segment_is_invariant = ( VCFfile->previous_seg_state != MISSING ) ; // previous_segment_is_invariant could be ZERO_SEG or SEQ_INVARIANT
    bool current_segment_is_invariant = VCFfile->current_seg_state != MISSING; // current_segment_is_invariant could be ZERO_SEG or SEQ_INVARIANT
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
    
    /*!
     * P0-------p0------P1-----p1-----P2
     * 
     * VCF: p0 = P0, p1 = P1, P0, P1, and P2 are SNP
     * Previously updated till P0/p0, 
     * Action: extend from P0/p0 to P1, then update at P1.
     * 
     * GVCF: P0: SNP, P1: SNP, same as VCF
     *       P0: invariant, P1: SNP
     *           previouly updated till p0, for invariant between P0 and p0
     *           Action: extend missing from p0 to P1, then update at P1/p1
     *       P0: SNP, P1: invariant
     *           previously updated till P0/p0, 
     *           Action: extending invariant from P0 till p1, no update at P1 nor p1
     *       P0: invariant, P1: invariant
     *           previouly updated till p0, extend invariant till p1, no update at P1 nor p1
     * 
     * RGVCF: P0: SNP, P1: SNP, same as VCF
     *       P0: invariant, P1: SNP
     *           previouly updated till P0/p0, for missing between P0 and p0
     *           Action: extend invariant from p0 to P1, then update at P1/p1
     *       P0: SNP, P1: invariant
     *           previously updated till P0/p0, 
     *           Action: extending invariatn from P0/p0 till p1, no update at P1 nor p1
     *       P0: invariant, P1: invariant
     *           previouly updated till p0, extend missing till p1, no update at P1 nor p1
     * 
     */ 
    
    /*!
     * extend to P1, 
     * If P1 == p1
     *      update at p1
     * else 
     *      then extend to p1 (two types of extension, different from the previous extension)
     */ 
    //Extend ARGs and update weight for not seeing mutations along the equences
    cout << " extend ARG part I, previous_segment is " ;
    cout << (previous_segment_is_invariant? "":"Not" );
    cout << " invariant" <<endl;
    this->extend_ARGs( P1, mutation_rate, previous_segment_is_invariant );
    cout << " extend ARG part I finished "<<endl;
    if ( P1 == p1 ){
        //Update weight for seeing mutation at the position 
        cout << " Update state weight at a SNP "<<endl;
        this->update_state_weights_at_A_single_site( p1, mutation_rate, VCFfile->current_variant_state != SNP, VCFfile->vec_of_sample_alt_bool ); 
        }
    else {
        //Extend ARGs and update weight for not seeing mutations along the equences
    cout << " extend ARG part II "<<endl;
        this->extend_ARGs( p1, mutation_rate, current_segment_is_invariant );
        }
    
    //Update the cumulated probabilities, as well as computing the effective sample size
    this->update_cum_sum_array_find_ESS( weight_cum_sum );
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
        //this->particles[i]->setParticleWeight( this->random_generator()->sample() );
        this->particles[i]->setParticleWeight( this->particles[i]->random_generator()->sample() );
        }
    }

// We need to decide at the tail of the data, until the end of the sequence, whether to perform recombination or not, extend arg from the prior? or ?
void ParticleContainer::cumulate_recomb_opportunity_at_seq_end( double seqend ){
    for (size_t i = 0; i < this->particles.size(); i++){
        double opportunity_x = seqend - this->particles[i]->current_base();
        double opportunity_y = this->particles[i]->getLocalTreeLength();
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
