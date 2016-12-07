/*
 * smcsmc is short for particle filters for ancestral recombination graphs.
 * This is a free software for demographic inference from genome data with particle filters.
 *
 * Copyright (C) 2013, 2014 Sha (Joe) Zhu and Gerton Lunter
 *
 * This file is part of smcsmc.
 *
 * smcsmc is free software: you can redistribute it and/or modify
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
ParticleContainer::ParticleContainer(Model* model,
                                     MersenneTwister *rg,
                                     const vector<int>& record_event_in_epoch,
                                     size_t Num_of_states,
                                     double initial_position,
                                     bool emptyFile,
                                     vector <int> first_allelic_state) {
    this->random_generator_ = rg;
    this->set_ESS(0);
    this->ln_normalization_factor_ = 0;
    this->set_current_printing_base(0);
    this->model = model;
    bool make_copies_of_model_and_rg = false; // Set to true to use a random generator, and model per particle for multithreading (slower!)
    dout << " --------------------   Particle Initial States   --------------------" << std::endl;
    for ( size_t i=0; i < Num_of_states ; i++ ){
        ForestState* new_state = new ForestState( model, rg, record_event_in_epoch, make_copies_of_model_and_rg );  // create a new state, using scrm; scrm always starts at 0.
        // Initialize members of FroestState (derived class members)
        new_state->init_EventContainers( model );
        new_state->buildInitialTree();
        new_state->setSiteWhereWeightWasUpdated( initial_position );
        new_state->setParticleWeight( 1.0/Num_of_states );
        new_state->setDelayedWeight( 1.0/Num_of_states );
        this->particles.push_back(new_state);
        // If no data was given, the initial tree should not include any data
        if ( emptyFile ){
            new_state->include_haplotypes_at_tips(first_allelic_state);
        }
        #ifdef _SCRM
        cout << new_state->newick( new_state->local_root() ) << ";" << endl;
        #endif
    }
}


/*!
 * @ingroup group_pf_update
 * \brief Update the current state to the next state, at the given site, update all particles to it's latest genealogy state.  Also include the likelihood for no mutations.
 */
void ParticleContainer::extend_ARGs( double mutation_rate, double extend_to ){
    dout << endl<<" We are extending particles" << endl<<endl;

        for (size_t particle_i = 0; particle_i < this->particles.size(); particle_i++){
        dout << "We are updating particle " << particle_i << endl;
        /*!
         * For each particle, extend the current path until the the site such that the next genealogy change is beyond the mutation
         * Invariant: the likelihood is correct up to 'updated_to'
         */
        dout << " before extend ARG particle weight is " << this->particles[particle_i]->weight() << endl;
        particles[particle_i]->extend_ARG ( mutation_rate, extend_to );
        dout << " after extend ARG particle weight is " << this->particles[particle_i]->weight() << endl;
    }
    /*! normalize the probability upon until the mutation */
    //this->normalize_probability(); // This normalization doesn't seem to do much ...
}


bool ParticleContainer::next_haplotype( vector<int>& haplotype_at_tips, const vector<int>& data_at_tips ) const {

        // find the next haplotype that is compatible with the genotype encoded in data_at_tips
        for (int i=0; i<haplotype_at_tips.size(); i+=2) {

                if (data_at_tips[i] != 2) continue;   // phased site -- ignore
                if (haplotype_at_tips[i] == 0) {      // phase 0 -- change to phase 1, and done
                        haplotype_at_tips[i] = 1;
                        haplotype_at_tips[i+1] = 0;
                        return true;
                }
                haplotype_at_tips[i] = 0;             // phase 1 -- flip back, and continue
                haplotype_at_tips[i+1] = 1;

        }
        // we fell through, so we have considered all haplotypes; stop main loop
        return false;
}


void ParticleContainer::update_data_status_at_leaf_nodes( const vector<int>& data_at_tips ) {

    vector<int> haplotype_at_tips;
    calculate_initial_haplotype_configuration( data_at_tips, haplotype_at_tips );

    for (size_t particle_i = 0; particle_i < particles.size(); particle_i++) {

        particles[particle_i]->include_haplotypes_at_tips( haplotype_at_tips );

    }
}


int ParticleContainer::calculate_initial_haplotype_configuration( const vector<int>& data_at_tips, vector<int>& haplotype_at_tips ) const {

    // just copy any phased haplotype or homozygous genotype data across -- as this is provided in pairs, this is
    // automatically correct haplotype data.  For unphased het sites, which are encoded as pairs of '2's, assign
    // a 0/1 configuration arbitrarily.

    haplotype_at_tips = data_at_tips;
    int num_configurations = 1;
    for (int i=0; i < data_at_tips.size(); i+=2) {
        if (data_at_tips[i] == 2) {
            num_configurations *= 2;
            assert (data_at_tips[i+1] == 2);
            haplotype_at_tips[i] = 0;    // initialize phase vector
            haplotype_at_tips[i+1] = 1;
        } else {
            assert (data_at_tips[i+1] != 2);
        }
    }
    return num_configurations;
}


/*! \brief Update particle weight according to the haplotype or genotype data
 *      @ingroup group_pf_update
 */
void ParticleContainer::update_weight_at_site( double mutation_rate, const vector <int> &data_at_tips ){

    vector<int> haplotype_at_tips;
    int num_configurations = calculate_initial_haplotype_configuration( data_at_tips, haplotype_at_tips );
    double normalization_factor = 1.0 / num_configurations;

    // now update the weights of all particles, by calculating the likelihood of the data over the previous segment
    for (size_t particle_i = 0; particle_i < particles.size(); particle_i++) {
        double likelihood_of_haplotype_at_tips = 0;
        do {

            particles[particle_i]->include_haplotypes_at_tips( haplotype_at_tips );
            likelihood_of_haplotype_at_tips += particles[particle_i]->calculate_likelihood( );

        } while ((num_configurations > 1) && next_haplotype(haplotype_at_tips, data_at_tips));

        likelihood_of_haplotype_at_tips *= normalization_factor;

        dout << "updated weight =" << particles[particle_i]->weight() << "*" << likelihood_of_haplotype_at_tips << endl;
        // the following assertion checks appropriate importance factors have been accounted for
        //   if we ever change the emission factor of the sampling distribution, this assertion should fail
        this->particles[particle_i]->setParticleWeight( particles[particle_i]->weight()
                                                        * likelihood_of_haplotype_at_tips );
        this->particles[particle_i]->setDelayedWeight( particles[particle_i]->delayed_weight()
                                                       * likelihood_of_haplotype_at_tips );
        dout << "particle " << particle_i << " done" << endl;
    }

    this->store_normalization_factor();
    this->normalize_probability(); // It seems to converge slower if it is not normalized ...
    dout << endl;
}


/*! \brief Resampling step
 *  If the effective sample size is less than the ESS threshold, do a resample, currently using systemetic resampling scheme.
 */
void ParticleContainer::ESS_resampling(valarray<double> weight_cum_sum, valarray<int> &sample_count, int mutation_at, const PfParam &pfparam, int num_state)
{
    dout << "At pos " << mutation_at << " ESS is " <<  this->ESS() <<", number of particle is " <<  num_state << ", and ESSthreshold is " << pfparam.ESSthreshold <<endl;
    double ESS_diff = pfparam.ESSthreshold - this->ESS();
    if ( ESS_diff > 1e-6 ) { // resample if the effective sample size is small, to check this step, turn the if statement off
        resampledout<<" ESS_diff = " << ESS_diff<<endl;
        resampledout << " ### PROGRESS: ESS resampling" << endl;
        pfparam.append_resample_file( mutation_at, ESS() );
        this->systematic_resampling( weight_cum_sum, sample_count, num_state);
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
    dout << endl;
    resampledout << " Recreate particles" << endl;
    resampledout << " ****************************** Start making list of new states ****************************** " << std::endl;
    resampledout << " will make total of " << sample_count.sum()<<" particle states" << endl;
    size_t number_of_particles = sample_count.size();
    bool flush = this->particles[0]->segment_count() > 50;  // keep Forest::rec_bases_ vector down to reasonable size
    // calculate the sum of the delayed weights, this should be invariant
    double sum_of_delayed_weights = 0;
    if( model->biased_sampling ) {
        for (size_t old_state_index = 0; old_state_index < number_of_particles; old_state_index++) {
            sum_of_delayed_weights += particles[old_state_index]->delayed_weight();
        }
    }
    for (size_t old_state_index = 0; old_state_index < number_of_particles; old_state_index++) {
        if ( sample_count[old_state_index] > 0 ) {
            ForestState * current_state = this->particles[old_state_index];
            if (flush) current_state->flushOldRecombinations();
            // we need at least one copy of this particle; it keeps its own random generator
            resampledout << " Keeping  the " << std::setw(5) << old_state_index << "th particle" << endl;
            if( !model->biased_sampling ){
                            current_state->setParticleWeight( 1.0/number_of_particles );
                        } else {
                assert( sum_of_delayed_weights != 0 );
                            current_state->setDelayedWeight( 1.0/number_of_particles * sum_of_delayed_weights );
                            current_state->setParticleWeight( current_state->delayed_weight()*current_state->total_delayed_adjustment );
                    }
            this->particles.push_back(current_state);
            // create new copy of the resampled particle
            for (int ii = 2; ii <= sample_count[old_state_index]; ii++) {
                resampledout << " Making a copy of the " << old_state_index << "th particle ... " ;
                dout << std::endl;
                dout << "current end position for particle current_state " << current_state->current_base() << endl;
                for (size_t ii =0 ; ii < current_state ->rec_bases_.size(); ii++){
                    dout << current_state ->rec_bases_[ii] << " ";
                }
                dout <<endl;

                ForestState* new_copy_state = new ForestState( *this->particles[old_state_index] );
                dout <<"making particle finished" << endl; // DEBUG

                // Resample the recombination position if particle has not hit end of sequence, and give particle its own event history
                if ( new_copy_state->current_base() < new_copy_state->next_base() ){
                    new_copy_state->resample_recombination_position();
                }
                if( !model->biased_sampling ) {
                                    new_copy_state->setParticleWeight( 1.0/number_of_particles );
                                } else {
                    assert( sum_of_delayed_weights != 0 );
                                    new_copy_state->setDelayedWeight( 1.0/number_of_particles * sum_of_delayed_weights );
                                new_copy_state->setParticleWeight( new_copy_state->delayed_weight()*new_copy_state->total_delayed_adjustment );
                }
                this->particles.push_back(new_copy_state);
            }
        } else {
            resampledout << " Deleting the " << std::setw(5) << old_state_index << "th particle ... " ;
            delete this->particles[old_state_index];
        }

        this->particles[old_state_index]=NULL;
    }

    for (int i = 0; i < number_of_particles; i++) {
        this->particles[i] = this->particles[i + number_of_particles];
    }
    this -> particles.resize(number_of_particles);

    resampledout << "There are " << this->particles.size() << "particles in total." << endl;
    resampledout << " ****************************** End of making list of new particles ****************************** " << std::endl;
    }



/*!
 * ParticleContatiner destructor
 */
ParticleContainer::~ParticleContainer() {}


/*!
 * Proper particleContatiner destructor, remove pointers of the ForestState
 */
void ParticleContainer::clear(){
    for (size_t i = 0; i < this->particles.size(); i++){
        if (this->particles[i]!=NULL){
            delete this->particles[i];
            this->particles[i]=NULL;
        }
    }
    this->particles.clear();
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
        double w_i = model->biased_sampling ? this->particles[i]->delayed_weight() : this->particles[i]->weight();
        weight_cum_sum[i+1]=weight_cum_sum[i]+w_i;
        wi_sum = wi_sum + w_i;
        wi_sq_sum = wi_sq_sum + w_i * w_i;
        }

    //check for the cum weight
    if( !model->biased_sampling ){
            dout << "### particle weights ";
            for (size_t i=0;i<Num_of_states;i++){
                dout << this->particles[i]->weight()<<"  ";
                } dout << std::endl<<std::endl;

            dout << "### updated cum sum of particle weight ";
            for (size_t i=0;i<weight_cum_sum.size();i++){
                dout << weight_cum_sum[i]<<"  ";
                } dout << std::endl;
        } else {
            dout << "### delayed particle weights ";
            for (size_t i=0;i<Num_of_states;i++){
                dout << this->particles[i]->delayed_weight()<<"  ";
                } dout << std::endl<<std::endl;

            dout << "### updated cum sum of delayed particle weight ";
            for (size_t i=0;i<weight_cum_sum.size();i++){
                dout << weight_cum_sum[i]<<"  ";
                } dout << std::endl;
        }

    this->set_ESS(wi_sum * wi_sum / wi_sq_sum);
    }


/*!
 * Normalize the particle weight, inorder to prevent underflow problem
 */
void ParticleContainer::normalize_probability() {

    double total_probability = 0;
    for ( size_t particle_i = 0;particle_i < this->particles.size(); particle_i++ ) {
        total_probability += this->particles[particle_i]->weight();
    }
    for ( size_t particle_i = 0; particle_i < this->particles.size(); particle_i++ ) {
        this->particles[particle_i]->setParticleWeight( this->particles[particle_i]->weight() / total_probability);
        this->particles[particle_i]->setDelayedWeight(this->particles[particle_i]->delayed_weight() / total_probability); // should be conditional on biased_sampling
    }
}


void ParticleContainer::update_state_to_data( double mutation_rate, double loci_length, Segment * Segfile, valarray<double> & weight_cum_sum ){

    dout <<  " ******************** Update the weight of the particles  ********** " <<endl;
    dout << " ### PROGRESS: update weight at " << Segfile->segment_start()<<endl;

    // Assign presence/absence status of data to each of the leaf nodes of all the particles
    update_data_status_at_leaf_nodes( Segfile->allelic_state_at_Segment_end );

    //Extend ARGs and update weight for not seeing mutations along the sequences
    this->extend_ARGs( mutation_rate, (double)min(Segfile->segment_end(), loci_length) );

    //Update weight for seeing mutation at the position
    dout << " Update state weight at a SNP "<<endl;
    if (Segfile->segment_state() == SEGMENT_INVARIANT) {
        // ensure that if this segment is a partial segment, and does not end in a mutation
        // even though data is not missing (state == SEGMENT_INVARIANT_PARTIAL), the mutation
        // is not accounted for here, but only when the last partial segment is processed.
        this->update_weight_at_site( mutation_rate, Segfile->allelic_state_at_Segment_end );
    }

    dout << "Extended until " << this->particles[0]->current_base() <<endl;

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
        assert( sample_count.sum() == num_state );
    }


/*!
 * @ingroup group_systematic
 * \brief Use systematic resampling \cite Doucet2008 to generate sample count for each particle
 */
void ParticleContainer::systematic_resampling(std::valarray<double> cum_sum, std::valarray<int>& sample_count, int sample_size){
    size_t interval_j = 0;
    size_t sample_i = 0;
    size_t N = sample_size;
    //double u_j = rand() / double(RAND_MAX) / N;
    double u_j = this->random_generator()->sample() / N;
    double cumsum_normalization = cum_sum[cum_sum.size()-1];

    resampledout << "systematic sampling procedure on interval:" << std::endl;
    resampledout << " ";
    for (size_t i=0;i<cum_sum.size();i++){dout <<  (cum_sum[i]/cumsum_normalization )<<"  ";}dout << std::endl;

    sample_count[sample_i] = 0;
    while (sample_i < N) {
        resampledout << "Is " <<  u_j<<" in the interval of " << std::setw(10)<< (cum_sum[interval_j]/ cumsum_normalization) << " and " << std::setw(10)<< (cum_sum[interval_j+1]/ cumsum_normalization) << " ? ";
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

    resampledout << "systematic sampling procedue finished with total sample count " << sample_count.sum()<<std::endl;
    resampledout << "Sample counts: " ;
    for (size_t i=0;i<sample_count.size();i++){dout << sample_count[i]<<"  ";}  dout << std::endl;
    assert(sample_count.sum() == sample_size);
    }


// set_particles_with_random_weight() has not been updated for use with delayed IS
void ParticleContainer::set_particles_with_random_weight(){
    for (size_t i = 0; i < this->particles.size(); i++){
        this->particles[i]->setParticleWeight( this->particles[i]->random_generator()->sample() );
        }
    }


void ParticleContainer::print_particle_probabilities(){
    for (size_t i = 0; i < this->particles.size(); i++){
        clog<<"weight = "<<this->particles[i]->weight()<<endl;
        if( model->biased_sampling ) {clog<<"delayed weight = "<<this->particles[i]->delayed_weight()<<endl;}
        }
    }


void ParticleContainer::store_normalization_factor() {
        temp_sum_of_weights = 0;
        for( size_t i = 0; i < this->particles.size(); i++){
                temp_sum_of_weights += this->particles[i]->weight();
    }
        ln_normalization_factor_ += log( temp_sum_of_weights );
}


void ParticleContainer::print_ln_normalization_factor(){
        clog << "Our likelihood measure is " << ln_normalization_factor() << endl;
}
