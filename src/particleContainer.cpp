/*
 * smcsmc is short for particle filters for ancestral recombination graphs.
 * This is a free software for demographic inference from genome data with particle filters.
 *
 * Copyright (C) 2013-2017 Donna Henderson, Sha (Joe) Zhu and Gerton Lunter
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
ParticleContainer::ParticleContainer(Model *model,
                                     MersenneTwister *rg,
                                     PfParam *pfparam,
                                     size_t Num_of_states,
                                     double initial_position,
                                     bool emptyFile,
                                     vector <int> first_allelic_state) {
    random_generator_ = rg;
    ln_normalization_factor_ = 0;
    this->model = model;
    bool make_copies_of_model_and_rg = false; // Set to true to use rg and model per particle, for multithreading
    dout << " --------------------   Particle Initial States   --------------------" << std::endl;
    for ( size_t i=0; i < Num_of_states ; i++ ) {

        // create a new state, using scrm; scrm always starts at site 0.
        ForestState* new_state = new ForestState( model, rg, pfparam, make_copies_of_model_and_rg );  

        // Initialize members of FroestState (derived class members)
        new_state->init_EventContainers( model );
        new_state->buildInitialTree( true );
        new_state->record_recomb_extension();
        new_state->save_recomb_state();
        new_state->setSiteWhereWeightWasUpdated( initial_position );
        new_state->setPosteriorWeight( 1.0/Num_of_states );
        new_state->setPilotWeight( 1.0/Num_of_states );
        particles.push_back(new_state);

        #ifdef _SCRM
        cout << new_state->newick( new_state->local_root() ) << ";" << endl;
        #endif
    }
}

void ParticleContainer::print_recent_recombination_histogram() {

    /*
    vector<int> hist( 1 );
    cout << "======= Recent recombination count histogram: =========" << endl;
    for (size_t i=0; i<particles.size(); i++) {
        size_t count = particles[i]->recent_recombination_count;
        if (count >= hist.size()) hist.resize(count+1);
        hist[count]++;
    }
    for (size_t i=0; i<hist.size(); i++) {
        cout << i << "\t" << hist[i] << endl;
    }
    */
    double rec_acc = 0, rec_wt = 0;
    for (size_t i=0; i<particles.size(); i++) {
        rec_acc += particles[i]->recent_recombination_count * particles[i]->posteriorWeight();
        rec_wt  += particles[i]->posteriorWeight();
    }
    cout << "Posterior mean recent recombination count: " << rec_acc / rec_wt << endl;
    
}


/*!
 * @ingroup group_pf_update
 * \brief Update the current state to the next state, at the given site, update all particles to its latest genealogy state.  
 * \brief Also include the likelihood for no mutations.
 */
void ParticleContainer::extend_ARGs( double extend_to, const vector<int>& data_at_site ){
    dout << endl<<" We are extending particles" << endl<<endl;

    // Calculate number of leaves with missing data (for speeding up the likelihood calculation)
    int num_leaves_missing = 0;
    int leaf_status = 0;   // -1 all missing; 1 all present; 0 mixed
    for (int i=0; i<data_at_site.size(); i++) {
        num_leaves_missing += data_at_site[i] == -1;
    }
    if (num_leaves_missing == 0)
        leaf_status = 1;
    if (num_leaves_missing == data_at_site.size())
        leaf_status = -1;

    
    for (size_t particle_i = 0; particle_i < particles.size(); particle_i++){
        dout << "We are updating particle " << particle_i << endl;
        /*!
         * For each particle, extend the current path until the the site such that 
         * the next genealogy change is beyond the mutation
         * Invariant: the likelihood is correct up to 'updated_to'
         */
        dout << " before extend ARG particle weight is " << particles[particle_i]->posteriorWeight() << endl;
        particles[particle_i]->restore_recomb_state();
        particles[particle_i]->extend_ARG ( extend_to, leaf_status, data_at_site );
        particles[particle_i]->save_recomb_state();
        dout << " after extend ARG particle weight is " << particles[particle_i]->posteriorWeight() << endl;
    }
}


int ParticleContainer::calculate_initial_haplotype_configuration( const vector<int>& data_at_tips,
                                                                  vector<int>& haplotype_at_tips ) const {

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


/*! \brief Update particle weight according to the haplotype or genotype data
 *      @ingroup group_pf_update
 */
void ParticleContainer::update_weight_at_site( const vector <int> &data_at_tips, bool ancestral_aware ){

    vector<int> haplotype_at_tips;
    int num_configurations = calculate_initial_haplotype_configuration( data_at_tips, haplotype_at_tips );
    double normalization_factor = 1.0 / num_configurations;

    // now update the weights of all particles, by calculating the likelihood of the data over the previous segment
    for (size_t particle_i = 0; particle_i < particles.size(); particle_i++) {
        double likelihood_of_haplotype_at_tips = 0;
        do {

            likelihood_of_haplotype_at_tips += particles[particle_i]->calculate_likelihood( ancestral_aware, haplotype_at_tips );

        } while ((num_configurations > 1) && next_haplotype(haplotype_at_tips, data_at_tips));

        likelihood_of_haplotype_at_tips *= normalization_factor;

        // the following assertion checks appropriate importance factors have been accounted for
        //   if we ever change the emission factor of the sampling distribution, this assertion should fail
        this->particles[particle_i]->adjustWeights( likelihood_of_haplotype_at_tips );
    }
}


/*! \brief Resampling step
 *  If the effective sample size is less than the ESS threshold, do a resample, currently using systemetic resampling scheme.
 */
int ParticleContainer::resample(int update_to, const PfParam &pfparam)
{
    double Num_of_states = particles.size();
    valarray<double> partial_sums( Num_of_states + 1 );
    double wi_sum = 0;
    double wi_sq_sum = 0;

    //partial_sums[0] = 0;
    for (size_t i = 0; i < Num_of_states; i++) {
        partial_sums[i] = wi_sum;
        double w_i = particles[i]->pilotWeight();
        //partial_sums[i+1] = partial_sums[i] + w_i;
        wi_sum += w_i;
        wi_sq_sum += w_i * w_i;
    }
    partial_sums[ Num_of_states ] = wi_sum;
    double ess = wi_sum * wi_sum / wi_sq_sum;
    
    dout << "At pos " << update_to << " ESS is " << ess <<", number of particles is " << Num_of_states
         << " threshold is " << pfparam.ESSthreshold << " diff=" << pfparam.ESSthreshold - ess << endl;
    
    if ( ess < pfparam.ESSthreshold - 1e-6 ) {
        pfparam.append_resample_file( update_to, ess );
        valarray<int> sample_count( Num_of_states );
        systematic_resampling( partial_sums, sample_count, Num_of_states );
        implement_resampling( sample_count, wi_sum);
        return 1;
    } else {
        // The sampling procedure will flush old recombinations, but this won't happen in the absence of data.
        // Ensure that even in the absence of data, old recombinations do not build up and cause memory overflows
        bool flush = particles[0]->segment_count() > 50;
        if (flush) {
            for (size_t state_index = 0; state_index < particles.size(); state_index++) {
                particles[state_index]->flushOldRecombinations();
            }
        }
        return 0;
    }
}



/*! \brief Particle filtering Resampling step
 *
 *  Update particles according to the sample counts, and set the particle probabilities to 1
 *
 * \ingroup group_resample
 */
void ParticleContainer::implement_resampling(valarray<int> & sample_count, double sum_of_delayed_weights) {

    size_t number_of_particles = sample_count.size();
    bool flush = this->particles[0]->segment_count() > 50;  // keep Forest::rec_bases_ vector down to reasonable size

    for (size_t old_state_index = 0; old_state_index < number_of_particles; old_state_index++) {

        if ( sample_count[old_state_index] > 0 ) {

            ForestState * current_state = particles[old_state_index];
            if (flush) current_state->flushOldRecombinations();

            // calculate weight adjustment: the ratio of the desired weight, sum_of_delayed_weights / num_particles,
            // to current weight, pilotWeight()
            double adjustment = sum_of_delayed_weights / ( number_of_particles * current_state->pilotWeight() );
            current_state->adjustWeights( adjustment );
            particles.push_back(current_state);

            // create new copy of the resampled particle
            for (int ii = 1; ii < sample_count[old_state_index]; ii++) {

                ForestState* new_copy_state = new ForestState( *particles[old_state_index] );

                // Resample the recombination position if particle has not hit end of sequence,
                // and give particle its own event history
                if ( new_copy_state->current_base() < new_copy_state->next_base() ) {
                    new_copy_state->restore_recomb_state();
                    new_copy_state->resample_recombination_position();
                    new_copy_state->save_recomb_state();
                }

                // calculate weight adjustment. Since the old particle's weight was changed in-place above,
                // the weight adjustment is 1.0 now (if the resampling of recombination position above did not
                // change its pilot weight -- which it doesn't.)
                /*
                double adjustment = sum_of_delayed_weights / ( number_of_particles * new_copy_state->pilotWeight() );
                new_copy_state->adjustWeights( adjustment );
                */

                particles.push_back(new_copy_state);
            }
        } else {
            delete particles[old_state_index];
        }

        particles[old_state_index]=NULL;
    }

    for (int i = 0; i < number_of_particles; i++) {
        particles[i] = particles[i + number_of_particles];
    }
    particles.resize(number_of_particles);
}



/*!
 * ParticleContatiner destructor
 */
ParticleContainer::~ParticleContainer() {}


/*!
 * Proper particleContatiner destructor, remove pointers of the ForestState
 */
void ParticleContainer::clear(){
    for (size_t i = 0; i < particles.size(); i++){
        if (particles[i]!=NULL){
            delete particles[i];
            particles[i]=NULL;
        }
    }
    particles.clear();
}



/*!
 * Normalize the particle weight, inorder to prevent underflow problem
 */
void ParticleContainer::normalize_probability() {

    double total_probability = 0;
    for ( size_t particle_i = 0; particle_i < particles.size(); particle_i++ )
        total_probability += particles[particle_i]->posteriorWeight();

    if (total_probability <= 0) throw std::runtime_error("Zero or negative probabilities");
    
    // keep track of overall (log) posterior probability
    ln_normalization_factor_ += log( total_probability );

    // normalize
    for ( size_t particle_i = 0; particle_i < particles.size(); particle_i++ ) {
        particles[particle_i]->adjustWeights( 1.0 / total_probability );
    }
}


void ParticleContainer::update_state_to_data( Segment * Segfile, bool ancestral_aware ) {

    dout <<  " ******************** Update the weight of the particles  ********** " <<endl;
    dout << " ### PROGRESS: update weight at " << Segfile->segment_start()<<endl;

    // Allelic state; -1 = absent, 0,1 are ref/alt, 2 = part of het site
    const vector<int>& data_at_site = Segfile->allelic_state_at_Segment_end;

    //Extend ARGs and update weight for not seeing mutations along the sequences
    extend_ARGs( (double)min(Segfile->segment_end(), (double)model->loci_length() ), data_at_site );

    //Update weight for seeing mutation at the position
    dout << " Update state weight at a SNP "<<endl;
    if (Segfile->segment_state() == SEGMENT_INVARIANT) {
        // account for the mutation at the end of an INVARIANT segment (which ends in a mutation).
        // Do not account for a mutation in an INVARIANT_PARTIAL segment (which is an initial fragment
        // of an INVARIANT segment, and therefore does not end with a mutation).  Also don't account
        // for a mutation in a MISSING segment.
        update_weight_at_site( Segfile->allelic_state_at_Segment_end, ancestral_aware );

        // TODO: Consider adding a column to the segfile which specifies if the ancestral allele is known for this
        //   position. This would replace the blanket ancestral_aware with a Segfile->ancestral_aware.
        //   We could also include a likelihood of 0 being ancestral.
    }

    dout << "Extended until " << particles[0]->current_base() <<endl;

    // Update the accumulated probabilities, as well as computing the effective sample size
    normalize_probability();
}



/*!
 * @ingroup group_systematic
 * \brief Use systematic resampling \cite Doucet2008 to generate sample count for each particle
 */
void ParticleContainer::systematic_resampling(std::valarray<double>& partial_sum,
                                              std::valarray<int>& sample_count,
                                              size_t N){
    size_t interval_j = 0;
    size_t sample_i = 0;
    double u_j = random_generator()->sample() / N;
    double partial_sum_normalization = partial_sum[ N ];

    sample_count[sample_i] = 0;
    while (sample_i < N) {

        /* invariants: */
        assert( (partial_sum[interval_j] / partial_sum_normalization) < u_j );
        assert( sample_i < N );

        if ( (sample_i == N) || partial_sum[interval_j+1] / partial_sum_normalization > u_j ) {
            sample_count[interval_j] += 1;
            sample_i += 1;
            u_j += 1.0/double(N);
        } else {
            interval_j += 1;
            sample_count[ interval_j ] = 0;
        }
    }
    interval_j=interval_j+1;
    for (;interval_j<N;interval_j++){
        sample_count[ interval_j ] = 0;
    }
}


// set_particles_with_random_weight() has not been updated for use with delayed IS
void ParticleContainer::set_particles_with_random_weight(){
    for (size_t i = 0; i < particles.size(); i++){
        particles[i]->setPosteriorWeight( particles[i]->random_generator()->sample() );
    }
}

