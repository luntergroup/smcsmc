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
#include <omp.h>

//#ifdef _OPENMP
    //#include <omp.h>
    //#define TIMESCALE 1
//#else
    //#define omp_get_thread_num() 0
    //#define omp_get_num_procs() 0
    //#define omp_get_num_threads() 1
    //#define omp_set_num_threads(bob) 0
    //#define omp_get_wtime() clock()
    //#define TIMESCALE CLOCKS_PER_SEC
//#endif
/*! 
 * @ingroup group_pf_update
 * \brief Update the current state to the next state, at the given site, update all particles to it's latest genealogy state.  Also include the likelihood for no mutations.
 */
void ParticleContainer::extend_ARGs( double mutation_at, double mutation_rate, bool withdata ){
    dout << endl<<" We are extending particles" << endl<<endl;    
    // USE MULTITHREADING ...

    //#pragma omp for schedule(dynamic,100) nowait    
    #pragma omp parallel for schedule(dynamic) 
 	for (size_t particle_i=0; particle_i < this->particles.size(); particle_i++){
        
        //cout << "Hello World from thread " << omp_get_thread_num() << '\n';
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
    #pragma omp parallel for schedule(dynamic) 
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


void ParticleContainer::duplicate_particles ( valarray<int> & sample_count ){
    #pragma omp parallel for schedule(dynamic) 
    for ( size_t i = 0 ;  i < this->particles.size(); i++ ){
        this->particles[i]->making_copies( sample_count[i] ); 
        }
    }
