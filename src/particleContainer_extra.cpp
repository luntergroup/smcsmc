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

/*! 
 * @ingroup group_pf_update
 * \brief Update the current state to the next state, at the given site, update all particles to it's latest genealogy state.  Also include the likelihood for no mutations.
 */
void ParticleContainer::extend_ARGs( double mutation_rate, double extend_to, Segment_State segment_state ){
    dout << endl<<" We are extending particles" << endl<<endl;
    
 	for (size_t particle_i = 0; particle_i < this->particles.size(); particle_i++){
        dout << "We are updating particle " << particle_i << endl;
        /*! 
         * For each particle, extend the current path until the the site such that the next genealogy change is beyond the mutation
         * Invariant: the likelihood is correct up to 'updated_to'
         */
        double updated_to = this->particles[particle_i]->site_where_weight_was_updated();
        dout << "Particle current base is at " << this->particles[particle_i]->current_base() << " weight is updated to " << updated_to <<endl;
        assert (updated_to >= this->particles[particle_i]->current_base());
        while ( updated_to < extend_to ) {
            dout << "  Now at " <<this->particles[particle_i]->current_base()<< " updated_to " << updated_to << " and extending to " << extend_to << endl;            
            /*!
             * First, update the likelihood up to either extend_to or the end of this state
             */
            double update_to = min( extend_to, this->particles[particle_i]->next_base() );
            double length_of_local_tree = this->particles[particle_i]->getLocalTreeLength(); // in generations
            double likelihood_of_segment = ( segment_state == SEGMENT_INVARIANT ) ? exp( -mutation_rate * length_of_local_tree * (update_to - updated_to) ) : 1 ;// assume infinite site model
            dout << " Likelihood of no mutations in segment of length " << (update_to - updated_to) << " is " << likelihood_of_segment ;
            dout << ( ( segment_state == SEGMENT_INVARIANT ) ? ", as invariant.": ", as missing data" ) << endl;
            this->particles[particle_i]->setParticleWeight( this->particles[particle_i]->weight() * likelihood_of_segment);
            updated_to = update_to;  // rescues the invariant
            /*!
             * Next, if we haven't reached extend_to now, add a new state and iterate
             */
            if ( updated_to < extend_to ) {
                this->particles[particle_i]->sampleNextGenealogy();
                if ( this->heat_bool_ ){
                    TmrcaState tmrca( this->particles[particle_i]->site_where_weight_was_updated(), this->particles[particle_i]->local_root()->height() );
                    this->particles[particle_i]->TmrcaHistory.push_back ( tmrca );
                    }
                
                }
            
            }
        assert (updated_to == extend_to);        
        this->particles[particle_i]->setSiteWhereWeightWasUpdated( extend_to );
        }
    /*! normalize the probability upon until the mutation */
    //this->normalize_probability(); // This normalization doesn't seem to do much ...
    
    }



/*! \brief Update particle weight according to the haplotype data
 *	@ingroup group_pf_update
 */
void ParticleContainer::update_weight_at_site( double mutation_rate, vector <int> &haplotypes_at_tips ){
			
	// now update the weights of all particles, by calculating the likelihood of the data over the previous segment	
	for (size_t particle_i=0; particle_i < this->particles.size(); particle_i++){
		this->particles[particle_i]->include_haplotypes_at_tips(haplotypes_at_tips);

        double likelihood_of_haplotypes_at_tips = this->particles[particle_i]->calculate_likelihood( );
        dout << "updated weight =" << this->particles[particle_i]->weight()  << "*" <<  likelihood_of_haplotypes_at_tips <<endl;

        this->particles[particle_i]->setParticleWeight( this->particles[particle_i]->weight() * likelihood_of_haplotypes_at_tips);
		dout << "particle " <<  particle_i<<" done" << endl;
        }
    
    this->normalize_probability(); // It seems to converge slower if it is not normalized ...
	dout << endl;
    }


void ParticleContainer::duplicate_particles ( valarray<int> & sample_count ){
    for ( size_t i = 0 ;  i < this->particles.size(); i++ ){
        this->particles[i]->making_copies( sample_count[i] ); 
        }
    }
