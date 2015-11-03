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
        (void)this->particles[particle_i]->extend_ARG ( mutation_rate, extend_to, segment_state);
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

        this->particles[particle_i]->setParticleWeight( this->particles[particle_i]->weight() 
                                                      * likelihood_of_haplotypes_at_tips
                                                      * this->particles[particle_i]->importance_weight_predata());
        this->particles[particle_i]->reset_importance_weight_predata();
		dout << "particle " <<  particle_i<<" done" << endl;
        }
    
    this->normalize_probability(); // It seems to converge slower if it is not normalized ...
	dout << endl;
    }
