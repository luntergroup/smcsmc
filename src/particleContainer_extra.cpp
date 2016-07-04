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
        dout << " before extend ARG particle weight is " << this->particles[particle_i]->weight() << endl;
        (void)this->particles[particle_i]->extend_ARG ( mutation_rate, extend_to, segment_state);
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


/*! \brief Update particle weight according to the haplotype or genotype data
 *	@ingroup group_pf_update
 */
void ParticleContainer::update_weight_at_site( double mutation_rate, const vector <int> &data_at_tips ){
				
	// first check if there are any (unphased) genotypes in the data
	vector<int> haplotype_at_tips = data_at_tips;
	int num_haplotypes = 1;
	for (int i=0; i < data_at_tips.size(); i+=2) {
		if (data_at_tips[i] == 2) {
			num_haplotypes *= 2;
			assert (data_at_tips[i+1] == 2);
			haplotype_at_tips[i] = 0;    // initialize phase vector
			haplotype_at_tips[i+1] = 1;
		} else {
			assert (data_at_tips[i+1] != 2);
		}
	}
	double normalization_factor = 1.0 / num_haplotypes;

	// now update the weights of all particles, by calculating the likelihood of the data over the previous segment	
	for (size_t particle_i = 0; particle_i < particles.size(); particle_i++) {
		double likelihood_of_haplotype_at_tips = 0;
		do {
			particles[particle_i]->include_haplotypes_at_tips( haplotype_at_tips );
			likelihood_of_haplotype_at_tips += particles[particle_i]->calculate_likelihood( );
		} while ((num_haplotypes > 1) && next_haplotype(haplotype_at_tips, data_at_tips));
		likelihood_of_haplotype_at_tips *= normalization_factor;
        dout << "updated weight =" << particles[particle_i]->weight() << "*" << likelihood_of_haplotype_at_tips << "*" << particles[particle_i]->importance_weight_predata() << endl;
        // the following assertion checks appropriate importance factors have been accounted for
        //   if we ever change the emission factor of the sampling distribution, this assertion should fail
        assert( particles[particle_i]->importance_weight_predata() == 1 );
        this->particles[particle_i]->setParticleWeight( particles[particle_i]->weight() 
                                                      * likelihood_of_haplotype_at_tips
                                                      * particles[particle_i]->importance_weight_predata() );
        this->particles[particle_i]->setDelayedWeight( particles[particle_i]->delayed_weight() 
                                                      * likelihood_of_haplotypes_at_tips
                                                      * particles[particle_i]->importance_weight_predata() );
        this->particles[particle_i]->reset_importance_weight_predata();
		dout << "particle " << particle_i << " done" << endl;
	}    
    
    this->store_normalization_factor();
    this->normalize_probability(); // It seems to converge slower if it is not normalized ...
	dout << endl;
}
    
void ParticleContainer::store_normalization_factor() {
	temp_sum_of_weights = 0;
	for( size_t i = 0; i < this->particles.size(); i++){
		temp_sum_of_weights += this->particles[i]->weight();
    }
	ln_normalization_factor_ += log( temp_sum_of_weights );
}
