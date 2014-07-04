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

#include"count.hpp"
void CountModel::extract_and_update_count(ParticleContainer &Endparticles, double current_base, bool end_data ) {
    // loop over all epochs
    for (size_t epoch_idx = 0; epoch_idx < this->change_times_.size(); epoch_idx++) {

        // calculate the required lagging for this epoch; don't use lagging for the final interval
        double lagging = end_data ? 0 : lags[epoch_idx];
        double x_end = current_base - lagging;

//cout << "["<<epoch_idx<< "] " << this->change_times_[epoch_idx] << " x_end = " << x_end<<endl;
		// Check that we're updating over more than minimal_lag_update_ratio * lagging nucleotides.
		// (If this is the last update, lagging will be 0, and we will do the update)
        if ( (x_end - this->counted_to[epoch_idx]) <= lagging * this->const_minimal_lag_update_ratio_ ) {
            continue;
			}	
		// loop over all particles
		for (size_t i = 0; i < Endparticles.particles.size(); i++) {

			ForestState* thisState = Endparticles.particles[i];
			double weight = thisState->weight();

            // update counts, remove pointers to events that are processed, and remove events when reference count goes to 0
            this->update_star_count( thisState->CoaleventContainer[ epoch_idx ],   weight, x_end, this->total_coal_count[ epoch_idx ],   this->total_weighted_coal_opportunity[ epoch_idx ] );
            this->update_star_count( thisState->RecombeventContainer[ epoch_idx ], weight, x_end, this->total_recomb_count[ epoch_idx ], this->total_weighted_recomb_opportunity[ epoch_idx ] );
            if (this->population_number() > 1){
                this->update_migration_count( thisState->MigreventContainer[ epoch_idx ], weight, x_end, epoch_idx );
                }
            }
        this->counted_to[epoch_idx] = x_end;
        }
    }
