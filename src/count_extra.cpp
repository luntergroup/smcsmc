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

#include"count.hpp"

/*
void CountModel::extract_and_update_count(ParticleContainer &Endparticles, double current_base, bool end_data ) {

    // loop over all particles
	for (size_t i = 0; i < Endparticles.particles.size(); i++) {

		ForestState* thisState = Endparticles.particles[i];
		double weight = thisState->weight();

		// loop over all epochs
		for (size_t epoch_idx = 0; epoch_idx < this->change_times_.size(); epoch_idx++) {
        
			// calculate the required lagging for this epoch; don't use lagging for the final interval
			double lagging = end_data ? 0 : lags[epoch_idx];
			double x_end = current_base - lagging;
			if ( x_end < this->counted_to[epoch_idx] ) continue;
			// Check that we're updating over more than minimal_lag_update_ratio * lagging nucleotides.
			// (If this is the last update, lagging will be 0, and we will do the update)
			// Lagging is used to update the current count with probabilities that are lagging-distanced in the future
			// const_minimal_lag_update_ratio is for optimise purpose, so the count will not updated as frequent
			if ( (x_end - this->counted_to[epoch_idx]) < lagging * this->const_minimal_lag_update_ratio_ ) {
				continue;
			}

			dout << "At epoch "<< epoch_idx<< ", and current_base is at "<< current_base << ", and counting between "<< this->counted_to[epoch_idx] << " to "<< x_end <<endl;
		
            // update counts, remove pointers to events that are processed, and remove events when reference count goes to 0
			this->update_coalescent_count( thisState->eventContainer[ epoch_idx ], weight, x_end, this->total_coal_count[ epoch_idx ],   this->total_weighted_coal_opportunity[ epoch_idx ], epoch_idx );
            if (this->population_number() > 1) {
                this->update_migration_count( thisState->eventContainer[ epoch_idx ], weight, x_end, epoch_idx );
            }
            this->update_recombination_count( thisState->eventContainer[ epoch_idx ], weight, this->counted_to[epoch_idx], x_end, this->total_recomb_count[ epoch_idx ], this->total_weighted_recomb_opportunity[ epoch_idx ], epoch_idx );

            // update counted_to pointers, after we're done with the last particle (and we decided to update this epoch)
            if (i == Endparticles.particles.size() - 1) {
				this->counted_to[epoch_idx] = x_end;
			}
        }
	}
}
*/

void CountModel::extract_and_update_count(ParticleContainer &Endparticles, double current_base, bool end_data ) {

	//
	// calculate update positions per epoch, and identify first epoch to update
	//
	vector<double> update_to;
	size_t first_epoch_to_update = change_times_.size();
	
	for (size_t epoch_idx = 0; epoch_idx < change_times_.size(); epoch_idx++) {

		// calculate the required lagging for this epoch; don't use lagging for the final interval
		double lagging = end_data ? 0 : lags[epoch_idx];
		
		// calculate position to which to update
		double x_end = current_base - lagging;

		// Check that we're updating over more than minimal_lag_update_ratio * lagging nucleotides.
		// (If this is the last update, lagging will be 0, and we will do the update)
		// (If x_end <= counted_to[epoch_idx], the first term will be <= 0 and we will skip this update)
		// Always update if earlier epochs are updating, to ensure that the update boundary is convex.
		// (Note to self: this may cause inferences to be slightly different.  We can in the first instance remove this.)
		if ( (x_end - counted_to[epoch_idx]) < lagging * const_minimal_lag_update_ratio_  && 
		      epoch_idx < first_epoch_to_update ) {
			// no update
			update_to.push_back( counted_to[epoch_idx] );
			//cout << " no update; first_epoch_to_update=" << first_epoch_to_update << endl;
		} else {		
			// update
			update_to.push_back( x_end );
			first_epoch_to_update = min( first_epoch_to_update, epoch_idx );
			//cout << " update to " << x_end << "; first_epoch_to_update=" << first_epoch_to_update << endl;
		}
	}

	first_epoch_to_update = 0;  // update all; for debugging, to make sure the codes are comparable
	
	//
	// update counts for all particles
	//
	for (size_t i = 0; i < Endparticles.particles.size(); i++) {

		ForestState* thisState = Endparticles.particles[i];
		double weight = thisState->weight();
		// loop over the per-epoch containers; this will go in time
		for (size_t epoch_idx = first_epoch_to_update; epoch_idx < change_times_.size(); epoch_idx++) {
			
			//DEBUG -- use original code to evaluate coalescent opp and counts
			double x_end = update_to[ epoch_idx ];
			//double opp1 = update_coalescent_count( thisState->eventContainer[ epoch_idx ], weight, x_end, total_coal_count[ epoch_idx ], total_weighted_coal_opportunity[ epoch_idx ], epoch_idx );
			//update_recombination_count( thisState->eventContainer[ epoch_idx ], weight, counted_to[epoch_idx], x_end, total_recomb_count[ epoch_idx ], total_weighted_recomb_opportunity[ epoch_idx ], epoch_idx );	

			double opp2 = update_all_counts( thisState->eventContainer[ epoch_idx ], weight, update_to, first_epoch_to_update );
			//cout << "Difference " << opp1 - opp2 << endl;
		}
	}
	
	//
    // update counted_to pointers, after we're done with the last particle
    //
    for (size_t epoch_idx = 0; epoch_idx < change_times_.size(); epoch_idx++) {
		counted_to[epoch_idx] = update_to[ epoch_idx ];
	}
}


double CountModel::update_all_counts( deque<EvolutionaryEvent*>& eventContainer, double weight, vector<double>& update_to, size_t first_epoch_to_update ) {

	// Update the relevant counts for all events within the lagging window
	// Consider all epochs; not relevant for current version, but will work later when actual large events are used
    
    //cout << "update_all_counts" << endl;

	double total_opp = 0.0;

	for (int idx=0; idx < eventContainer.size(); idx++) {

		EvolutionaryEvent* event = eventContainer[idx];
		//cout << "  "; event->print_event();

		// consider updating counts for this event.  This is necessary iff the top-left corner is in the update region, since the update region is convex
		if ( event->end_height_epoch() >= first_epoch_to_update && event->start_base() < update_to[ event->end_height_epoch() ] ) {

			//cout << "Considering to update from x=" << event->start_base() << " because this is < update_to = " << update_to[ event->end_height_epoch() ] << endl;

			total_opp += update_all_counts_single_evolevent( event, weight, update_to, first_epoch_to_update );
			
			// if the bottom-right corner has contributed its count, the whole event has, and it can be deleted
			if ( idx == 0 && 
				 event->start_height_epoch() >= first_epoch_to_update && 
				 event->end_base() < update_to[ event->start_height_epoch() ] ) {

				//cout << " -  WILL now remove event " << (long long)event << endl;
				if (event->decrease_refcount_is_zero()) {
					delete event;					
				}
				eventContainer.pop_front();
				idx--;  // account for increment in for loop
			}

		} // evolevent contributes

	} // loop over evolevents
	return total_opp;
}


double CountModel::update_all_counts_single_evolevent( EvolutionaryEvent* event, double weight, vector<double>& update_to, size_t first_epoch_to_update ) {
		
	double total_opp = 0.0;
			
	// loop over all the relevant epochs
	for (size_t epoch_idx = max( first_epoch_to_update, event->start_height_epoch() );
		        epoch_idx <= event->end_height_epoch();
		        epoch_idx++ ) {

		double x_start = counted_to[ epoch_idx ];  // counts have been updated to here
		double x_end = update_to[ epoch_idx ];     // and should be updated to here							

		//cout << "  -- update to " << update_to[epoch_idx] << " (and counted to " << counted_to[epoch_idx] << ")" << endl;		

		// do counts in this epoch require updating?
		if ( event->start_base() < x_end )  {						

			double epoch_start = change_times_[ epoch_idx ];
			double epoch_end = epoch_idx+1 < change_times_.size() ? change_times_[ epoch_idx + 1 ] : DBL_MAX;	

			//cout << "  considerig epoch " << epoch_idx << ":  " << epoch_start << " - " << epoch_end << endl;

			// consider coalescences and migration
			if (event->is_coalmigr()) {

				// any events occur always at the top of the time segment (which cannot be the top of the epoch segment)
				// also ensure that events are not counted twice, i.e. they fall in [x_start,x_end).
				if ( event->end_height < epoch_end && x_start <= event->start_base() ) {
					if (event->is_coal_event()) {
						total_coal_count[ epoch_idx ][event->get_population()].add( weight * event->coal_event_count() );
						//cout << "+CNT = " << weight * event->coal_event_count();
					} else { 
						//cout << "+CNT = 0";
					}
					if (event->is_migr_event()) {   // if !is_migr_event(), event->get_migr_to_population() isn't valid								
						total_mig_count[epoch_idx][event->get_population()][event->get_migr_to_population()].add( weight * event->migr_event_count() );
					}
				} else {
					//cout << "+CNT = 0";
				}
				// account for the opportunity in the segment [epoch_start, epoch_end)
				total_weighted_coal_opportunity[ epoch_idx ][event->get_population()].add( weight * event->coal_opportunity_between( epoch_start, epoch_end ) );
				//cout << " +OPP = " << weight * event->coal_opportunity_between( epoch_start, epoch_end ) << endl;
				//total_opp += weight * event->coal_opportunity_between( epoch_start, epoch_end );
				total_weighted_mig_opportunity[ epoch_idx ][event->get_population()].add( weight * event->migr_opportunity_between( epoch_start, epoch_end ) );
			}											

			// consider recombinations
			if (event->is_recomb()) {

				// the recombination event is considered to be at the top-right corner (arbitrarily -- we don't keep track of recombination locations)
				// also, the recombination event (and opportunity) is arbitrarily assigned to population 0
				if ( x_start <= event->end_base() && event->end_base() < x_end && event->end_height < epoch_end ) {
					total_recomb_count[epoch_idx][0].add( weight * event->recomb_event_count() );
				}
				// account for the opportunity in the appriopriate overlapping rectangle
				total_weighted_recomb_opportunity[epoch_idx][0].add( weight * event->recomb_opportunity_between( epoch_start, epoch_end, x_start, x_end ) );

			}
		} // counts need updating
	} // loop over epochs								
	return total_opp;
}

