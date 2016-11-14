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

#include "count.hpp"
#include "descendants.hpp"


void CountModel::init() {

    this->init_coal_and_recomb();
    this->init_migr();
    this->init_local_recomb();
    this->init_lags();

    this->update_param_interval_  = 5e6; // ONLINE EM
    this->update_param_threshold_ = 1e7; // ONLINE EM

}


void CountModel::reset_model_parameters(double current_base, Model * model, bool useCap, double cap, bool online, bool force_update, bool print){

    bool update = false;
    if ( online && current_base > this->update_param_threshold_ ) {
        update = true;
        this->update_param_threshold_ += this->update_param_interval_;
    }

    if ( update || force_update  ) {
        clog<<" MODEL IS RESET at base " << current_base <<endl;
        this->reset_recomb_rate ( model );
        this->reset_Ne ( model, useCap, cap );
        this->reset_mig_rate ( model );
        model->finalize();
    } else {
        if (print) {
            this->print_coal_count();
        }
    }
}


void CountModel::log_counts( PfParam& param ) {
    // log coalescent counts
    for (size_t epoch_idx = 0; epoch_idx < change_times_.size(); epoch_idx++ ) {
        for (size_t pop_idx = 0; pop_idx < this->population_number(); pop_idx++ ) {
            param.appendToOutFile( param.EMcounter(),
                                   (int)epoch_idx,
                                   change_times_[epoch_idx],
                                   ( epoch_idx == change_times_.size()-1 ) ? 1e+99 : change_times_[epoch_idx+1],
                                   "Coal",
                                   pop_idx,
                                   -1,
                                   this->total_coal_opportunity[epoch_idx][pop_idx].final_answer(),
                                   this->total_coal_count[epoch_idx][pop_idx].final_answer(),
                                   this->total_coal_weight[epoch_idx][pop_idx].final_answer());
        }
    }

    // log recombination counts
    double recomb_opportunity = 0;
    double recomb_count = 0;
    double recomb_weight = 0;
    for (size_t epoch_idx = 0; epoch_idx < change_times_.size(); epoch_idx++ ) {
        recomb_opportunity += total_recomb_opportunity[epoch_idx][0].final_answer();
        recomb_count       += total_recomb_count[epoch_idx][0].final_answer();
        recomb_weight      += total_recomb_weight[epoch_idx][0].final_answer();
        if ( false ){ // when false, do not print out recomb events at every epoch.
            param.appendToOutFile( param.EMcounter(),
                                   (int)epoch_idx,
                                   change_times_[epoch_idx],
                                   ( epoch_idx == change_times_.size()-1 ) ? 1e+99 :change_times_[epoch_idx+1],
                                   "Recomb",
                                   -1,
                                   -1,
                                   this->total_recomb_opportunity[epoch_idx][0].final_answer(),
                                   this->total_recomb_count[epoch_idx][0].final_answer(),
                                   this->total_recomb_weight[epoch_idx][0].final_answer());
        }
    }
    param.appendToOutFile( param.EMcounter(),
                           -1,
                           0.0,
                           1e+99,
                           "Recomb",
                           -1,
                           -1,
                           recomb_opportunity,
                           recomb_count,
                           recomb_weight);

    // log migration counts
    for (size_t epoch_idx = 0; epoch_idx < change_times_.size(); epoch_idx++ ) {
        for (size_t from_pop_idx = 0; from_pop_idx < this->population_number(); from_pop_idx++ ) {
            for (size_t to_pop_idx = 0; to_pop_idx < this->population_number(); to_pop_idx++ ) {
                if (from_pop_idx != to_pop_idx) {
                    param.appendToOutFile( param.EMcounter(),
                                           (int)epoch_idx,
                                           change_times_[epoch_idx],
                                           ( epoch_idx == change_times_.size()-1 ) ? 1e+99 : change_times_[epoch_idx+1],
                                           "Migr",
                                           from_pop_idx,
                                           to_pop_idx,
                                           this->total_mig_opportunity[epoch_idx][from_pop_idx].final_answer(),
                                           this->total_mig_count[epoch_idx][from_pop_idx][to_pop_idx].final_answer(),
                                           this->total_mig_weight[epoch_idx][from_pop_idx].final_answer());
                }
            }
        }
    }
}


void CountModel::init_coal_and_recomb() {

    total_coal_count.clear();
    total_coal_opportunity.clear();
    total_coal_weight.clear();
    total_recomb_count.clear();
    total_recomb_opportunity.clear();
    total_recomb_weight.clear();

    resetTime();
    for (size_t epoch_idx = 0 ; epoch_idx < change_times_.size(); epoch_idx++) {
        // move to the next epoch
        if (epoch_idx > 0)
            increaseTime();
        // populate coalescent and recombination event counters
        vector <TwoDoubles> tmp_count(this->population_number(), TwoDoubles(0));
        total_coal_count.push_back(tmp_count);
        total_recomb_count.push_back(tmp_count);
        // enter initial value
        for (size_t pop_i = 0 ; pop_i < this->population_number(); pop_i++ ) {
            total_coal_count[epoch_idx][pop_i] = 1 / ( 2 * this->population_size( pop_i ) );
            /*! Note that this uses recombination rate at position -1 */
            total_recomb_count[epoch_idx][pop_i] = this->recombination_rate();
        }
        // populate and enter initial value for opportunity
        vector <TwoDoubles> tmp_opportunity(this->population_number(), TwoDoubles(1));
        total_coal_opportunity.push_back(tmp_opportunity);
        total_recomb_opportunity.push_back(tmp_opportunity);
        // and same for weight counters
        total_coal_weight.push_back(tmp_opportunity);
        total_recomb_weight.push_back(tmp_opportunity);
    }
}


void CountModel::init_migr() {

    total_mig_count.clear();
    total_mig_opportunity.clear();
    total_mig_weight.clear();

    // set initial counts/rates for all epochs
    resetTime();
    for (size_t epoch_idx = 0 ; epoch_idx < change_times_.size(); epoch_idx++) {

        // move to the next epoch
        if (epoch_idx > 0)
            increaseTime();

        // populate and set up initial values for the event count, opportunity, and inferred rate vectors, for one epoch
        vector < vector < TwoDoubles > > tmp_count_Time_i;               // Event counts for migrations pop_i -> pop_j
        vector < TwoDoubles >            tmp_opp_Time_i;                 // Opportunity  for migrations from pop_i
        vector < vector < double > >      tmp_rate_Time_i_double;         // Rates        for migrations pop_i -> pop_j
        for (size_t pop_i = 0 ; pop_i < this->population_number(); pop_i++ ){
            tmp_count_Time_i.      push_back( vector<TwoDoubles>( this->population_number(), TwoDoubles( 0.0 ) ) );
            tmp_opp_Time_i.        push_back( TwoDoubles(1) );
            tmp_rate_Time_i_double.push_back( vector<double>( this->population_number(), 0 ) );
            for (size_t pop_j = 0 ; pop_j < this->population_number(); pop_j++) {
                tmp_count_Time_i[ pop_i ][ pop_j ]       = this->migration_rate( pop_i, pop_j );
                tmp_rate_Time_i_double[ pop_i ][ pop_j ] = this->migration_rate( pop_i, pop_j );
            }
        }
        total_mig_count.      push_back(tmp_count_Time_i);
        total_mig_opportunity.push_back(tmp_opp_Time_i);
        total_mig_weight.     push_back(tmp_opp_Time_i);
    }
}


void CountModel::init_lags(){

    this->counted_to.clear();
    this->lags.clear();
    for (size_t epoch_idx = 0 ; epoch_idx < change_times_.size(); epoch_idx++){
        this->counted_to.push_back( (double)0 );
        double top_t = epoch_idx == (change_times_.size() -1) ? change_times_[change_times_.size()-1] : change_times_[epoch_idx+1];
        double lag_i = this->const_lag_ > 0 ? this->const_lag_ : double(4) / (this->recombination_rate() * top_t) ;
        this->lags.push_back( lag_i );
    }
}

void CountModel::init_local_recomb() {

    local_recomb_opportunity.clear();
    local_recomb_counts.clear();
    for (int sample = 0; sample < sample_size_; ++sample) {
        local_recomb_counts.push_back( vector<double>() );
    }

}


// Do we want to allow for different lags in diff populations?
void CountModel::reset_lag ( std::vector<double> survival, double lag_fraction ){
    for (size_t epoch_idx = 0; epoch_idx < change_times_.size(); epoch_idx++) {
        lags.at(epoch_idx) = survival.at(epoch_idx) * lag_fraction;
    }
}

void CountModel::reset_Ne ( Model *model, bool useCap, double cap){

    for (size_t epoch_idx = 0; epoch_idx < change_times_.size(); epoch_idx++) {
        for (size_t pop_j = 0 ; pop_j < this->population_number(); pop_j++ ) {
            double coal_opp    = total_coal_opportunity[epoch_idx][pop_j].final_answer();
            double coal_count  = total_coal_count[epoch_idx][pop_j].final_answer();
            double coal_weight = total_coal_weight[epoch_idx][pop_j].final_answer();
            double coal_rate  = coal_count / coal_opp;
            double pop_size   = 1.0 / (2.0 * coal_rate);
            if ( useCap ){
				if ( pop_size >= cap ){
					model->addPopulationSize(change_times_[epoch_idx], pop_j, cap ,false, false);
				} else {
					model->addPopulationSize(change_times_[epoch_idx], pop_j, pop_size ,false, false);
				}
			} else {
				model->addPopulationSize(change_times_[epoch_idx], pop_j, pop_size ,false, false);
			}

            clog << " Setting size of population " << pop_j << " @ " << setw(8) << change_times_[epoch_idx] << " to "
                 << setw(8) << pop_size
                 << " ( 0.5 * " << total_coal_opportunity[epoch_idx][pop_j].final_answer() << " / " << this->total_coal_count[epoch_idx][pop_j].final_answer() << "; post-lag ESS "
                 << 1.0 / (coal_weight / coal_opp) << " )" << endl;
        }
    }
    this->check_model_updated_Ne( model );
}


void CountModel::reset_recomb_rate ( Model *model ){

    double recomb_opportunity = 0;
    double recomb_count = 0;
    double recomb_weight = 0;

    for ( size_t epoch_idx = 0; epoch_idx < change_times_.size(); epoch_idx++ ) {
        recomb_opportunity += total_recomb_opportunity[epoch_idx][0].final_answer();
        recomb_count       += total_recomb_count[epoch_idx][0].final_answer();
        recomb_weight      += total_recomb_weight[epoch_idx][0].final_answer();
    }
    double inferred_recomb_rate = recomb_count / recomb_opportunity;

    model->setRecombinationRate( inferred_recomb_rate, false, false, 0);
    clog << " Setting recombination rate to " << model->recombination_rate(0)
         << " ( " << recomb_count<<" / " << recomb_opportunity << "; post-lag ESS "
         << 1.0 / (recomb_weight / recomb_opportunity) << " )" << endl;
}


void CountModel::clear_2d_vector ( vector <vector<double>> & rates_list ) {
    for (size_t i = 0; i < rates_list.size(); i++ ) {
        if (!rates_list[i].empty()) {
            for (size_t j = 0; j < rates_list[i].size() ; j++) {
                rates_list[i][j] = 0;
            }
        }
    }
}


void CountModel::reset_mig_rate ( Model *model ) {

    if (!has_migration()) return;

    clear_2d_vector( model->mig_rates_list_ );
    clear_2d_vector( model->total_mig_rates_list_ );

    vector < vector < vector<double> > > inferred_mig_rate;
    for (size_t epoch_idx = 0; epoch_idx<change_times_.size(); epoch_idx++) {
        for (size_t pop_i = 0 ; pop_i < this->population_number(); pop_i++ ) {
            for (size_t pop_j = 0 ; pop_j < this->population_number(); pop_j++ ) {
                if ( pop_i == pop_j) continue;
                double migration_rate =
                    total_mig_count[epoch_idx][pop_i][pop_j].final_answer() /
                    total_mig_opportunity[epoch_idx][pop_i].final_answer();

                model->addMigrationRate(change_times_[epoch_idx], pop_i, pop_j, migration_rate, false, false);
            }
        }
    }

    this->check_model_updated_mig (model);

    assert( print_mig_rate (model->mig_rates_list_) );
    assert( print_mig_rate (model->total_mig_rates_list_) );
}


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
        // Always update if earlier epochs are updating, to ensure that the update boundary is nondecreasing.
        // NOTE: now that we're goign back to per-epoch storign of events, this is no longer necessary.
        if ( (x_end - counted_to[epoch_idx]) < lagging * const_minimal_lag_update_ratio_  &&
             (first_epoch_to_update > epoch_idx ) ) {
            // no update
            update_to.push_back( counted_to[epoch_idx] );
        } else {
            // update
            update_to.push_back( x_end );
            first_epoch_to_update = min( first_epoch_to_update, epoch_idx );
        }
    }

    //
    // update counts for all particles
    //
    for (int i = Endparticles.particles.size()-1; i>=0; i--) {

        ForestState* thisState = Endparticles.particles[i];

        for (int epoch_idx = change_times_.size()-1; epoch_idx >= (int)first_epoch_to_update; epoch_idx--) {

            update_all_counts( &thisState->eventTrees[ epoch_idx ], thisState->weight(), update_to, epoch_idx );

        }
    }

    //
    // update counted_to pointers, after we're done with the last particle
    //
    for (size_t epoch_idx = 0; epoch_idx < change_times_.size(); epoch_idx++) {
        counted_to[epoch_idx] = update_to[ epoch_idx ];
    }
}



                    /*! \verbatim
                            xstart
                            .                      xend                         VCFfile->site()
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
                     \endverbatim
                     *
                     * Count the coalescent events between position xstart and xend.
                     *
                     * At the beginning of this function, the tail ForestState is at
                     * state 6, whose weight represents the weight for the entire particle
                     * As lagging is applied, we need to skip a few states before start counting.
                     *
                     * In this example, only count the coalescent events occured on states 1 and 2.
                     */



void CountModel::update_all_counts( EvolutionaryEvent** event_ptr, double posterior_weight, vector<double>& update_to, size_t epoch_idx ) {

    // Recursively traverse the tree, updating complete nodes, until an incomplete node or the root is found

    EvolutionaryEvent* event;
    // Find first non-deleted event, updating *event_ptr (but not event_ptr) as necessary; exit if root was found
    while ((event = purge_events( event_ptr, epoch_idx ))) {

        if (!event->update_posterior_is_done( posterior_weight )) {
            // incomplete node encountered -- stop, and wait until other updates complete this node
            break;
        }

        // destructively read accumulated posterior weight
        posterior_weight = event->get_and_reset_posterior();

        // update counters if top-left corner is in update region
        if ( event->start_base() < update_to[ epoch_idx ] ) {

            update_all_counts_single_evolevent( event, posterior_weight, update_to, epoch_idx );

            // if the bottom-right corner has contributed its count, the whole event has, and it can be deleted
            if ( event->end_base() < update_to[ epoch_idx ] ) {

                event->mark_as_removed();
                if (!remove_event( event_ptr, epoch_idx )) {
                    // event was not removed, but a new pointer now points to the parent, from an event
                    // that has already been updated and therefore will not update the parent.
                    // Do an empty update on the parent to account for this.  (The alternative is to mark event
                    //  as 'removed' and purge it on the next iteration, which is more straightforward and easier
                    //  to understand; however the event is in the cache so might as well make use of that.)
                    event = *event_ptr;
                    if (event && !event->is_removed()) {
                        event->update_posterior_is_done( 0.0 );
                    }
                }
                continue;
            }
        }

        // move event_ptr to point to (pointer to parent of) this event.
        event_ptr = &event->parent();

    }
}



void CountModel::update_all_counts_single_evolevent( EvolutionaryEvent* event, double weight, vector<double>& update_to, size_t epoch_idx ) {

    double x_start = counted_to[ epoch_idx ];  // counts have been updated to here
    double x_end = update_to[ epoch_idx ];     // and should be updated to here

    // do counts in this epoch require updating?
    if ( event->start_base() < x_end )  {

        double epoch_start = change_times_[ epoch_idx ];
        double epoch_end = epoch_idx+1 < change_times_.size() ? change_times_[ epoch_idx + 1 ] : DBL_MAX;

        // consider coalescences and migration
        if (event->is_coalmigr()) {

            // any events occur always at the top of the time segment (which cannot be the top of the epoch segment)
            // also ensure that events are not counted twice, i.e. they fall in [x_start,x_end).
            int pop = event->get_population();
            if ( x_start <= event->start_base() ) {
                if (event->is_coal_event()) {   // not stricly necessary, as if !is_coal_event, coal_event_count()==0
                    total_coal_count[ epoch_idx ][ pop ].add( weight * event->coal_event_count() );
                }
                if (event->is_migr_event()) {   // if !is_migr_event(), event->get_migr_to_population() is not valid
                    total_mig_count[epoch_idx][ pop ][ event->get_migr_to_population() ].add( weight * event->migr_event_count() );
                }
            }
            // account for the opportunity in the segment [epoch_start, epoch_end)
            double coal_opp = event->coal_opportunity_between( epoch_start, epoch_end );
            double migr_opp = event->migr_opportunity_between( epoch_start, epoch_end );
            total_coal_opportunity[ epoch_idx ][ pop ].add( weight * coal_opp );
            total_coal_weight[ epoch_idx ][ pop ].add( weight * weight * coal_opp );
            total_mig_opportunity[ epoch_idx ][ pop ].add( weight * migr_opp );
            total_mig_weight[ epoch_idx ][ pop ].add( weight * weight * migr_opp );
        }

        // consider recombinations
        if (event->is_recomb()) {
            bool isEndOfSeq = this->loci_length() == x_end;
            // the recombination event (and opportunity) is arbitrarily assigned to population 0
            double recomb_opp = event->recomb_opportunity_between( epoch_start, epoch_end, x_start, x_end );
            total_recomb_count[epoch_idx][ 0 ].add( weight * event->recomb_event_count_between( epoch_start, epoch_end, x_start, x_end, isEndOfSeq ) );
            total_recomb_opportunity[epoch_idx][ 0 ].add( weight * recomb_opp );
            total_recomb_weight[ epoch_idx ][ 0 ].add( weight * weight * recomb_opp );
            double local_x_start = max( x_start, event->get_start_base() );
            double local_x_end = min( x_end, event->get_end_base() );
            record_local_recomb_events( local_x_start, local_x_end, weight, recomb_opp,
                                        event->recomb_event_base(), event->get_descendants() );
        }
    }
}



void CountModel::record_local_recomb_events( double x_start, double x_end, double weight, double opportunity, double event_base,
                                             Descendants_t descendants ) {

    size_t first_index = size_t( x_start / local_recording_interval_ );
    size_t last_index = size_t( 1 + x_end / local_recording_interval_ );
    double first_interval = min( (first_index+1)*local_recording_interval_, x_end ) - x_start;
    double last_interval = x_end - max( (last_index-1)*local_recording_interval_, x_start );
    double opp_density = weight * opportunity / (x_end - x_start);

    // ensure the count vectors are large enough; if not add space for 10Mb worth of samples
    if (local_recomb_opportunity.size() <= last_index) {
        size_t new_size = local_recomb_opportunity.size() + (size_t)(1e7 / local_recording_interval_);
        local_recomb_opportunity.resize( new_size );
        for (int sample = 0; sample < sample_size_; ++sample) {
            local_recomb_counts[ sample ].resize( new_size );
        }
    }

    // store differential opportunity
    if (first_index == last_index - 1) {
        local_recomb_opportunity[ first_index ]   += first_interval * opp_density;
        local_recomb_opportunity[ first_index+1 ] -= first_interval * opp_density;
    } else {
        local_recomb_opportunity[ first_index ]   += first_interval * opp_density;
        local_recomb_opportunity[ first_index+1 ] += (local_recording_interval_ - first_interval) * opp_density;
        local_recomb_opportunity[ last_index-1 ]  += (last_interval - local_recording_interval_) * opp_density;
        local_recomb_opportunity[ last_index ]    -= last_interval * opp_density;
    }

    // record recombination event, if any, if it falls within the segment
    if (x_start <= event_base && event_base <= x_end) {

        // count the number of descendants, for normalisation
        int descendant = 0, num_descendants = 0;
        Descendants_t data = descendants;
        while (get_next_descendant( data, descendant )) {
            ++num_descendants;
        }
        
        // store the count.  Note: the 'descendant' values are 1-based!
        size_t index = size_t( event_base / local_recording_interval_ );
        descendant = 0;
        data = descendants;
        while (get_next_descendant( data, descendant )) {
            assert (descendant > 0);
            assert (descendant < sample_size_);
            local_recomb_counts[ descendant-1 ][ index ] += weight / num_descendants;
        }
    }
};


void CountModel::dump_local_recomb_logs( ostream& stream, double locus_length, int iteration ) {

    // find last nonzero entry
    size_t last_idx = min( local_recomb_opportunity.size(), (size_t)(locus_length / local_recording_interval_) );

    // write header
    if (iteration == 0) {
        stream << "iter\tlocus\tsize\topp_per_nt";
        for (int sample = 0; sample < sample_size_; ++sample) {
            stream << "\t" << sample+1;
        }
        stream << "\n";
    }

    // write data; convert local recomb opportunity CHANGES into absolute
    // recomb opportunity on-the-fly.  (Todo: change name of local_rec_opp.)
    size_t idx = 0;
    double current_opportunity = 0.0;
    // invariant: current_opportunity = sum (0 <= i < idx) local_recomb_opportunity[i]
    while (idx < last_idx) {
        int num_indices = 0, first_idx = idx;
        double total_count = 0;
        do {
            current_opportunity += local_recomb_opportunity[ idx ];
            for (int sample = 0; sample < sample_size_; ++sample) {
                total_count += local_recomb_counts[ sample ][ idx ];
            }
            ++num_indices;
            ++idx;
        } while (total_count == 0.0 &&
                 idx < last_idx &&
                 abs(local_recomb_opportunity[ idx ]) < 1e-4 * current_opportunity);
        // summarize streak of null counts
        int streak = num_indices - 1;
        if (total_count == 0.0) {
            // loop was ended because of changing opportunity, or end-of-data,
            // rather than a nonzero count, so the zero-count streak includes idx-1
            ++streak;
        }
        if (streak > 0) {
            stream << iteration
                   << "\t"
                   << fixed << setprecision(0)
                   << first_idx * local_recording_interval_
                   << "\t"
                   << streak * local_recording_interval_
                   << "\t"
                   << scientific << setprecision(5)
                   << current_opportunity / local_recording_interval_;
            for (int sample = 0; sample < sample_size_; ++sample) {
                stream << "\t" << 0.0;
            }
            stream << "\n";
            // point to first index that needs processing
            while (idx < first_idx + streak) {
                current_opportunity += local_recomb_opportunity[ idx ];
                ++idx;
            }
        } else {
            // streak == 0 && total_count > 0 && idx == first_idx + 1
            // (i.e. we found a single nonzero-count line)
            stream << iteration
                   << "\t" 
                   << fixed << setprecision(0)
                   << first_idx * local_recording_interval_
                   << "\t"
                   << local_recording_interval_
                   << "\t"
                   << scientific << setprecision(5)
                   << current_opportunity / local_recording_interval_;
            for (int sample = 0; sample < sample_size_; ++sample) {
                stream << "\t" << local_recomb_counts[ sample ][ first_idx ] / local_recording_interval_;
            }
            stream << "\n";
        } 
    }
}
