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

void CountModel::init() {

    this->init_coal_and_recomb();
    this->init_migr();
    this->init_lags();
        
    this->update_param_interval_  = 5e6; // ONLINE EM
    this->update_param_threshold_ = 1e7; // ONLINE EM

}


void CountModel::init_coal_and_recomb() {

    this->total_coal_count.clear();
    this->total_weighted_coal_opportunity.clear();
    this->total_recomb_count.clear();
    this->total_weighted_recomb_opportunity.clear();

    for (size_t epoch_idx = 0 ; epoch_idx < change_times_.size(); epoch_idx++) {
		// populate coalescent and recombination event counters
        vector <Two_doubles> tmp_count(this->population_number(), Two_doubles(0));
        this->total_coal_count.push_back(tmp_count);
        this->total_recomb_count.push_back(tmp_count);
        // enter initial value
        for (size_t pop_i = 0 ; pop_i < this->population_number(); pop_i++ ) {
            this->total_coal_count[epoch_idx][pop_i] = 1 / ( 2 * population_size() );
        }
        // populate and enter initial value for opportunity
        /*! \todo Need to populate initial count for recombination with correct value */
        vector <Two_doubles> tmp_opportunity(this->population_number(), Two_doubles(1));
        this->total_weighted_coal_opportunity.push_back(tmp_opportunity);
        this->total_weighted_recomb_opportunity.push_back(tmp_opportunity);
    }
}


void CountModel::init_migr() {

    this->total_mig_count.clear();
    this->total_weighted_mig_opportunity.clear();    

	// populate and set up initial values for the event count, opportunity, and inferred rate vectors, for one epoch
	/*! \todo Need to populate with correct initial values */
	vector < vector < Two_doubles > > tmp_count_Time_i;               // Event counts for migrations pop_i -> pop_j
	vector < Two_doubles >            tmp_opp_Time_i;                 // Opportunity  for migrations from pop_i
	vector < vector < double > >      tmp_rate_Time_i_double;         // Rates        for migrations pop_i -> pop_j
	for (size_t pop_i = 0 ; pop_i < this->population_number(); pop_i++ ){
		tmp_count_Time_i.      push_back( vector<Two_doubles>( this->population_number(), Two_doubles( 0.0 ) ) );
		tmp_opp_Time_i.        push_back( Two_doubles(1) );
		tmp_rate_Time_i_double.push_back( vector<double>( this->population_number(), 0 ) );
	}
    
	// set initial counts/rates for all epochs
    for (size_t epoch_idx = 0 ; epoch_idx < change_times_.size(); epoch_idx++) {
        this->total_mig_count.                push_back(tmp_count_Time_i);
        this->total_weighted_mig_opportunity. push_back(tmp_opp_Time_i);
        this->inferred_mig_rate.              push_back(tmp_rate_Time_i_double);
    }
}


void CountModel::init_lags(){

    this->counted_to.clear();
    this->lags.clear();        
    for (size_t epoch_idx = 0 ; epoch_idx < change_times_.size(); epoch_idx++){
        this->counted_to.push_back( (double)0 );
        double top_t = epoch_idx == (change_times_.size() -1) ? change_times_[change_times_.size()-1] : change_times_[epoch_idx+1];
        //double lag_i =  double(4) / this->recombination_rate() / top_t ; 
        double lag_i = this->const_lag_ > 0 ? this->const_lag_ : double(4) / (this->recombination_rate() * top_t) ; 
        this->lags.push_back( lag_i );
    }    
}


void CountModel::initialize_mig_rate ( vector <vector<double>*> & rates_list ){
    for (size_t i = 0; i < rates_list.size(); i++ ){
        if (rates_list[i]){
            for (size_t j = 0; j < rates_list[i]->size() ; j++){
                rates_list[i]->at(j) = 0;
            }
        }    
    }
}
    
    
void CountModel::reset_Ne ( Model *model ){

    for (size_t epoch_idx = 0; epoch_idx < change_times_.size(); epoch_idx++){
        for (size_t pop_j = 0 ; pop_j < this->population_number(); pop_j++ ){
            //model->addPopulationSize(this->change_times_[epoch_idx], pop_j, this->total_weighted_coal_opportunity[epoch_idx][pop_j] / this->total_coal_count[epoch_idx][pop_j] /2 ,false, false);    
            this->total_weighted_coal_opportunity[epoch_idx][pop_j].compute_final_answer();
            this->total_coal_count[epoch_idx][pop_j].compute_final_answer();
            double tmp_pop_size = this->total_weighted_coal_opportunity[epoch_idx][pop_j].final_answer() / this->total_coal_count[epoch_idx][pop_j].final_answer() / (double)2;
            //tmp_pop_size = roundf(tmp_pop_size * (double)1e6)/(double)1e6; // Rounding, Gerton: this is not a good idea
            model->addPopulationSize(this->change_times_[epoch_idx], pop_j, tmp_pop_size ,false, false);    
            cout << " popsize is equal to " << tmp_pop_size << " ( "<<this->total_weighted_coal_opportunity[epoch_idx][pop_j].final_answer()<<" / "<< this->total_coal_count[epoch_idx][pop_j].final_answer() << "/2)" <<endl;
        }
    }
    this->check_model_updated_Ne( model );
}


void CountModel::reset_recomb_rate ( Model *model ){
    #ifdef _RecombRecordingOff // if the recombination recording off macro is defined, then return without recording the event
        return;
    #endif
    this->compute_recomb_rate();
    //this->inferred_recomb_rate = roundf ( this->inferred_recomb_rate * (double)1e16)/(double)1e16; // Rounding, Gerton: this is not a good idea
    model->setRecombinationRate( this->inferred_recomb_rate , false, false, 0);
    cout << " set recombination rate " << model->recombination_rate(0) << "("<<this->recomb_count_<<" / "<< this->recomb_opportunity_ << ")" <<endl;
}


void CountModel::reset_mig_rate ( Model *model ) {

    if (!this->has_migration()) 
		return;
		
    this->compute_mig_rate();
    
    assert( this->print_mig_rate (model->mig_rates_list_) );
    assert( this->print_mig_rate (model->total_mig_rates_list_) );

    this->initialize_mig_rate ( model->mig_rates_list_ );
    this->initialize_mig_rate ( model->total_mig_rates_list_ );

    assert( this->print_mig_rate (model->mig_rates_list_) );
    assert( this->print_mig_rate (model->total_mig_rates_list_) );

    for (size_t epoch_idx = 0; epoch_idx < change_times_.size(); epoch_idx++){
        for (size_t pop_j = 0 ; pop_j < this->population_number(); pop_j++ ){
            for (size_t pop_k = 0 ; pop_k < this->population_number(); pop_k++ ) {
                if ( pop_j == pop_k) continue;
                model->addMigrationRate(this->change_times_[epoch_idx], pop_j, pop_k, this->inferred_mig_rate[epoch_idx][pop_j][pop_k], false, false);
            }
        }
    }

    this->check_model_updated_mig (model);
    
    assert( this->print_mig_rate (model->mig_rates_list_) );
    assert( this->print_mig_rate (model->total_mig_rates_list_) );
}


void CountModel::reset_model_parameters(double current_base, Model * model, bool online, bool force_update, bool print){
    
    bool update = false;
    if ( online && current_base > this->update_param_threshold_ ){
        update = true;
        this->update_param_threshold_ += this->update_param_interval_;
        }
                    
    if ( update || force_update  ){ 
        cout<<" MODEL IS RESET at base " << current_base <<endl;
        this->reset_recomb_rate ( model );
        this->reset_Ne ( model );
        this->reset_mig_rate ( model );
        //this->reset_single_mig_rate ( model );
        model->finalize();
        } 
    else {
        if (print){
            this->print_coal_count();            
            }
        return;
        }
    }


void CountModel::log_counts( PfParam& param ) {
	// log coalescent counts
	for (size_t epoch_idx = 0; epoch_idx < change_times_.size(); epoch_idx++ ) {
		for (size_t pop_idx = 0; pop_idx < this->population_number(); pop_idx++ ) {
			param.append_to_count_file( epoch_idx, "Coal", pop_idx, -1, this->total_weighted_coal_opportunity[epoch_idx][pop_idx].final_answer(), 
																	    this->total_coal_count[epoch_idx][pop_idx].final_answer() );
		}
	}
	// log recombination counts
	for (size_t epoch_idx = 0; epoch_idx < change_times_.size(); epoch_idx++ ) {
		param.append_to_count_file( epoch_idx, "Recomb", -1, -1, this->total_weighted_recomb_opportunity[epoch_idx][0].final_answer(), 
																 this->total_recomb_count[epoch_idx][0].final_answer() );
	}
	// log migration counts
	for (size_t epoch_idx = 0; epoch_idx < change_times_.size(); epoch_idx++ ) {
		for (size_t from_pop_idx = 0; from_pop_idx < this->population_number(); from_pop_idx++ ) {
			for (size_t to_pop_idx = 0; to_pop_idx < this->population_number(); to_pop_idx++ ) {
				if (from_pop_idx != to_pop_idx) {
					param.append_to_count_file( epoch_idx, "Migr", from_pop_idx, to_pop_idx, this->total_weighted_mig_opportunity[epoch_idx][from_pop_idx].final_answer(),
																						     this->total_mig_count[epoch_idx][from_pop_idx][to_pop_idx].final_answer() );
				}
			}
		}
	}
}

    

void CountModel::compute_mig_rate(){

    for (size_t epoch_idx = 0; epoch_idx<change_times_.size(); epoch_idx++){
        for (size_t pop_i = 0 ; pop_i < this->population_number(); pop_i++ ){
			
			this->total_weighted_mig_opportunity[epoch_idx][pop_i].compute_final_answer();
            
            for (size_t pop_j = 0 ; pop_j < this->population_number(); pop_j++ ){
				this->total_mig_count[epoch_idx][pop_i][pop_j].compute_final_answer();
                this->inferred_mig_rate                 [epoch_idx][pop_i][pop_j] = 
					this->total_mig_count               [epoch_idx][pop_i][pop_j].final_answer() / 
					this->total_weighted_mig_opportunity[epoch_idx][pop_i].final_answer();
                //cout<<"this->inferred_mig_rate["<<epoch_idx<<"]["<<pop_i<<"]["<<pop_j<<"] = " << this->inferred_mig_rate[epoch_idx][pop_i][pop_j]<<endl;
            }
        }
    }
}


/*! 
 * \brief Compute the recombination rate once we have sweeped through all the data, and recorded the recomb_opportunity and recomb_counts
 */ 
void CountModel::compute_recomb_rate () {

    this->recomb_opportunity_ = 0;
    this->recomb_count_ = 0;
    for ( size_t epoch_idx = 0; epoch_idx < change_times_.size(); epoch_idx++ ){
		this->total_weighted_recomb_opportunity[epoch_idx][0].compute_final_answer();
		this->total_recomb_count[epoch_idx][0].compute_final_answer();
        this->recomb_opportunity_ += this->total_weighted_recomb_opportunity[epoch_idx][0].final_answer() ;
        this->recomb_count_ += this->total_recomb_count[epoch_idx][0].final_answer();
    }
    
    this->inferred_recomb_rate = this->recomb_count_ / this->recomb_opportunity_;
    //cout << "Recombination rate is " << this->inferred_recomb_rate << " ( " << this->recomb_count_ << " / " << this->recomb_opportunity_ << " )"<<endl;
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
                    

double CountModel::update_coalescent_count( deque<EvolutionaryEvent*>& eventContainer_i, double weight, double x_end, vector<Two_doubles>& total_coal_count, vector<Two_doubles>& total_coal_opportunity, size_t epoch_idx ){
    // Go through the events, starting from the leftmost and going up to x_end, and add events (weighted by weight) to the appropriate counters
    // When processed remove the event pointer from the deque; remove the event itself if its reference count becomes 0
    int idx = 0;
    double total_opp = 0.0;
    //cout << "update_coalescent_count -- update to " << x_end << " (and counted to " << counted_to[epoch_idx] << ")" << endl;
	for (idx=0; idx < eventContainer_i.size(); ) {
		
		EvolutionaryEvent* event = eventContainer_i[idx];
		if (event->is_coalmigr()) {
			if (event->end_base() >= x_end) break;
			//cout << "  "; event->print_event();
			//cout << "-CNT = " << weight * event->coal_event_count() << " -OPP = " << weight * event->coal_opportunity() << endl;
			total_opp += weight * event->coal_opportunity();
			total_coal_opportunity[event->get_population()].add( weight * event->coal_opportunity() );
			total_coal_count[event->get_population()].add( weight * event->coal_event_count() );
		}
		++idx;
	}
	return total_opp;
}

void CountModel::update_migration_count( deque<EvolutionaryEvent*>& eventContainer_i, double weight, double x_end, size_t epoch_idx ) {
    // Go through the events, starting from the leftmost and going up to x_end, and add events (weighted by weight) to the appropriate counters
    // When processed remove the event pointer from the deque; remove the event itself if its reference count becomes 0
	for (int idx=0; idx < eventContainer_i.size(); ) {
		
		EvolutionaryEvent* event = eventContainer_i[idx];
		if (event->is_coalmigr()) {
			if (event->end_base() >= x_end) break;
			total_weighted_mig_opportunity[epoch_idx][event->get_population()].add( weight * event->migr_opportunity() );
			if (event->is_migr_event()) {
				// This is not just an efficiency measure: if !is_migr_event(), event->get_migr_to_population() isn't valid
				total_mig_count[epoch_idx][event->get_population()][event->get_migr_to_population()].add( weight * event->migr_event_count() );
			}
		}
		++idx;

	}    
}



void CountModel::update_recombination_count( deque<EvolutionaryEvent*>& eventContainer_i, double weight, double x_start, double x_end, vector<Two_doubles>& total_recomb_count, vector<Two_doubles>& total_recomb_opportunity, size_t epoch_idx ){
    #ifdef _RecombRecordingOff // if the recombination recording off macro is defined, then return without recording the event
        return;
    #endif
         
    // Go through the events, starting from the leftmost and going up to x_end, and add events (weighted by weight) to the appropriate counters
    // When processed remove the event pointer from the deque; remove the event itself if its reference count becomes 0

    // First process all recombination event segments that lie fully to the left of x_end
    int idx;
	for (idx=0; idx < eventContainer_i.size(); ) {
		
		EvolutionaryEvent* event = eventContainer_i[idx];
		if (event->is_recomb()) {
			if (event->start_base() >= x_end) break;
			if ( event->start_base() >= x_start && event->start_base() < x_end ) {
				total_recomb_count[0].add( weight * event->recomb_event_count() );
			}
			double t_start = event->start_height;
			double t_end = event->end_height;
			total_recomb_opportunity[0].add( weight * event->recomb_opportunity_between( t_start, t_end, x_start, x_end ) );
		}

		// now, all types are processed, so remove events (of any type) that are completely done.
        if (idx == 0 && event->end_base() < x_end) {
			if (event->decrease_refcount_is_zero()) {
				delete event;
			}
			eventContainer_i.pop_front();
			//cout << " - WOULD now remove event " << (long long)event << " at (consecutive..?!) index " << idx << endl; ++idx;
		} else {
			++idx;
		}
	}
}


