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

void CountModel::init(){
    this->init_coal_and_recomb();
    this->init_migr();
    this->init_lags();
        
    this->update_param_interval_  = 5e6; // ONLINE EM
    this->update_param_threshold_ = 1e7; // ONLINE EM
    return ;
    }


void CountModel::init_coal_and_recomb(){
    this->total_coal_count.clear();
    this->total_weighted_coal_opportunity.clear();
    this->total_recomb_count.clear();
    this->total_weighted_recomb_opportunity.clear();

    this->resetTime();    
    for (size_t epoch_idx = 0 ; epoch_idx < change_times_.size(); epoch_idx++){
        vector <double> tmp_count(this->population_number(), 0);
        this->total_coal_count.push_back(tmp_count);
        this->total_recomb_count.push_back(tmp_count);
        for (size_t pop_i = 0 ; pop_i < this->population_number(); pop_i++ ){
            this->total_coal_count[epoch_idx][pop_i] = 1 / ( 2 * population_size() );
        }
        
        vector <double> tmp_opportunity(this->population_number(), 1);
        this->total_weighted_coal_opportunity.push_back(tmp_opportunity);
        this->total_weighted_recomb_opportunity.push_back(tmp_opportunity);
        }
    this->resetTime();        
    }


void CountModel::init_migr(){ /*! \todo This requires more work*/
    this->total_mig_count.clear();
    this->total_weighted_mig_opportunity.clear();    

    this->resetTime();    
    for (size_t epoch_idx = 0 ; epoch_idx < change_times_.size(); epoch_idx++){
        vector < vector < double > > tmp_count_Time_i;
        vector < double > tmp_opp_Time_i;
        for (size_t pop_i = 0 ; pop_i < this->population_number(); pop_i++ ){
            vector <double> tmp_count(this->population_number(), 0);
            tmp_count_Time_i.push_back(tmp_count);
            tmp_opp_Time_i.push_back( 1 );
            }        
        this->total_mig_count.push_back(tmp_count_Time_i);
        this->inferred_mig_rate.push_back(tmp_count_Time_i);
        this->total_weighted_mig_opportunity.push_back(tmp_opp_Time_i);
        }
    }


void CountModel::init_lags(){
    this->previous_base.clear();
    this->lags.clear();        
    this->resetTime();    
    for (size_t epoch_idx = 0 ; epoch_idx < change_times_.size(); epoch_idx++){
        this->previous_base.push_back( (double)0 );
        double top_t = epoch_idx == (change_times_.size() -1) ? change_times_[change_times_.size()-1] : change_times_[epoch_idx+1];
        //double lag_i =  double(4) / this->recombination_rate() / top_t ; 
        double lag_i = this->const_lag_ > 0 ? this->const_lag_ : double(4) / this->recombination_rate() / top_t ; 
        cout<<"lag_i = " << lag_i<<endl;
        this->lags.push_back( lag_i );
        }    
    this->resetTime();
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
    model->resetTime();
    this ->resetTime();
    for (size_t epoch_idx = 0; epoch_idx < change_times_.size(); epoch_idx++){
        for (size_t pop_j = 0 ; pop_j < this->population_number(); pop_j++ ){
            model->addPopulationSize(this->change_times_[epoch_idx], pop_j, this->total_weighted_coal_opportunity[epoch_idx][pop_j] / this->total_coal_count[epoch_idx][pop_j] /2 ,false, false);    
            }
        }
    this->check_model_updated_Ne( model );
    }


void CountModel::reset_recomb_rate ( Model *model ){
    this->compute_recomb_rate();
    model->setRecombinationRate( this->inferred_recomb_rate , false, false, 0);
    cout << " set recombination rate " << model->recombination_rate(0) << "("<<this->recomb_count_<<" / "<< this->recomb_opportunity_ << ")" <<endl;
    }

void CountModel::reset_mig_rate ( Model *model ) {
    if (this->has_migration() == false) return;
    this->compute_mig_rate();
    model->has_migration_ = false;
    
    assert( this->print_mig_rate (model->mig_rates_list_) );
    assert( this->print_mig_rate (model->total_mig_rates_list_) );
    this->initialize_mig_rate ( model->mig_rates_list_ );
    this->initialize_mig_rate ( model->total_mig_rates_list_ );
    assert( this->print_mig_rate (model->mig_rates_list_) );
    assert( this->print_mig_rate (model->total_mig_rates_list_) );

    model->resetTime();
    this ->resetTime();
    for (size_t time_i = 0; time_i < change_times_.size(); time_i++){
        for (size_t pop_j = 0 ; pop_j < this->population_number(); pop_j++ ){
            for (size_t pop_k = 0 ; pop_k < this->population_number(); pop_k++ ) {
                if ( pop_j == pop_k) continue;
                model->addMigrationRate(this->change_times_[time_i], pop_j, pop_k, this->inferred_mig_rate[time_i][pop_j][pop_k], false, false);
                }
            }
        }
    this->check_model_updated_mig (model);
    model->finalize();
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
        } 
    else {
        if (print){
            this->print_Time_count_pop();            
            }
        return;
        }
    }
    

void CountModel::compute_mig_rate(){
    this->resetTime();
    for (size_t epoch_idx = 0; epoch_idx<change_times_.size(); epoch_idx++){
        for (size_t pop_i = 0 ; pop_i < this->population_number(); pop_i++ ){
            for (size_t pop_j = 0 ; pop_j < this->population_number(); pop_j++ ){
                this->inferred_mig_rate[epoch_idx][pop_i][pop_j] = this->total_mig_count[epoch_idx][pop_i][pop_j] / this->total_weighted_mig_opportunity[epoch_idx][pop_i];
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
        this->recomb_opportunity_ += this->total_weighted_recomb_opportunity[epoch_idx][0] ;
        this->recomb_count_ += this->total_recomb_count[epoch_idx][0];
        }
    this->inferred_recomb_rate = this->recomb_count_ / this->recomb_opportunity_;
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
            
void CountModel::extract_and_update_count(ParticleContainer &Endparticles, double current_base, bool end_data ) {
    // loop over all particles
    for (size_t i = 0; i < Endparticles.particles.size(); i++) {
        ForestState* thisState = Endparticles.particles[i];
        double weight = thisState->weight();
        //cout<< " weight is "<<weight<<endl;
        // loop over all epochs
        for (size_t epoch_idx = 0; epoch_idx < this->change_times_.size(); epoch_idx++) {
            // calculate the required lagging for this epoch; don't use lagging for the final interval
            double lagging = end_data ? 0 : lags[epoch_idx];
            double x_end = current_base - lagging;
            // update counts, remove pointers to events that are processed, and remove events when reference count goes to 0
            this->update_star_count( thisState->CoaleventContainer[ epoch_idx ],   weight, x_end, this->total_coal_count[ epoch_idx ],   this->total_weighted_coal_opportunity[ epoch_idx ] );
            this->update_star_count( thisState->RecombeventContainer[ epoch_idx ], weight, x_end, this->total_recomb_count[ epoch_idx ], this->total_weighted_recomb_opportunity[ epoch_idx ] );
            this->update_migration_count( thisState->MigreventContainer[ epoch_idx ], weight, x_end, epoch_idx );
            }
        }
    }
        

void CountModel::update_star_count( deque < Starevent *> & StareventContainer_i, double weight, size_t x_end, vector<double>& total_star_count, vector<double>& total_star_opportunity ) {
    // Go through the events, starting from the leftmost and going up to x_end, and add events (weighted by weight) to the appropriate counters
    // When processed remove the event pointer from the deque; remove the event itself if its reference count becomes 0
    while (StareventContainer_i.size() > 0 && StareventContainer_i[0]->base() <= x_end) { // DEBUG changed "<" to "<=" ???
        Starevent * current_Starevent = StareventContainer_i[0];
        total_star_count[current_Starevent->pop_i()]       += weight * current_Starevent->num_event();
        total_star_opportunity[current_Starevent->pop_i()] += weight * current_Starevent->opportunity();            
        current_Starevent->pointer_counter_ --;
        if (current_Starevent->pointer_counter_ == 0) {
            delete current_Starevent;
            }
        StareventContainer_i.pop_front();
        }
    }

void CountModel::update_migration_count( deque < Migrevent *> & MigreventContainer_i, double weight, size_t x_end, size_t epoch_idx ) {
    // Go through the events, starting from the leftmost and going up to x_end, and add events (weighted by weight) to the appropriate counters
    // When processed remove the event pointer from the deque; remove the event itself if its reference count becomes 0
    while (MigreventContainer_i.size() > 0 && MigreventContainer_i[0]->base() < x_end) {
        Migrevent * current_Migrevent = MigreventContainer_i[0];
        if (current_Migrevent->event_state() == EVENT) {
            total_mig_count[epoch_idx][current_Migrevent->pop_i()][current_Migrevent->mig_pop()] += weight * current_Migrevent->num_event();
        }
        total_weighted_mig_opportunity[epoch_idx][current_Migrevent->pop_i()] += weight * current_Migrevent->opportunity();            
        current_Migrevent->pointer_counter_ --;
        if (current_Migrevent->pointer_counter_ == 0) {
            delete current_Migrevent;
            }
        MigreventContainer_i.pop_front();
        }
    }

