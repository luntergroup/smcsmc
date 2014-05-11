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
//#include <limits>       // std::numeric_limits


void CountModel::init(){
    this->init_coal_and_recomb();
    this->init_migr();
    this->init_lags();
    return ;
    }


void CountModel::init_coal_and_recomb(){
    this->total_coal_count.clear();
    this->total_weighted_coal_opportunity.clear();
    this->total_recomb_count.clear();
    this->total_weighted_recomb_opportunity.clear();

    this->resetTime();    
    for (size_t time_layer_i = 0 ; time_layer_i < change_times_.size(); time_layer_i++){
        vector <double> tmp_count(this->population_number(), 0);
        this->total_coal_count.push_back(tmp_count);
        this->total_recomb_count.push_back(tmp_count);
        for (size_t pop_i = 0 ; pop_i < this->population_number(); pop_i++ ){
            this->total_coal_count[time_layer_i][pop_i] = 1 / ( 2 * population_size() );
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

    for (size_t pop_i = 0 ; pop_i < this->population_number(); pop_i++ ){
        vector <double> tmp_count(this->population_number(), 0);
        this->total_mig_count.push_back(tmp_count);
        this->inferred_mig_rate.push_back(tmp_count);
        
        vector <double> tmp_opportunity(this->population_number(), 1);
        this->total_weighted_mig_opportunity.push_back(tmp_opportunity);
        }        
    }


void CountModel::init_lags(){
    this->previous_base.clear();
    this->lags.clear();        
    this->resetTime();    
    for (size_t time_layer_i = 0 ; time_layer_i < change_times_.size(); time_layer_i++){
        this->previous_base.push_back( (double)0 );
        double top_t = time_layer_i == (change_times_.size() -1) ? change_times_[change_times_.size()-1] : change_times_[time_layer_i+1];
        //double lag_i =  double(4) / this->recombination_rate() / top_t ; 
        double lag_i = this->const_lag_ > 0 ? this->const_lag_ : double(4) / this->recombination_rate() / top_t ; 
        cout<<"lag_i = " << lag_i<<endl;
        this->lags.push_back( lag_i );
        }    
    this->resetTime();
    }


void CountModel::reset_Ne ( Model *model ){
    model->resetTime();
    this ->resetTime();
    for (size_t time_i = 0; time_i < change_times_.size(); time_i++){
        for (size_t pop_j = 0 ; pop_j < this->population_number(); pop_j++ ){
            model->addPopulationSize(this->change_times_[time_i], pop_j, this->total_weighted_coal_opportunity[time_i][pop_j] / this->total_coal_count[time_i][pop_j] /2 ,false, false);    
            }
        }
    this->check_model_updated_Ne( model );
    }


void CountModel::reset_recomb_rate ( Model *model ){
    this->compute_recomb_rate();
    model->rec_rate_ = this->inferred_recomb_rate;
    cout << " set recombination rate " << model->rec_rate_ << "("<<this->recomb_count_<<")" <<endl;
    }


void CountModel::initialize_mig_rate ( vector <vector<double>*> & rates_list ){
    for (size_t i = 0; i < rates_list.size(); i++ ){
        for (size_t j = 0; j < rates_list[i]->size() ; j++){
            rates_list[i]->at(j) = 0;
            }
        }    
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
                model->addMigrationRate(this->change_times_[time_i], pop_j, pop_k, this->inferred_mig_rate[pop_j][pop_k], false, false);
                }
            }
        }
        
    this->check_model_updated_mig (model);

    model->finalize();
    assert( this->print_mig_rate (model->mig_rates_list_) );
    assert( this->print_mig_rate (model->total_mig_rates_list_) );
    
    }



void CountModel::reset_model_parameters(Model * model, bool online, bool print){
    if ( !online ){ 
        if (print){
            this->print_Time_count_pop();            
            }
        return;
        } 
    else {
        cout<<" MODEL IS RESET "<<endl;
        this->reset_recomb_rate ( model );
        this->reset_Ne ( model );
        this->reset_mig_rate ( model );
        }
    }
    

void CountModel::compute_mig_rate(){
    for (size_t pop_i = 0 ; pop_i < this->population_number(); pop_i++ ){
        for (size_t pop_j = 0 ; pop_j < this->population_number(); pop_j++ ){
            this->inferred_mig_rate[pop_i][pop_j] = this->total_mig_count[pop_i][pop_j] / this->total_weighted_mig_opportunity[pop_i][pop_j];
            }
        }
    }


/*! 
 * \brief Compute the recombination rate once we have sweeped through all the data, and recorded the recomb_opportunity and recomb_counts
 */ 
void CountModel::compute_recomb_rate () {
    this->resetTime();
    double scaling_pop_size_N0 = this->population_size();
    double recomb_opportunity = 0;
    //double recomb_count = 0;
    this->recomb_count_ = 0;
    
    for (size_t time_layer_i = 0; time_layer_i<change_times_.size(); time_layer_i++){
        for (size_t pop_j = 0 ; pop_j < this->population_number(); pop_j++ ){
            double pop_ratio = this->population_size(pop_j) / scaling_pop_size_N0;          
            recomb_opportunity += this->total_weighted_recomb_opportunity[time_layer_i][pop_j] / pop_ratio ;
            this->recomb_count_ += this->total_recomb_count[time_layer_i][pop_j] / pop_ratio ;
            }
        // Advance to the next interval level
        if ( current_time_idx_ == change_times_.size() - 1) break;  
        this->increaseTime(); 
        }
    this->inferred_recomb_rate = this->recomb_count_ / recomb_opportunity;
    // reset the inferred recombination rate in the current Model and CountModel?
    }


double CountModel::extract_and_update_count(ParticleContainer &Endparticles, double current_base, bool end_data ){
    //for ( size_t time_i = 0 ; time_i < this->change_times_.size(); time_i ++){
    for ( size_t time_i = this->change_times_.size() - 1 ; (int)time_i >= 0 ; time_i --){
        //cout << "at time level " << time_i << "current_base " << current_base << " this->lags[time_i] " << this->lags[time_i] <<endl;
        double x_end =  (double)this->previous_base[time_i] < ( current_base - this->lags[time_i] ) ? ( current_base - this->lags[time_i] ) : (double)this->previous_base[time_i] ;
        if (end_data){
            x_end = current_base;
            }
        for (size_t pop_j = 0 ; pop_j < this->population_number(); pop_j++ ){
            this->count_events_in_one_interval(Endparticles, time_i, pop_j, previous_base[time_i], x_end );                
            }
        double previous_base_tmp = current_base - lags[time_i] ;
        previous_base[time_i] = previous_base_tmp > 0 ? previous_base_tmp : (double)0;
        }  
    
    double remove_particle_before_site = previous_base[0];
    for (size_t i = 0 ; i < previous_base.size() ; i++ ){
        remove_particle_before_site = remove_particle_before_site < previous_base[i] ? remove_particle_before_site : previous_base[i];
        }
    dout<< "remove_particle_before_site = "<<remove_particle_before_site<<endl;
    return remove_particle_before_site;
    }




void CountModel::count_events_in_one_interval(ParticleContainer &Endparticles, size_t time_i, size_t pop_j, double x_start, double x_end){

    if ( x_start == x_end ){ return; }
    
    for (size_t i = 0; i < Endparticles.particles.size(); i++){
        
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
        
        ForestState* counting_state = Endparticles.particles[i];
        
        double weight = counting_state->weight();
        
        
        // Skip a few state between lagging until the most updated case
        while (counting_state -> current_base() >= x_end && counting_state->previous_state){
            counting_state = counting_state->previous_state;                
        }
        
        // Start counting
        while ( counting_state -> current_base() >= x_start ) {  // Making sure there is coalescent events between the interval
            for ( size_t event_i = 0 ; event_i < counting_state->CoaleventContainer.size() ; event_i++ ){
                Coalevent * current_Coalevent = counting_state->CoaleventContainer[event_i];
                /*! Cumulate the coalescent events if the event is within the interval 
                 */
                if ( time_i == current_Coalevent->change_time_i() && pop_j == current_Coalevent->pop_i() && !current_Coalevent->counted() ){
                    this->total_coal_count[time_i][current_Coalevent->pop_i()]                    += weight * current_Coalevent->num_event();
                    this->total_weighted_coal_opportunity[time_i][current_Coalevent->pop_i()]     += weight * current_Coalevent->opportunity();            
                    current_Coalevent->set_counted(true);
                    }
                } //  < counting_state->CoaleventContainer.size() 

            for ( size_t event_i = 0 ; event_i < counting_state->RecombeventContainer.size() ; event_i++ ){
                Recombevent * current_Recombevent = counting_state->RecombeventContainer[event_i];
                /*! Cumulate the recombination events if the event is within the interval 
                 */
                if ( time_i == current_Recombevent->change_time_i() && pop_j == current_Recombevent->pop_i() && !current_Recombevent->counted() ){
                    this->total_recomb_count[time_i][current_Recombevent->pop_i()]                    += weight * current_Recombevent->num_event();
                    this->total_weighted_recomb_opportunity[time_i][current_Recombevent->pop_i()]     += weight * current_Recombevent->opportunity();            
                    current_Recombevent->set_counted(true);
                    }
                } //  < counting_state->RecombeventContainer.size() 
                
            for ( size_t j = 0; j < counting_state->MigreventContainer.size(); j++){
                Migrevent * current_Migrevent = counting_state->MigreventContainer[j];

                /*! Cumulate the recombination events if the event is within the interval 
                 */
                if (current_Migrevent->event_state() == EVENT){
                    this->total_mig_count[current_Migrevent->pop_i()][current_Migrevent->mig_pop()] += current_Migrevent->num_event() * weight;
                    } 
                for (size_t potential_pop = 0; potential_pop < this->total_weighted_mig_opportunity[current_Migrevent->pop_i()].size(); potential_pop++){
                    this->total_weighted_mig_opportunity[current_Migrevent->pop_i()][potential_pop] += current_Migrevent->opportunity() * weight;    
                    }    
                } //  < counting_state->MigreventContainer.size() 
                                
                counting_state = counting_state->previous_state;      
                if (!counting_state) break;                      
                
            }  // End of while loop: counting_state -> current_base() >= x_start
        } //  End of for loop: < Endparticles.particles.size()
        return;
    } // 


void CountModel::count_events_in_one_interval_alt(ParticleContainer &Endparticles, size_t time_i, size_t pop_j, double x_start, double x_end){

    if ( x_start == x_end ){ return; }
    
    for (size_t i = 0; i < Endparticles.particles.size(); i++){
        
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
        
        ForestState* counting_state = Endparticles.particles[i];
        
        double weight = counting_state->weight();
        
        
        // Skip a few state between lagging until the most updated case
        while (counting_state -> current_base() >= x_end && counting_state->previous_state){
            counting_state = counting_state->previous_state;                
        }
        
        // Start counting
        while ( counting_state -> current_base() >= x_start ) {  // Making sure there is coalescent events between the interval
            while (!counting_state->CoaleventContainer.empty()){
                Coalevent * current_Coalevent = counting_state->CoaleventContainer.back();
                /*! Cumulate the coalescent events if the event is within the interval 
                 */
                if ( time_i == current_Coalevent->change_time_i() && pop_j == current_Coalevent->pop_i() ){
                    this->total_coal_count[time_i][current_Coalevent->pop_i()]                    += weight * current_Coalevent->num_event();
                    this->total_weighted_coal_opportunity[time_i][current_Coalevent->pop_i()]     += weight * current_Coalevent->opportunity();            
                    delete current_Coalevent;
                    counting_state->CoaleventContainer.pop_back();
                    }
                else {
                    break;
                    } 
                } //  < counting_state->CoaleventContainer.size() 

            while (!counting_state->RecombeventContainer.empty()){
                Recombevent * current_Recombevent = counting_state->RecombeventContainer.back();
                /*! Cumulate the recombination events if the event is within the interval 
                 */
                if ( time_i == current_Recombevent->change_time_i() && pop_j == current_Recombevent->pop_i() ){
                    this->total_recomb_count[time_i][current_Recombevent->pop_i()]                    += weight * current_Recombevent->num_event();
                    this->total_weighted_recomb_opportunity[time_i][current_Recombevent->pop_i()]     += weight * current_Recombevent->opportunity();            
                    delete current_Recombevent;
                    counting_state->RecombeventContainer.pop_back();
                    }
                else {
                    break;
                    } 
                } //  < counting_state->RecombeventContainer.size() 
                
            for ( size_t j = 0; j < counting_state->MigreventContainer.size(); j++){
                Migrevent * current_Migrevent = counting_state->MigreventContainer[j];

                /*! Cumulate the recombination events if the event is within the interval 
                 */
                if (current_Migrevent->event_state() == EVENT){
                    this->total_mig_count[current_Migrevent->pop_i()][current_Migrevent->mig_pop()] += current_Migrevent->num_event() * weight;
                    } 
                for (size_t potential_pop = 0; potential_pop < this->total_weighted_mig_opportunity[current_Migrevent->pop_i()].size(); potential_pop++){
                    this->total_weighted_mig_opportunity[current_Migrevent->pop_i()][potential_pop] += current_Migrevent->opportunity() * weight;    
                    }    
                } //  < counting_state->MigreventContainer.size() 
                                
                counting_state = counting_state->previous_state;      
                if (!counting_state) break;                      
                
            }  // End of while loop: counting_state -> current_base() >= x_start
        } //  End of for loop: < Endparticles.particles.size()
        return;
    } // 

