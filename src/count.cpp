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
    this->total_coal_count.clear();
    this->total_weighted_coal_opportunity.clear();
    this->total_recomb_count.clear();
    this->total_weighted_recomb_opportunity.clear();
    this->total_mig_count.clear();
    this->total_weighted_mig_opportunity.clear();
    
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

    for (size_t pop_i = 0 ; pop_i < this->population_number(); pop_i++ ){
        vector <double> tmp_count(this->population_number(), 0);
        this->total_mig_count.push_back(tmp_count);
        
        vector <double> tmp_opportunity(this->population_number(), 1);
        this->total_weighted_mig_opportunity.push_back(tmp_opportunity);
        }
    //this->print_Time_count_pop();
    
    //cout << endl;    
    this->resetTime();
    return ;
    }


void CountModel::reset_model_Ne(Model * model, bool online, bool print){
    if ( !online ){ 
        if (print){
            this->print_Time_count_pop();            
            }
        return;
        } 
    else {
        model->resetTime();
        this ->resetTime();

        for (size_t time_i = 0; time_i < change_times_.size(); time_i++){
            for (size_t pop_j = 0 ; pop_j < this->population_number(); pop_j++ ){
                //this->count_events_in_one_interval(Endparticles, time_i, pop_j, x_start, x_end);
                model->addPopulationSize(this->change_times_[time_i], pop_j, this->total_weighted_coal_opportunity[time_i][pop_j] / this->total_coal_count[time_i][pop_j] /2 ,false, false);    
                }
            }
        //this->check_model_updated_Ne((this)); //For checking only
         cout<<" MODEL IS RESET "<<endl;
                
        this->compute_recomb_rate();
        this->compute_mig_rate();
        this->check_model_updated_Ne(model);
        }
    }
    


/*! 
 * Compute the recombination rate once we have sweeped through all the data, and recorded the recomb_opportunity and recomb_counts
 */ 
void CountModel::compute_recomb_rate () {
    this->resetTime();
    double scaling_pop_size_N0 = this->population_size();
    double recomb_opportunity = 0;
    double recomb_count = 0;
    
    for (size_t time_layer_i = 0; time_layer_i<change_times_.size(); time_layer_i++){
        for (size_t pop_j = 0 ; pop_j < this->population_number(); pop_j++ ){
            double pop_ratio = this->population_size(pop_j) / scaling_pop_size_N0;          
            recomb_opportunity += this->total_weighted_recomb_opportunity[time_layer_i][pop_j] / pop_ratio ;
            recomb_count += this->total_recomb_count[time_layer_i][pop_j] / pop_ratio ;
            }
        // Advance to the next interval level
        if ( current_time_idx_ == change_times_.size() - 1) break;  
        this->increaseTime(); 
        }
    this->inferred_recomb_rate = recomb_count / recomb_opportunity;
    // reset the inferred recombination rate in the current Model and CountModel?
    }

                
size_t CountModel::find_time_interval (double start_height, double end_height){
    double pop_start_height = 0; 
    for (size_t time_i = 0; time_i < change_times_.size() ; time_i++){
        double pop_end_height = (time_i == this->change_times_.size()-1) ? FLT_MAX : change_times_[time_i+1];
            if ( pop_start_height <= start_height && end_height <= pop_end_height ) {
                return time_i;
                }
        pop_start_height = pop_end_height;
        } // End of for loop: < change_times.size()  
    } // End of function: CountModel::find_time_interval( ... )



void CountModel::extract_and_update_count(ParticleContainer &Endparticles, double x_start, double x_end){
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
        while (counting_state -> current_base() >= x_end){
            counting_state = counting_state->previous_state;                
        }
        
        // Start counting
        while ( counting_state -> current_base() >= x_start ) {  // Making sure there is coalescent events between the interval

            for ( size_t j = 0; j < counting_state->CoaleventContainer.size(); j++){
                Coalevent * current_Coalevent = counting_state->CoaleventContainer[j];

                /*! Cumulate the coalescent events if the event is within the interval 
                 */
                size_t time_i = this->find_time_interval (current_Coalevent->start_height(),  current_Coalevent->end_height());
                 
                this->total_coal_count[time_i][current_Coalevent->pop_i()]                    += weight * current_Coalevent->num_event();
                this->total_weighted_coal_opportunity[time_i][current_Coalevent->pop_i()]     += weight * current_Coalevent->opportunity();            
                    
                } //  < counting_state->CoaleventContainer.size()                        

            //for ( size_t j = 0; j < counting_state->CoaleventContainer.size(); j++){
                //Coalevent * current_Coalevent = counting_state->CoaleventContainer[j];

                ///*! Cumulate the coalescent events if the event is within the interval 
                 //*/
                //size_t time_i = this->find_time_interval (current_Coalevent->start_height(),  current_Coalevent->end_height());
                 
                //if (current_Coalevent->event_state() == COAL_EVENT || current_Coalevent->event_state() == COAL_NOEVENT){                        
                    //this->total_coal_count[time_i][current_Coalevent->pop_i()]                    += weight * current_Coalevent->num_coal();
                    //this->total_weighted_coal_opportunity[time_i][current_Coalevent->pop_i()]     += weight * current_Coalevent->opportunity();            
                    //}
                //if (current_Coalevent->event_state() == REC_EVENT || current_Coalevent->event_state() == REC_NOEVENT){
                    //this->total_recomb_count[time_i][current_Coalevent->pop_i()]                  += weight * current_Coalevent->num_recomb();
                    //this->total_weighted_recomb_opportunity[time_i][current_Coalevent->pop_i()]   += weight * current_Coalevent->opportunity();
                    //}
                //if ((current_Coalevent->event_state() == MIGR_EVENT || current_Coalevent->event_state() == MIGR_NOEVENT)){                        
                    //if (current_Coalevent->event_state() == MIGR_EVENT){
                        //this->total_mig_count[current_Coalevent->pop_i()][current_Coalevent->mig_pop()] += current_Coalevent->num_mig() * weight;
                        //} 
                    //for (size_t potential_pop = 0; potential_pop < this->total_weighted_mig_opportunity[current_Coalevent->pop_i()].size(); potential_pop++){
                        //this->total_weighted_mig_opportunity[current_Coalevent->pop_i()][potential_pop] += current_Coalevent->opportunity() * weight;    
                        //}
                    //}
                    
                //} //  < counting_state->CoaleventContainer.size()                        

                
                counting_state = counting_state->previous_state;      
                if (!counting_state) break;                      
                
            }  // End of while loop: counting_state -> current_base() >= x_start
        } //  End of for loop: < Endparticles.particles.size()
        return;
    } // 

