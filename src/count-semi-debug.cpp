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


void CountModel::print_Time_count_pop(){
    cout << " ### " ;
    for (size_t pop_i = 0 ; pop_i < this->population_number(); pop_i++ ){
        cout << " | " << std::setw(15) << "time" << std::setw(15) << "count" << std::setw(15) << "popsize" ;
        } cout<<endl;
    this->resetTime();
    for (size_t time_layer_i = 0; time_layer_i<change_times_.size(); time_layer_i++){
        cout << " ### " ;
        for (size_t pop_i = 0 ; pop_i < this->population_number(); pop_i++ ){
             cout << " | " << std::setw(15) << this->total_weighted_coal_opportunity[time_layer_i][pop_i]
                           << std::setw(15) << this->total_coal_count[time_layer_i][pop_i]
                           << std::setw(15) << this->population_size(pop_i);
            }  cout<<endl;
        // move the time_interval temp variable to the next level 
        if ( current_time_idx_ == change_times_.size() - 1) break;  
        this->increaseTime(); 
        }   cout<<endl;
    }



void CountModel::check_model_updated_Ne(Model * model){
    cout << " check_model_updated_Ne " << endl;
    model->resetTime();
    size_t time_i = 0;
    while (model->getNextTime() < FLT_MAX){
        cout << "Updated Ne at time " << setw(8) <<  model->getCurrentTime() <<" : " ;
        for (size_t pop_j = 0 ; pop_j < model->population_number() ; pop_j++){
            cout << " | " << setw(12) << model->population_size(pop_j) << " ("<< setw(12) << this->total_coal_count[time_i][pop_j] <<")" ;
        }  cout<<endl;
        model->increaseTime();
        time_i++;
        }
    
    cout << "Updated Ne at time " << setw(8) <<  model->getCurrentTime() <<" : " ;
    for (size_t pop_j = 0 ; pop_j < model->population_number() ; pop_j++){
        cout << " | " << setw(12) << model->population_size(pop_j) << " ("<< setw(12) << this->total_coal_count[time_i][pop_j] <<")" ;
        } cout<<endl;
    
    return ;                
    }



