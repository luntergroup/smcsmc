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

void CountModel::print_change_time(){
    this->resetTime();
    cout<<"change time size = "<< this->change_times_.size()<<endl;
    cout << "change_times at: ";
    while ( this->getNextTime() < FLT_MAX ){
        cout << setw(8) << this->getCurrentTime() << setw(3) << "(" << this->current_time_idx_ << ") --- ";
        this->increaseTime();
        }
    cout << setw(8) << this->getCurrentTime() << setw(3) << "(" << this->current_time_idx_ << ") --- Infinity" << endl;
    }


void CountModel::print_pop_size(){
    for (size_t pop_i = 0; pop_i < this->population_number(); pop_i++ ){
        this->resetTime();
        cout << "Popsize (" << pop_i+1 << "):     "  ;
        do {            
            cout << setw(8) << population_size(pop_i) << "      --- ";
            if ( current_time_idx_ == change_times_.size() - 1) break;              
            this->increaseTime();                    
            } while (this->getCurrentTime() <= FLT_MAX);        
        cout<<endl;
        }
    }








