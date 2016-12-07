/*
 * smcsmc is short for particle filters for ancestral recombination graphs. 
 * This is a free software for demographic inference from genome data with particle filters. 
 * 
 * Copyright (C) 2013-2017 Donna Henderson, Sha (Joe) Zhu and Gerton Lunter
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


#include "usage.hpp"    

//void process(struct rusage *p, char *when)
void process(struct rusage *p, double site)
{
    printf("\n At site %f ,", site );
    printf("maximum resident set size: %8zu kilobytes, ",  p->ru_maxrss          );
}


void printCurrentUsage(int &who, struct rusage *p, double averageNumNode, int data_base_at){
    getrusage(who,p);
    cout << "Usage at " << data_base_at << " : ";
    printf("page reclaims: %8zu  kilobytes, ",  p->ru_minflt          );
    printf("maximum resident set size: %8zu kilobytes, ",  p->ru_maxrss          );    
    cout << "average number of nodes: " <<  averageNumNode << endl;
}

pfTime::pfTime(time_t initial_time){
    this->initial_time_ = initial_time;
    for (size_t i=0;i<3;i++){
        this->set_time(0, i);
    }
}

void pfTime::stopwatch_start(){
    this->tmp_time_ = time(0);
}

void pfTime::stopwatch_end(size_t position){
    this->add_time(time(0) - this->tmp_time_, position);
}

