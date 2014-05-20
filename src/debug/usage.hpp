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


#include <stdlib.h>
#include <sys/resource.h>
#include <unistd.h>
#include <stdio.h>
#include <sys/time.h>
#include <iostream>

using namespace std;

#ifndef USAGE
#define USAGE

/*! @ingroup group_resource */
//void process(struct rusage *p, char *when);
void process(struct rusage *p, double site);
/*! @ingroup group_resource */
void printCurrentUsage(int &who, struct rusage *p, double averageNumNode,int data_base_at);

/*! @ingroup group_resource */
class pfTime{
    public:

        /*!
         * Constructors and Destructors
         */      
        pfTime();
        pfTime(time_t initial_time);
        ~pfTime(){}
        
        /*!
         * Methods
         */ 
        void add_time(int t, size_t position){this->timing_[position] += t;};    
        void stopwatch_start();
        void stopwatch_end(size_t position);
    
        /*!
         *  Setters and getters:
         */ 
        void set_time(int t, size_t position){this->timing_[position] = t;};    
        void the_end(){this->end_time_ = time(0);}

        /*!
         * Members
         */     
        int timing_[3];
            //timing_[0] time for building particle initial states
            //timing_[1] time for resampling
            //timing_[2] time for update 
    
    private:
    
        /*!
         * Members
         */     
        time_t initial_time_;
        time_t end_time_;
        time_t tmp_time_;
};

#endif
