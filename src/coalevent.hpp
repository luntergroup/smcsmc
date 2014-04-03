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


#include "scrm/time_interval.h"

#ifndef COALEVENT
#define COALEVENT

/*!
 * \class Coalevent is used for recording the number and the time intervals of the recombination and coalescent events between two ARGs 
 */

class Coalevent{
    public:
        
        /*!
         * Constructors and Destructors
         */  
        Coalevent();
        Coalevent(Coalevent* previous_Coalevent);
        Coalevent(size_t pop_i, 
                  size_t mig_pop,
                  double start_time,
                  double end_time, 
                  double opportunity,
                  eventCode event_code );
        ~Coalevent(){};

        /*!
         * Methods
         */ 
        void init();
        void clear();
        
        /*! 
         * Getters and setters
         */ 
        size_t pop_i() const {return this->pop_i_; } 
        void set_pop_i( size_t i ){ this->pop_i_ = i; }
        
        size_t mig_pop() const {return this->mig_pop_; } 
        void set_mig_pop_i( size_t i ){ this->mig_pop_ = i; }
         
        double start_height() const { return this->start_height_; };
        void set_start_height(double start_height) { this->start_height_ = start_height; };
        
        double end_height() const { return this->end_height_; };
        void set_end_height(double end_height) { this->end_height_ = end_height; };
    
        size_t num_coal() const { return this->num_coal_ ;};
        void set_num_coal(size_t num_coal) { this->num_coal_ = num_coal; };
        
        size_t num_recomb() const { return this->num_recomb_ ;};
        void set_num_recomb(size_t num_recomb) { this->num_recomb_ = num_recomb; };
        
        size_t num_mig() const { return this->num_mig_ ;};
        void set_num_mig(size_t num_mig) { this->num_mig_ = num_mig; };
        
        double opportunity() const {return this->opportunity_;};
        void set_opportunity( double value ) { this->opportunity_ = value;}
                
        eventCode event_state() const { return this->event_state_ ;};
        void set_event_state(eventCode state){ this->event_state_ = state; };
        
    private:

        /*!
         * Members
         */ 
        size_t pop_i_;
        size_t mig_pop_;
        
        double start_height_;
        double end_height_; 
                
        double opportunity_;
        
        size_t num_coal_;
        size_t num_recomb_;
        size_t num_mig_;
        
        eventCode event_state_;
};

#endif
