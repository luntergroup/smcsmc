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


#include"coalevent.hpp"
/*!
 * another constructor for rec_point, to be called from sample the next genealogy
 */ 

Coalevent::Coalevent(){
    this->init();
    }
    
    
Coalevent::Coalevent(Coalevent* previous_Coalevent){ /*!< Copy the Coalevent from the previous ForestState*/
    this->pop_i_        = previous_Coalevent->pop_i();
    this->mig_pop_      = previous_Coalevent->mig_pop();
    this->start_height_ = previous_Coalevent->start_height();
    this->end_height_   = previous_Coalevent->end_height(); 
    this->num_coal_     = previous_Coalevent->num_coal();
    this->num_mig_      = previous_Coalevent->num_mig();
    this->num_recomb_   = previous_Coalevent->num_recomb();
    this->opportunity_  = previous_Coalevent->opportunity();
    this->event_state_  = previous_Coalevent->event_state();;
    }
    

Coalevent::Coalevent( 
           size_t pop_i,
           size_t mig_pop,
           double start_time,
           double end_time, 
           double opportunity,
           eventCode event_code ) {
    this->init();
    
    this->set_pop_i(pop_i);
    this->set_mig_pop_i(mig_pop);
    this->set_start_height( start_time );
    this->set_end_height ( end_time );
    this->set_opportunity( opportunity );
    this->set_event_state(event_code);

    switch (this->event_state()) {
        case COAL_EVENT:        this->set_num_coal(1);   break;
        case REC_EVENT:         this->set_num_recomb(1); break;
        case MIGR_EVENT:        this->set_num_mig(1);    break;
        default: /*no events!*/ this->set_num_coal(0);
                                this->set_num_coal(0);
                                this->set_num_coal(0);   break;    
        }

    }


void Coalevent::init(){
    this->set_pop_i(0);
    this->set_mig_pop_i(-1);
    this->set_start_height(0);
    this->set_end_height(0);    
    this->set_num_recomb(0);
    this->set_num_mig(0);
    this->set_num_coal(0);
    this->opportunity_ = 0;
    this->set_event_state(INIT_NULL);
    }
