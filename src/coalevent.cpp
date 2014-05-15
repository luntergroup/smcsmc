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
 * \class Starevent
 *  another constructor for rec_point, to be called from sample the next genealogy
 */ 


Starevent::Starevent(){
    this->init();
    }
    
/*! \brief Copy the Starevent from the previous ForestState*/    
Starevent::Starevent(Starevent* previous_Starevent){ 
    this->pop_i_        = previous_Starevent->pop_i();
    //this->start_height_ = previous_Starevent->start_height();
    //this->end_height_   = previous_Starevent->end_height(); 
    this->num_event_    = previous_Starevent->num_event();
    this->opportunity_  = previous_Starevent->opportunity();
    this->event_state_  = previous_Starevent->event_state();
    
    this->change_time_i_ = previous_Starevent->change_time_i();
    //this->counted_ = previous_Starevent->counted();
    this->base_ = previous_Starevent->base();
    }
    

/*! \brief  Initialize Starevent */    
Starevent::Starevent( 
           size_t pop_i,
           //double start_time,
           //double end_time, 
           double opportunity,
           eventCode event_code ) {
    this->init();    
    this->set_pop_i(pop_i);
    //this->set_start_height( start_time );
    //this->set_end_height ( end_time );
    this->set_opportunity( opportunity );
    this->set_event_state(event_code);
    this->set_num_event ( this->event_state() == EVENT ? 1 : 0);
    }


void Starevent::init(){
    this->set_pop_i(0);
    //this->set_start_height(0);
    //this->set_end_height(0);    
    this->set_num_event(0);
    this->set_opportunity(0);
    this->set_event_state(INIT_NULL);
    
    this->set_change_time_i ( 0 );
    //this->set_counted ( false );
    this->set_base ( 0 );

    }
