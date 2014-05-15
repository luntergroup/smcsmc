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


#include "scrm/forest.h"

#ifndef STAREVENT
#define STAREVENT

/*!
 * \brief Used for recording the number and the time intervals of events between two ForestState 
 */
class Starevent{ 
    friend class CountModel;
    friend class ForestState;
    friend class ParticleContainer;
    friend class Coalevent;
    friend class Migrevent;
    friend class Recombevent;
    public:
        //
        // Constructors and Destructors
        //    
        Starevent();
        Starevent(Starevent* previous_Starevent);
        Starevent(size_t pop_i, 
                  //double start_time,
                  //double end_time, 
                  double opportunity,
                  eventCode event_code );
        ~Starevent(){};
    private:
        //
        // Methods
        //   
        void init();
        void clear();

        //
        // Getters and setters
        //
        size_t pop_i() const {return this->pop_i_; } 
        void set_pop_i( size_t i ){ this->pop_i_ = i; }
                 
        //double start_height() const { return this->start_height_; };
        //void set_start_height(double start_height) { this->start_height_ = start_height; };
        
        //double end_height() const { return this->end_height_; };
        //void set_end_height(double end_height) { this->end_height_ = end_height; };
        
        size_t num_event () const { return this->num_event_;}
        void  set_num_event ( size_t num ){ this->num_event_ = num; }    
        
        double opportunity() const {return this->opportunity_;};
        void set_opportunity( double value ) { this->opportunity_ = value;}
                
        eventCode event_state() const { return this->event_state_ ;}
        void set_event_state(eventCode state){ this->event_state_ = state; }
        
        size_t change_time_i() const { return this->change_time_i_; }
        void set_change_time_i( size_t i ){ this->change_time_i_ = i ; }
        
        //bool counted() const { return this->counted_ ; }
        //void set_counted (bool counted) { this->counted_ = counted; }
        
        double base() const { return this->base_; }
        void set_base ( double base ) { this->base_ = base; }
        
        //
        // Members
        //   
        size_t change_time_i_;
        size_t pop_i_;                
        //double start_height_;
        //double end_height_;                 
        double opportunity_;        
        size_t num_event_;        
        eventCode event_state_; 
        //bool counted_;
        double base_;
    };


/*!
 * \brief Derived class of Starevent, recording the number and the time intervals of Coalescent events between two ForestState 
 */
class Coalevent : public Starevent{
    //friend class CountModel;
    //friend class ForestState;
    //friend class ParticleContainer;
    public:
        Coalevent(size_t pop_i, 
                  //double start_time,
                  //double end_time, 
                  double opportunity,
                  eventCode event_code ) 
                 : Starevent ( pop_i, 
                    //start_time, 
                    //end_time, 
                    opportunity, 
                    event_code ){ };  
        ~Coalevent(){};

    };
    
/*!
 * \brief Derived class of Starevent, recording the number and the time intervals of Recombination events between two ForestState 
 */    
class Recombevent : public Starevent{
    //friend class CountModel;
    //friend class ForestState;
    //friend class ParticleContainer;
    public:
        Recombevent(size_t pop_i, 
                    //double start_time,
                    //double end_time, 
                    double opportunity,
                    eventCode event_code ) 
                   : Starevent ( pop_i, 
                      //start_time, 
                      //end_time, 
                      opportunity, 
                      event_code ){ };
        ~Recombevent(){};
    };
    

/*!
 * \brief Derived class of Starevent, recording the number and the time intervals of Migration events between two ForestState 
 */    
class Migrevent : public Starevent{
    //friend class CountModel;
    //friend class ForestState;
    //friend class ParticleContainer;
    public:
        Migrevent(size_t pop_i, 
                  size_t mig_pop,
                  //double start_time,
                  //double end_time, 
                  double opportunity,
                  eventCode event_code )
                 : Starevent ( pop_i, 
                    //start_time, 
                    //end_time, 
                    opportunity, 
                    event_code ){ this->set_mig_pop_i (mig_pop); };
        ~Migrevent(){};
        
        size_t mig_pop() const {return this->mig_pop_; } 
        void set_mig_pop_i( size_t i ){ this->mig_pop_ = i; }
    
        size_t mig_pop_;
    };

#endif
