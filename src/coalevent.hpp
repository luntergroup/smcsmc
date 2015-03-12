/*
 * smcsmc is short for particle filters for ancestral recombination graphs. 
 * This is a free software for demographic inference from genome data with particle filters. 
 * 
 * Copyright (C) 2013, 2014 Sha (Joe) Zhu and Gerton Lunter
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


#include "forest.h"

#ifndef EventRecorder
#define EventRecorder

#ifndef NDEBUG
#define EventRecorderdout (std::cout << "    EventRecorder ")
#else
#pragma GCC diagnostic ignored "-Wunused-value"
#define EventRecorderdout 0 && (std::cout << "    EventRecorder ")
#endif


extern double recomb_opp; //DEBUG


class EvolutionaryEvent {
	// Constructors
	EvolutionaryEvent( double start_height, double end_height, double start_base, double end_base, int weight ) :
	                   start_height(start_height),
	                   end_height(end_height),
	                   start_base(start_base),
	                   end_base(end_base),
	                   weight(weight) {  // signal recombination opportunity
	                      assert((start_height <= end_height) && (start_base <= end_base) && (start_base >= 0) ); 
	                   };
	EvolutionaryEvent( double start_height, double end_height, double end_base, int population_index, int weight ) :
			   start_height(start_height),
			   end_height(end_height),
			   start_base(-1),        // signal coalescent/migration opportunity
			   end_base(end_base),
			   weight(weight),
			   a( population_index ) { 
                              assert(start_height <= end_height); 
                           };

	// Methods
	void set_recomb_event_pos( double recomb_x_position ) {
		assert( this->is_recomb() );
		assert( this->no_event() );
		assert( recomb_x_position <= end_base );
		assert( start_base <= recomb_x_position );
		event_data = -1;
		a.recomb_pos = recomb_x_position; }
	void set_recomb_event_time( double recomb_t_position ) {
		assert( this->is_recomb() );
		assert( this->no_event() );
		assert( recomb_t_position <= end_height );
		assert( start_height <= recomb_t_position );
		event_data = 0;
		a.recomb_pos = recomb_t_position; }
	void set_coal_event() {
		assert( this->is_coalmigr() );
		assert( this->no_event() );
		event_data = -1; }
	void set_migr_event( int to_population ) {
		assert( this->is_coalmigr() );
		assert( this->no_event() );
		assert( to_population != a.coal_migr_population );
		assert( to_population >= 0 );
		event_data = to_population; }
	bool is_recomb() const { return start_base >= 0; }
	bool is_coalmigr() const { return start_base < 0; }
	bool no_event() const { return event_data == -2; }
	bool is_recomb_event() const { assert(is_recomb()); return event_data > -2; }
	bool is_coal_event() const { assert(is_coalmigr()); return event_data == -1; }
	bool is_migr_event() const { assert(is_coalmigr()); return event_data >= 0; }
	double coalescence_opportunity() const {
		assert (is_coalmigr());
		return weight * (end_height - start_height); }
	double migration_opportunity() const {
		assert (is_coalmigr());
		return end_height - start_height; }
	double recombination_opportunity() const {
		assert (is_recomb());
		return (end_height - start_height) * (end_base - start_base); }
	int get_population() const {
		assert (is_coalmigr());
		return a.coal_migr_population; }
	int get_migr_to_population() {
		assert (is_coalmigr());
		assert (event_data >= 0);
		return event_data; }
	bool decrease_refcount_is_zero() {
		assert( ref_counter-- > 0 );
		return (ref_counter == 0); }
	void increase_refcount() { ref_counter++; }
	bool print_event();
        bool append_event( const EvolutionaryEvent& e );
	              
	// Members
private:
	double start_height;
	double end_height;
        double start_base;     // Recombinations: determines (w/end_base) the x-extent of recomb. opportunity.  For coal/migr, <0
	double end_base;
        int event_data {-2};   // -2 == no event; otherwise type-specific meaning:
	                       // recomb:    -1 == event at top edge (time-wise sampling)
	                       //             0 == event at right-hand edge (sequence-wise sampling)
	                       // coal/migr: -1 == coalescent event;
	                       //             0..: migration to this population index
	int weight;            // number of lineages (for recombination and coalescence, not migration) contributing to opportunity
        int ref_counter {1};
	union A {
	  int    coal_migr_population; // index of population of branch that is migrating or coalescing
	  double recomb_pos;           // position of recombination event, along right-hand (end_base) or top (end_height) edge
	  A( int p ): coal_migr_population( p ) {}
	  A(): coal_migr_population( 0 ) {}
	} a;
};


/*!
 * \brief Used for recording the number and the time intervals of events between two ForestState 
 */
class Coalevent{ 
    friend class CountModel;
    friend class ForestState;
    friend class Migrevent;
    friend class Recombevent;
    
    public:
        ~Coalevent(){ };
    //
    protected:
        // Methods
        void init();
        void clear();

        // Getters and setters
        size_t pop_i() const {return this->pop_i_; } 
        void set_pop_i( size_t i ){ this->pop_i_ = i; }
                 
        //double start_height() const { return this->start_height_; };
        //void set_start_height(double start_height) { this->start_height_ = start_height; };
        
        //double end_height() const { return this->end_height_; };
        //void set_end_height(double end_height) { this->end_height_ = end_height; };
        
        size_t num_event () const { return this->num_event_;}
        void  set_num_event ( size_t num ){ this->num_event_ = num; }    
        
        //double opportunity() const {return this->opportunity_;};
        //void set_opportunity( double value ) { this->opportunity_ = value;}
        
        double opportunity_y() const { return this->opportunity_y_; }
        void set_opportunity_y ( double value ) { this->opportunity_y_ = value; }

        eventCode event_state() const { return this->event_state_ ;}
        void set_event_state(eventCode state){ this->event_state_ = state; }
        
        size_t epoch_index() const { return this->epoch_index_; }
        void set_epoch_index( size_t i ){ this->epoch_index_ = i ; }
                
        double end_base() const { return this->end_base_; }
        void set_end_base ( double base ) { this->end_base_ = base; }        
    private:
        Coalevent(size_t pop_i, 
                  //double start_time,
                  //double end_time, 
                  double opportunity_y, eventCode event_code, double end_base );
                  
        double opportunity() const {return this->opportunity_y_; };
        bool print_event();
        // Members
        size_t epoch_index_;
        size_t pop_i_;                
        //double start_height_;
        double end_height_;
        //double opportunity_;
        double opportunity_y_;
        size_t num_event_;        
        eventCode event_state_; 
        double end_base_;

        int pointer_counter_; // Number of ForestStates are pointing at this event
    };


class Recombevent : public Coalevent{
    friend class CountModel;
    friend class ForestState;
    //public:
    Recombevent( size_t pop_i, double opportunity_y, eventCode event_code, double start_base, double end_base)
             : Coalevent ( pop_i, 
                //start_time, 
                //end_time, 
                opportunity_y, event_code, end_base){ 
                this->set_start_base (start_base); 
                recomb_opp += this->opportunity();
                };
    
    double opportunity() const { return this->opportunity_between ( this->start_base_, this->end_base_); };
    double opportunity_between( double start, double end ) const {
        assert ( end <=  this->end_base_);
        assert ( start >= this->start_base_);
        return this->opportunity_y_* ( end - start); };
    bool print_event();
    double start_base() const { return this->start_base_; }
    void set_start_base ( double base ) { this->start_base_ = base; }

    double start_base_;
};

/*!
 * \brief Derived class of Coalevent, recording the number and the time intervals of Migration events between two ForestState 
 */    
class Migrevent : public Coalevent{
    friend class CountModel;
    friend class ForestState;
    Migrevent( size_t pop_i, double opportunity_y, eventCode event_code, size_t mig_pop, double end_base )
             : Coalevent ( pop_i, 
                //start_time, 
                //end_time, 
                opportunity_y, event_code, end_base ){ 
                    this->set_mig_pop_i (mig_pop); 
                    };

    Migrevent(const Migrevent & previous_Coalevent) : Coalevent (previous_Coalevent) { 
        this->set_mig_pop_i ( previous_Coalevent.mig_pop() ); 
        };
    ~Migrevent(){};
    
    double opportunity() const {return this->opportunity_y_; };
    bool print_event();
    size_t mig_pop() const {return this->mig_pop_; } 
    void set_mig_pop_i( size_t i ){ this->mig_pop_ = i; }    
    
    size_t mig_pop_;
};

#endif
