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


/* An EvolutionaryEvent describes the opportunity for recombination, migration and coalescent events,
 * their location, and optionally one of these events.
 * 
 * It describes either: migration and coalescent opportunity, or recombination opportunity.
 * In the case of migration/coalescence, a (novel or existing) lineage in an existing genealogy is 
 * traversed backwards in time over an interval, at a particular sequence position.  The data structure
 * describes a subset of the interval, ending when
 *  (1) an event happens (migration or coalescence);
 *  (2) an existing coalescence or migration event is passed (changing the number of contemporary branches)
 * The data structure does straddle epochs, between which rate parameters may change.
 * 
 * In the case of recombination, the data structure describes a two-dimensional (sequence * time) area
 * of opportunity.  The algorithm identifies these areas in one of two ways:
 *  (1) by extending the coalescent manifold in the sequence direction, until a recombination is
 *      encountered on the local tree;
 *  (2) by 'filling in' an area between a non-local branch into which a new, local ranch has coalesced.
 *      This process evolves backwards in time, and is ended by a new recombination event, or encountering
 *      a coalescent event on the non-local branch (into either a local or non-local branch). 
 * In the first case, the data structure spans time intervals with constant number of contemporary branches.
 * 
 * The 'start_base' member is -1 if the data structure refers to coalescence/migration opportunity
 * The 'weight' member is equal to the number of contemporary branches (and weighs coalescence and 
 *   recombination opportunity, not migration opportunity.)
 * 
 */

class EvolutionaryEvent {
public:
	// Constructor for a recombination opportunity
	explicit EvolutionaryEvent( double start_height, double end_height, double start_base, double end_base, int weight ) :
	                   start_height(start_height),
	                   end_height(end_height),
	                   start_base_(start_base),
	                   end_base_(end_base),
	                   weight(weight),
	                   a( 0 ) {
	                      assert((start_height <= end_height) && (start_base <= end_base) && (start_base >= 0) ); 
	                   };
	// Constructor for migration/coalescence opportunity
	explicit EvolutionaryEvent( double start_height, double end_height, double end_base, size_t population_index, int weight ) :
			           start_height(start_height),
			           end_height(end_height),
			           start_base_(-1),
			           end_base_(end_base),
			           weight(weight),
			           a( population_index ) { 
                          assert(start_height <= end_height); 
                       };

	// Methods
	bool is_recomb() const       { return start_base_ >= 0; }
	bool is_coalmigr() const     { return start_base_ < 0; }
	bool is_no_event() const     { return event_data == -2; }
	bool is_recomb_event() const { assert(is_recomb()); return event_data > -2; }
	bool is_coal_event() const   { assert(is_coalmigr()); return event_data == -1; }
	bool is_migr_event() const   { assert(is_coalmigr()); return event_data >= 0; }
	int recomb_event_count() const { return is_recomb_event(); }
	int coal_event_count() const   { return is_coal_event(); }
	int migr_event_count() const   { return is_migr_event(); }
	double start_base() const      { return start_base_; }
	double end_base() const        { return end_base_; }
	void set_recomb_event_pos( double recomb_x_position ) {
		assert( this->is_recomb() );
		assert( this->is_no_event() );
		assert( recomb_x_position <= end_base_ );
		assert( start_base_ <= recomb_x_position );
		event_data = -1;
		a.recomb_pos = recomb_x_position; }
	void set_recomb_event_time( double recomb_t_position ) {
		assert( this->is_recomb() );
		assert( this->is_no_event() );
		assert( recomb_t_position <= end_height );
		assert( start_height <= recomb_t_position );
		event_data = 0;
		a.recomb_pos = recomb_t_position; }
	void set_coal_event() {
		assert( this->is_coalmigr() );
		assert( this->is_no_event() );
		event_data = -1; }
	void set_migr_event( int to_population ) {
		assert( this->is_coalmigr() );
		assert( this->is_no_event() );
		assert( to_population != a.coal_migr_population );
		assert( to_population >= 0 );
		event_data = to_population; }
	double coal_opportunity() const {
		assert (is_coalmigr());
		return weight * (end_height - start_height); }
	double migr_opportunity() const {
		assert (is_coalmigr());
		return end_height - start_height; }
	double recomb_opportunity() const {
		assert (is_recomb());
		return weight * (end_height - start_height) * (end_base_ - start_base_); }
	double recomb_opportunity_between( double height0, double height1, double base0, double base1) const {
		return weight * max(0.0, min(height1,end_height) - max(height0,start_height)) * max(0.0, min(base1,end_base_) - max(base0,start_base_)); }
	size_t get_population() const {
		assert (is_coalmigr());
		return a.coal_migr_population; }
	size_t get_migr_to_population() const {
		assert (is_coalmigr());
		assert (is_migr_event());
		return (size_t)event_data; }
	bool decrease_refcount_is_zero() {
		assert( ref_counter-- > 0 );
		return (ref_counter == 0); }
	void increase_refcount() { ref_counter++; }
	bool print_event();
        bool append_event( const EvolutionaryEvent& e );
	              
	// Members
private:
public: // for CountModel::update_recombination_count; temporarily
	double start_height;
	double end_height;
private:
    double start_base_;    // Recombinations: determines (w/end_base) the x-extent of recomb. opportunity.  For coal/migr, <0
	double end_base_;
    int event_data {-2};   // -2 == no event; otherwise type-specific meaning:
	                       // recomb:    -1 == event at top edge (time-wise sampling)
	                       //             0 == event at right-hand edge (sequence-wise sampling)
	                       // coal/migr: -1 == coalescent event;
	                       //             0..: migration to this population index
	int weight;            // number of lineages (for recombination and coalescence, not migration) contributing to opportunity
    int ref_counter {1};
	union A {
	  size_t coal_migr_population; // index of population of branch that is migrating or coalescing
	  double recomb_pos;           // position of recombination event, along right-hand (end_base) or top (end_height) edge
	  A( int p ): coal_migr_population( p ) {}
	  A(): coal_migr_population( 0 ) {}
	} a;
};


#endif
