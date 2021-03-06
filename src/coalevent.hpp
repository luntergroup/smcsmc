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


#include "scrm/src/node.h"
#include "scrm/src/forest.h"
#include "arena.hpp"
#include "general.hpp"
#include "descendants.hpp"


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
    explicit EvolutionaryEvent( double start_height, size_t start_height_epoch, double end_height,
                                size_t end_height_epoch, double start_base, double end_base, int weight,
                                EvolutionaryEvent* parent = NULL ) :
        start_height(start_height),
        end_height(end_height),
        start_base_(start_base),
        end_base_(end_base),
        a( 0 ),
        parent_(parent),
        weight(weight),
        descendants( NO_DESCENDANTS )
    {
        if (parent_) parent_->increase_refcount();
        assert((start_height <= end_height) && (start_base <= end_base) && (start_base >= 0) );
        this->init();
    };
    // Constructor for migration/coalescence opportunity
    explicit EvolutionaryEvent( double start_height, size_t start_height_epoch, double end_height, size_t end_height_epoch, double end_base, size_t population_index, int weight, EvolutionaryEvent* parent = NULL ) :
        start_height(start_height),
        end_height(end_height),
        start_base_(-1),
        end_base_(end_base),
        a( population_index ),
        parent_(parent),
        weight(weight),
        descendants( NO_DESCENDANTS )
    {
        if (parent_) parent_->increase_refcount();
        assert(start_height <= end_height);
        this->init();
    };
    // Copy constructor
    EvolutionaryEvent( const EvolutionaryEvent& obj ) :
        start_height(obj.start_height),
        end_height(obj.end_height),
        start_base_(obj.start_base_),
        end_base_(obj.end_base_),
        a( obj.a.coal_migr_population ),
        parent_(obj.parent_),
        weight(obj.weight),
        event_data_(obj.event_data_),
        descendants( obj.descendants )
    {
        if(parent_) parent_->increase_refcount();
        this->init();
        event_data_ = obj.event_data_;
        children_updated_ = obj.children_updated_;
    }
    // Copy constructor that does not steal the parent pointer -- used for building tree modifying events
    EvolutionaryEvent( const int _dummy, const EvolutionaryEvent& obj ) :
        start_height(obj.start_height),
        end_height(obj.end_height),
        start_base_(obj.start_base_),
        end_base_(obj.end_base_),
        a( obj.a.coal_migr_population ),
        parent_(NULL),
        weight(obj.weight),
        event_data_(obj.event_data_),
        descendants( obj.descendants )
    {
        this->init();
        event_data_ = obj.event_data_;
        children_updated_ = obj.children_updated_;
    }
    // Destructor.  Should not be used when Arena and placement new is used for allocation
    ~EvolutionaryEvent() {
        std::cout << "Destructor is called -- problem!" << std::endl;
        if (parent_ && parent_->decrease_refcount_is_zero()) delete parent_;
    }
    // Delete function for use with Arena
    void deletethis( size_t epoch_idx ) {
        if (parent_ && parent_->decrease_refcount_is_zero()) parent_->deletethis( epoch_idx );
        // don't call destructor -- no need to, since class only contains simple types
        Arena::deallocate( (void*)this, epoch_idx );
    }

    // Methods
    bool is_recomb() const            { return start_base_ >= 0; }
    bool is_coalmigr() const          { return start_base_ < 0; }
    bool is_no_event() const          { return event_data_ == -2; }
    bool is_recomb_event() const      { assert(is_recomb()); return event_data_ > -2; }
    bool is_time_recomb_event() const { assert(is_recomb()); return event_data_ == 0; }
    bool is_coal_event() const        { assert(is_coalmigr()); return event_data_ == -1; }
    bool is_migr_event() const        { assert(is_coalmigr()); return event_data_ >= 0; }
    int recomb_event_count() const    { return is_recomb_event(); }
    int coal_event_count() const      { return is_coal_event(); }
    int migr_event_count() const      { return is_migr_event(); }
    double start_base() const         { return is_recomb() ? start_base_ : end_base_; }
    double end_base() const           { return end_base_; }
    void set_recomb_event_pos( double recomb_x_position ) {
        assert( this->is_recomb() );
        assert( this->is_no_event() );
        assert( recomb_x_position <= end_base_ );
        assert( start_base_ <= recomb_x_position );
        event_data_ = -1;
        a.recomb_pos = recomb_x_position; }
    void set_recomb_event_time( double recomb_t_position ) {
        assert( this->is_recomb() );
        assert( this->is_no_event() );
        assert( recomb_t_position <= end_height );
        assert( start_height <= recomb_t_position );
        event_data_ = 0;
        a.recomb_pos = recomb_t_position; }
    void set_coal_event() {
        assert( this->is_coalmigr() );
        assert( this->is_no_event() );
        event_data_ = -1; }
    void set_migr_event( int to_population ) {
        assert( this->is_coalmigr() );
        assert( this->is_no_event() );
        assert( to_population != a.coal_migr_population );
        assert( to_population >= 0 );
        event_data_ = to_population; }
    void repr( std::ostream& stream ) const {
        if (is_recomb()) {
            stream << "Rw" << std::setw(1) << weight << "[" << std::fixed << std::setprecision(2) << std::setw(10) 
                   << start_base_ << "," << std::setw(10) << end_base_ << "]x["
                   << std::setw(10) << start_height << "," << std::setw(10) << end_height << "]";
            if (!is_no_event()) {
                if (event_data_ == -1) {
                    stream << " x=" << a.recomb_pos << " *";
                } else {
                    stream << " t=" << a.recomb_pos << " *";
                }
            }
        } else if (is_coalmigr()) {
            stream << "T" << std::setw(1) << a.coal_migr_population << "w" << weight << "("
                   << std::fixed << std::setprecision(2) << std::setw(10) << end_base_ << ")x["
                   << std::setw(10) << start_height << "," << std::setw(10) << end_height << "]";
            if (is_coal_event())  stream << " c*";
            if (is_migr_event())  stream << " m*";
        } else {
            stream << "?? nothing ??";
        }
    }
    double coal_opportunity() const {
        assert (is_coalmigr());
        return weight * (end_height - start_height); }
    double coal_opportunity_between( double height0, double height1 ) const {
        assert (is_coalmigr());
        return weight * std::max(0.0, (std::min(height1,end_height) - std::max(height0,start_height))); }
    double migr_opportunity() const {
        assert (is_coalmigr());
        return end_height - start_height; }
    double migr_opportunity_between( double height0, double height1 ) const {
        assert (is_coalmigr());
        return std::max(0.0, (std::min(height1,end_height) - std::max(height0,start_height))); }
    double recomb_opportunity() const {
        assert (is_recomb());
        return weight * (end_height - start_height) * (end_base_ - start_base_); }
    double recomb_opportunity_between( double height0, double height1, double base0, double base1) const {
        return weight * std::max(0.0, std::min(height1,end_height) - std::max(height0,start_height))
                      * std::max(0.0, std::min(base1,end_base_) - std::max(base0,start_base_)); }
    bool recomb_event_overlaps_opportunity_t( double recomb_t_position ) const {
        return ( start_height <= recomb_t_position && recomb_t_position <= end_height ); }
    bool recomb_event_overlaps_opportunity_x( double recomb_x_position ) const {
        return ( start_base_ <= recomb_x_position && recomb_x_position <= end_base_ ); }
    int recomb_event_count_between( double height0, double height1, double base0, double base1, bool isEndOfSeq ) const {
        if (!is_recomb_event()) return 0;
        double base_ = event_data_ == -1 ? a.recomb_pos : start_base_;
        double height = event_data_ == 0 ? a.recomb_pos : end_height;
        return (base0 <= base_) && ( (base_ < base1) || isEndOfSeq ) && (height0 <= height) && (height < height1); }
    double recomb_event_base() const {
        if (!is_recomb_event()) return -1;
        // recombination events are recorded in the sequence segment following the event (arbitrary choice).  Therefore,
        // the sequence position of the recombination event is the start of the current event record's sequence segment
        // (for events that occur when sweeping in the x direction)
        return event_data_ == -1 ? a.recomb_pos : start_base_; }
    double recomb_event_time() const {
        if (!is_recomb_event()) return -1;
        return event_data_ == 0 ? a.recomb_pos : end_height; }
    void set_descendants( const Descendants_t descendants_ ) {
        descendants = descendants_; }
    Descendants_t get_descendants() const {
        return this->descendants; }
    double get_start_base() const {
        return start_base_; }
    double get_end_base() const {
        return end_base_; }
    double get_start_height() const {
        return start_height; }
    double get_end_height() const {
        return end_height; }
    double get_recomb_pos() const {
        return a.recomb_pos; }
    int get_weight() const {
        return weight;
    }
    size_t get_population() const {
        assert (is_coalmigr());
        return a.coal_migr_population; }
    size_t get_migr_to_population() const {
        assert (is_coalmigr());
        assert (is_migr_event());
        return (size_t)event_data_; }
    bool decrease_refcount_is_zero() {
        assert( ref_counter_ > 0 );
        ref_counter_--;
        return (ref_counter_ == 0); } // Return TRUE, if there is nothing pointing here.
    void increase_refcount() { ref_counter_++; }
    bool print_event() const;
    bool append_event( const EvolutionaryEvent& e );
    EvolutionaryEvent* parent() const { return parent_; }
    EvolutionaryEvent*& parent() { return parent_; }

    // methods for the update algorithm

    /* check whether this event is marked as 'removed' (and should subsequently be spliced out of the tree whenever it's encountered) */
    bool is_removed() const { return children_updated_ < 0; }
    /* marks as deleted; actual deletion is done by calling purge_events */
    void mark_as_removed() {
        assert (children_updated_ <= 0); // ensure events are not removed in mid-flow
        children_updated_ = -1; }
    /* Adds (newly made, refcount==1) this to tree */
    void add_leaf_to_tree( EvolutionaryEvent** eventptr_location, bool always_treat_as_tree=false ) {
        if (parent_ || always_treat_as_tree) {
            // this is a tree; replace existing tree off eventptr_location with this
            bool eventIsLost = (*eventptr_location)->decrease_refcount_is_zero();
            // decrease_refcount_is_zero returns TRUE, when ref count is zero, in which case, eventIsLost is TRUE
            // decrease_refcount_is_zero returns FALSE, when ref count is not zero, in which case, eventIsLost is FALSE

            // eventIsLost should always be false! which means ref count should never be zero!
            // assertion fails, when eventIsLost is TRUE, i.e. the ref_count is zero.
            assert (!eventIsLost);  // *eventptr_location should be referenced elsewhere
            _unused(eventIsLost);
        } else {
            // this is a single node; splice into existing tree
            parent_ = *eventptr_location;
        }
        *eventptr_location = this; }
    /* update posterior, and checks whether the posterior is now complete */
    bool update_posterior_is_done( double additional_posterior ) {
        assert (!is_removed());
        posterior_ += additional_posterior;
        children_updated_ = (children_updated_ + 1) % ref_counter_;
        return (children_updated_ == 0); }
    double get_and_reset_posterior() {
        assert (!is_removed());
        assert (children_updated_ == 0);
        double p = posterior_;
        posterior_ = 0.0;
        return p;
    }

    // friend functions, for managing the tree
    friend bool remove_event( EvolutionaryEvent** eventptr_location, size_t epoch_idx );               /* Removes event after we're done updating the relevant counters */
    friend EvolutionaryEvent* purge_events( EvolutionaryEvent** eventptr_location, size_t epoch_idx ); /* Purges previously removed events, and returns first active parent (or NULL) */
    friend class ForestState;                                                                          /* for void ForestState::resample_recombination_position(void) */
    friend class ParticleContainer;  // for debugging

    // Members
private:
    double start_height;
    double end_height;
    double start_base_;    // Recombinations: determines (w/end_base) the x-extent of recomb. opportunity.  For coal/migr, <0
    double end_base_;
    union A {
      size_t coal_migr_population; // index of population of branch that is migrating or coalescing
      double recomb_pos;           // position of recombination event, along right-hand (end_base) or top (end_height) edge
      A( size_t p ): coal_migr_population( p ) {}
      A(): coal_migr_population( 0 ) {}
    } a;

    // helper variables for the update algorithm (put here for packing / alignment)
    EvolutionaryEvent* parent_;
    double posterior_;
    int children_updated_;
    int ref_counter_;

    // remainder of core variables
    short weight;         // number of lineages (for recombination and coalescence, not migration) contributing to opportunity
    short event_data_ ;   // -2 == no event; otherwise type-specific meaning:
                          // recomb:    -1 == event at top edge (time-wise sampling)
                          //             0 == event at right-hand edge (sequence-wise sampling)
                          // coal/migr: -1 == coalescent event;
                          //             0..: migration to this population index
    Descendants_t descendants;     // descendants of lineage where recombination occurred
    
    void set_posterior ( const double posterior ){ this->posterior_ =  posterior; }
    double posterior () const { return this->posterior_; }
    void set_children_updated ( const int children ){ this->children_updated_ = children; }
    int children_updated() const { return this->children_updated_; }
    void set_event_data ( const short event_data ){ this->event_data_ = event_data; }
    short event_data() const { return this->event_data_; }
    void set_end_base( double end_base ) { end_base_ = end_base; }

    void set_ref_counter ( const int ref_counter ){ this->ref_counter_ = ref_counter; }
    int ref_counter() const { return this->ref_counter_; }

    void init(){
        this->set_posterior(0.0);
        this->set_children_updated( 0 );
        this->set_event_data ( (short)-2 );
        this->set_ref_counter( 1 );
    }
};


#endif
