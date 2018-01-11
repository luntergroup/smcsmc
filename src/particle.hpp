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


#include "forest.h"
#include "coalevent.hpp"
#include "pfparam.hpp"
#include "segdata.hpp"

#include <deque>
#include <valarray>
#include <queue>
#include <cmath>

#ifndef NDEBUG
#define ForestStatedout (std::cout << "    ForestState ")
#else
#pragma GCC diagnostic ignored "-Wunused-value"
#define ForestStatedout 0 && (std::cout << "    ForestState ")
#endif

#ifndef NDEBUG
#define newForestdout (std::cout << "    NEW FOREST ")
#else
#pragma GCC diagnostic ignored "-Wunused-value"
#define newForestdout 0 && (std::cout << "    NEW FOREST ")
#endif

#pragma GCC diagnostic ignored "-Wsign-compare"

#ifndef PARTICLE
#define PARTICLE

extern int new_forest_counter;
extern int delete_forest_counter;

// machinery for delayedIS

class DelayedFactor {
    
public:
    // single-activation instance
    DelayedFactor(double pos, double factor) : application_position(pos), importance_factor(factor), delta(0), k(1) {}

    // k-activation instance
    DelayedFactor(double final_pos, double factor, double cur_pos, int k) : k(k) {
        if (k < 1) throw std::runtime_error("Cannot have DelayedFactor with k < 1");
        if (k == 1) {
            delta = 0; // not used
            application_position = final_pos;
            importance_factor = factor;
        } else {
            delta = (final_pos-cur_pos) / ((1<<k)-1);
            application_position = cur_pos + delta;
            importance_factor = pow(factor, 1.0/k);
        }
    }
    // internal: next activation
    DelayedFactor( const DelayedFactor& df, bool _ ) : application_position( df.application_position + 2*df.delta ),
                                                       importance_factor( df.importance_factor ),
                                                       delta( 2*df.delta ),
                                                       k( df.k-1 ) {}

    double application_position;
    double importance_factor;
    // for piecewise constant (approximately continuous) application:
    double delta;
    int k;

    void print_info() const {
        std::clog << "\nThe application position of this DF is " << this->application_position << std::endl;
        std::clog << "The importance factor of this DF is " << this->importance_factor << std::endl;
    }
};

class CompareDFs {
public:
    bool operator() (const DelayedFactor &lhs, const DelayedFactor &rhs) const {
        return lhs.application_position > rhs.application_position;
    }
};


// quantiles of the terminal branch lengths
class TerminalBranchLengthQuantiles {
public:
    TerminalBranchLengthQuantiles() {}
    vector<double> quantiles;        // between 0 and 1, ordered
    vector<vector<double>> lengths;  // [lineage][quantile]
    double mean_total_branch_length; // used in the approximate likelihood
};


/*!
 * \brief Derived class from Forest.
 * A particle is a special case of a Forest, with special members
 */
class ForestState : public Forest{
#ifdef UNITTEST
    friend class TestForestState;
#endif
    friend class ParticleContainer;
    friend class CountModel;

    // All members and methods are private
private:
    // Constructor, called at initial stage //
    ForestState(Model* model, RandomGenerator* random_generator, PfParam* pfparam,
                bool own_model_and_random_generator);    /*!< \brief Used to create a ForestState for the first time */
    ForestState(const ForestState &current_state);       /*!< \brief Copy constructor, used when resampling */

    // Destructors //
    ~ForestState();
    void clear_eventContainer();

    // Spawn off new particle from this with multiplicity > 1, leaving this with multiplicity = 1
    ForestState* spawn() {
        assert( multiplicity() > 1 );
        ForestState* new_particle = new ForestState(*this);
        new_particle->setMultiplicity( multiplicity() - 1 );
        setMultiplicity(1);
        return new_particle;
    }
    
    // Resampling //
    void init_EventContainers( Model * model );
    void copyEventContainers(const ForestState & copied_state );
    void resample_recombination_position(void);
    
    // Update weight
    double calculate_likelihood( bool ancestral_aware, const vector<int>& haplotype ); /*!< \brief Calculate the likelihood of the genealogy */
    void cal_partial_likelihood_infinite(Node * node, const vector<int>& haplotype, double marginal_likelihood[2]); /*!< Calculate the marginal likelihood of each node */
    double trackLocalTreeBranchLength( const vector<int>& data_at_site );
    double trackSubtreeBranchLength ( Node * currentNode, const vector<int>& data_at_site );
    
    // Extend
    void extend_ARG ( double extend_to, int leaf_status, const vector<int>& data_at_site,
		      vector<ForestState*>* pParticleContainer = NULL);
    
    // Record events
    void record_recomb_extension ( );
    void record_recomb_event( double event_height );
    void record_all_event(TimeIntervalIterator const &ti);
    
    // Setters and getters:
    void setSiteWhereWeightWasUpdated( double site ){ this->site_where_weight_was_updated_=site; }
    double site_where_weight_was_updated() const { return site_where_weight_was_updated_; }
    void setPosteriorWeight(double weight) { assert (weight > 0.0); this->posterior_weight_ = weight; }
    void setPilotWeight(double weight) { this->pilot_weight_ = weight; }
    double posteriorWeight() const { return this->posterior_weight_; }
    double pilotWeight() const { return this->pilot_weight_; }
    int multiplicity() const { return multiplicity_; }
    void setMultiplicity( int multiplicity ) { multiplicity_ = multiplicity; }
    void set_max_epoch_to_record( int max_epoch_to_record ) { max_epoch_to_record_ = max_epoch_to_record; }
    bool tree_modified_after_last_mutation() const { 
        // this works because the likelihood cache is updated at every mutation, and reset to -1.0 after any change
        return (likelihood_cache < 0.0); 
    }

    // Save and restore recombination rate index
    void save_recomb_state() { _current_seq_idx = model().get_position_index(); }
    void restore_recomb_state() { writable_model()->resetSequencePosition( _current_seq_idx ); }

    // Auxiliary particle filter
    void includeLookaheadLikelihood( const Segment& segment, const TerminalBranchLengthQuantiles& terminal_branch_lengths );
    void removeLookaheadLikelihood() { pilot_weight_ /= lookahead_weight_; lookahead_weight_ = 1.0; }
    
    // Biased sampling
    void adjustWeights(double adjustment) {
        posterior_weight_ *= adjustment;
        pilot_weight_ *= adjustment;
    }
    void adjustWeightsWithDelay(double adjustment, double delay, int k=1) {
        posterior_weight_ *= adjustment;
        if ((adjustment > 0.99 && adjustment < 1.01) || (delay <= 1)) {
            // don't delay marginal adjustments
            pilot_weight_ *= adjustment;
        } else {
            total_delayed_adjustment_ *= adjustment;
            delayed_adjustments.push( DelayedFactor ( current_base() + delay, adjustment, current_base(), k ) );
        }
    }
    void applyDelayedAdjustment() {
        const DelayedFactor& df = delayed_adjustments.top();
        pilot_weight_ *= df.importance_factor;
        total_delayed_adjustment_ /= df.importance_factor;
        // if necessary, add new, updated adjustment
        if (df.k > 1) {
            delayed_adjustments.push( DelayedFactor( df, true ) );
        }
        delayed_adjustments.pop();
        assert( std::abs(posterior_weight_ - pilot_weight_ * total_delayed_adjustment_) <= .001 * posterior_weight_ );
    }
    
    // biased sampling    
    double find_delay( double coal_height );
    double importance_weight_over_segment( double previously_updated_to, double update_to );
    double sampleOrMeasureWeightedTree( const Node* node, double& length_left, double& local_weight, const bool recomb_bias, const bool recomb_guide );
    void inline accumulateBranchLengths( double rate, bool recomb_bias, const Node* node, double& length, double& weight );
    double getWeightedLocalTreeLength( const bool recomb_bias, const bool recomb_guide );
    Node* trackLocalNode(Node *node) const;

    // trial sampling
    double sampleNextGenealogyWithoutImplementing();
    void dontImplementFixedTimeEvent(TimeIntervalIterator &ti);
    void dontImplementNoEvent(const TimeInterval &ti, bool &coalescence_finished);
    void dontImplementRecombination(const Event &event, TimeIntervalIterator &ti);
    void dontImplementCoalescence(const Event &event, TimeIntervalIterator &tii);
    Node* virtualPossiblyMoveUpwards(Node* node, const TimeInterval &time_interval);

    // below are overloaded
    virtual double sampleNextBase( bool record_and_bias );
    virtual double samplePoint( bool record_and_bias );
    
    // Debugging tools
    std::string newick(Node *node) ;
        
    // Members
    vector < EvolutionaryEvent* > eventTrees;
    double site_where_weight_was_updated_;
    double posterior_weight_;
    double pilot_weight_;
    double lookahead_weight_;
    int    multiplicity_;
    int    _current_seq_idx;                         // stores variable model.h, so that each particle looks at correct recomb rate
    double first_event_height_;                      // NOTE: this is a temporary variable; move elsewhere?
    double first_coal_height_;
    double last_coal_height_;
    double total_local_branch_length_;               // NOTE: this is a temporary variable; move elsewhere?
    PfParam& pfparam;                                // to give access to record_event_in_epoch and recomb_bias.  NOTE: move elsewhere?
    double recombination_bias_importance_weight_;    // factor of importance weight due to recombination biasing; delay is treated specially for this IW
    std::priority_queue<DelayedFactor, std::vector<DelayedFactor>, CompareDFs > delayed_adjustments;
    bool owning_model_and_random_generator;
    int max_epoch_to_record_;
    double likelihood_cache;

    // DEBUG
    //int recent_recombination_count;
    double total_delayed_adjustment_;
    
};
#endif
