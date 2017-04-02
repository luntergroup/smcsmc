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
#include <deque>
#include <valarray>
#include "segdata.hpp"
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
extern int recombination_counter; //DEBUG
extern int recombination_event_called; //DEBUG

struct BranchLengthData {
    BranchLengthData (double partialBranchLength = 0, double subtreeBranchLength = -1)
    {
        this->partialBranchLength = partialBranchLength;
        this->subtreeBranchLength = subtreeBranchLength;
    }
    ~BranchLengthData(){};
    // this holds the branch length of the partial tree consisting of all leaf nodes
    // that carry data, up to the current node.
    double partialBranchLength;
    // this holds the branch length of the subtree subtending all leaf nodes that carry
    // data, but not including the branch from the subtree's root to the current node
    // Special case: if no descendants of the current node carry data, this is -1.
    double subtreeBranchLength;
};

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
    ForestState(Model* model, RandomGenerator* random_generator, const vector<int>& record_event_in_epoch,
                bool own_model_and_random_generator);    /*!< \brief Used to create a ForestState for the first time */
    ForestState(const ForestState &current_state);       /*!< \brief Copy constructor, used when resampling */

    // Destructors //
    ~ForestState();
    void clear_eventContainer();
    
    // Resampling //
    void init_EventContainers( Model * model );
    void copyEventContainers(const ForestState & copied_state );
    void resample_recombination_position(void);
    
    // Update weight
    void include_haplotypes_at_tips(vector <int> &haplotypes_at_tips); /*!< \brief Update data to the particle */
    double calculate_likelihood( bool ancestral_aware ); /*!< \brief Calculate the likelihood of the genealogy */
    valarray<double> cal_partial_likelihood_infinite(Node * node); /*!< Calculate the marginal likelihood of each node */
    double trackLocalTreeBranchLength();
    BranchLengthData trackSubtreeBranchLength ( Node * currentNode );
    
    // Extend
    double extend_ARG ( double extend_to );
    
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
    
    //// biased sampling    
    double importance_weight_over_segment( double previously_updated_to, double update_to );

    double find_delay( double coal_height );
    
    double sampleHeightOnWeightedBranch( Node* node, double length_left, double* bias_ratio) const;
    double sampleBiasedPoint_recursive( Node* node, double& length_left );
    
    double getWeightedLocalTreeLength() const;
    double getWeightedLengthBelow( Node* node ) const;
    double WeightedBranchLengthAbove( Node* node ) const;

    // What does this do?
    Node* trackLocalNode(Node *node) const;

    // below are overloaded
    virtual double sampleNextBase();
    virtual double samplePoint();
    
    // Debugging tools
    std::string newick(Node *node) ;
        
    // Members
    vector < EvolutionaryEvent* > eventTrees;
    double site_where_weight_was_updated_;
    double posterior_weight_;
    double pilot_weight_;
    double first_coalescence_height_;
    bool owning_model_and_random_generator;
    const vector < int >& record_event_in_epoch;
    std::priority_queue<DelayedFactor, std::vector<DelayedFactor>, CompareDFs > delayed_adjustments;
    double total_delayed_adjustment_;

    // DEBUG
    int recent_recombination_count;
    
};
#endif
