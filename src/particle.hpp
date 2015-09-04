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
#include "coalevent.hpp"
#include <deque>
#include "segdata.hpp"

#ifndef NDEBUG
#define ForestStatedout (std::cout << "    ForestState ")
#else
#pragma GCC diagnostic ignored "-Wunused-value"
#define ForestStatedout 0 && (std::cout << "    ForestState ")
#endif

#pragma GCC diagnostic ignored "-Wsign-compare"

#ifndef PARTICLE
#define PARTICLE

extern int new_forest_counter;
extern int delete_forest_counter;
extern int recombination_counter; //DEBUG

struct TmrcaState {
    TmrcaState (double base, double tmrca) {
        this->base = base;
        this->tmrca = tmrca;
        }
    ~TmrcaState(){};
    double base;
    double tmrca;
    };

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
        ForestState(Model* model, RandomGenerator* random_generator, const vector<int>& record_event_in_epoch, bool own_model_and_random_generator);    /*!< \brief ForestState constructer, used when initialize ForestState from absolutely the first time */
        // Constructor, called when resampling, making new copy of a particle
        ForestState(const ForestState &current_state); /*!< \brief ForestState constructer, used when copy particle from a given particle */
        // Destructors //
        ~ForestState();
        void clear_eventContainer();

        // Resampling //
        void init_EventContainers( Model * model );
        void copyEventContainers(const ForestState & copied_state );
        void resample_recombination_position(void);

        // Update weight
        void include_haplotypes_at_tips(vector <int> &haplotypes_at_tips); /*!< \brief Update data to the particle */
        double calculate_likelihood( ); /*!< \brief Calculate the likelihood of the genealogy */
        valarray<double> cal_partial_likelihood_infinite(Node * node); /*!< Calculate the marginal likelihood of each node */
        double trackLocalTreeBranchLength();
        BranchLengthData trackSubtreeBranchLength ( Node * currentNode );

        // Extend
        double extend_ARG ( double mutation_rate, double extend_to, Segment_State segment_state, bool updateWeight=true, bool recordEvents=true );

        // Record events
        void record_Recombevent_b4_extension ( );
        void record_Recombevent_atNewGenealogy ( double event_height );
        void record_all_event(TimeInterval const &ti, double & recomb_opp_x_within_scrm);
        void clear_recomb_opp_within_scrm(){ this->recomb_opp_x_within_scrm = 0 ; }

        // Setters and getters:
        void setSiteWhereWeightWasUpdated( double site ){ this->site_where_weight_was_updated_=site; }
        double site_where_weight_was_updated() const { return site_where_weight_was_updated_; }
        void setParticleWeight(double weight) { this->particle_weight_ = weight; }
        double weight() const { return this->particle_weight_; }

        // Members
        vector < EvolutionaryEvent* > eventTrees;

        double site_where_weight_was_updated_;
        double particle_weight_;
        vector < TmrcaState > TmrcaHistory;
        const vector < int >& record_event_in_epoch;
        bool owning_model_and_random_generator;

        // Debugging tools
        std::string newick(Node *node) ;
};
#endif
