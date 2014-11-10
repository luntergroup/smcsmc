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
        ForestState(Model* model, RandomGenerator* random_generator);    /*!< \brief ForestState constructer, used when initialize ForestState from absolutely the first time */    
        // Constructor, called when resampling, making new copy of a particle
        ForestState(const ForestState &current_state); /*!< \brief ForestState constructer, used when copy particle from a given particle */
        // Destructors //
        ~ForestState();
        void clear_CoaleventContainer();
        void clear_RecombeventContainer();
        void clear_MigreventContainer();

        // Resampling //
        void init_EventContainers( Model * model );
        void copyEventContainers(const ForestState & copied_state );

        void making_copies( int number_of_copies );

        // Update weight
        void include_haplotypes_at_tips(vector <int> &haplotypes_at_tips); /*!< \brief Update data to the particle */        
        double calculate_likelihood( ); /*!< \brief Calculate the likelihood of the genealogy */
        valarray<double> cal_marginal_likelihood_infinite(Node * node); /*!< Calculate the marginal likelihood of each node */
        
        
        vector <double> opportunity_y_s ; 

        // Record events
        void compute_opportunity_y_s ();
        void record_Recombevent_b4_extension ( );
        void record_Recombevent_atNewGenealogy ( double event_height );
        //void record_the_final_recomb_opportunity ( double loci_length );
        void record_all_event(TimeInterval const &ti, double & recomb_opp_x_within_scrm);
        void record_Coalevent(size_t pop_i,
                          //double start_time, 
                          //double end_time, 
                          double opportunity, 
                          eventCode event_code, double end_base);
                                  
        void record_Recombevent(size_t pop_i,
                          //double start_time, 
                          //double end_time, 
                          double opportunity, 
                          eventCode event_code, double start_base, double end_base);

        void record_Migrevent(size_t pop_i,                          
                          //double start_time, 
                          //double end_time, 
                          double opportunity, 
                          eventCode event_code, size_t mig_pop, double end_base);                          
        
        //void record_recomb_opp_within_scrm ( double recomb_rate) const {
            //cout << "recomb_rate " << recomb_rate <<endl;
            //recomb_opp_x_within_scrm = recomb_opp_x_within_scrm + recomb_rate / this->model().recombination_rate();
        //}
        void clear_recomb_opp_within_scrm(){ this->recomb_opp_x_within_scrm = 0 ; }
        
        // Setters and getters: //
        void setSiteWhereWeightWasUpdated( double site ){ this->site_where_weight_was_updated_=site; }
        double site_where_weight_was_updated() const { return site_where_weight_was_updated_; }
        void setParticleWeight(double weight) { this->particle_weight_ = weight; }
        double weight() const { return this->particle_weight_; }
        void setAncestor ( size_t ancestor ){ this->ancestor_ = ancestor; }
        size_t ancestor() const { return this->ancestor_; }

        // Members //
        vector < deque < Coalevent* > > CoaleventContainer;   /*!< \brief Coalescent events recorder */
        vector < deque < Recombevent* > > RecombeventContainer; /*!< \brief Recombination events recorder */
        vector < deque < Migrevent* > > MigreventContainer;   /*!< \brief Migration events recorder */
                
        double site_where_weight_was_updated_;
        double particle_weight_;
        size_t ancestor_;
        vector < TmrcaState > TmrcaHistory;
        
        vector < ForestState* > ForestState_copies; // NEW        

        // Debugging tools //
        bool print_Coalevent();
        bool print_Recombevent();
        bool print_Migrevent();
        
        
        //valarray<double> cal_marginal_likelihood_finite(Node * node); /*!< Calculate the marginal likelihood of each node */
    };
#endif
