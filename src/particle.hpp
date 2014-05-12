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


#include"scrm/forest.h"
#include"coalevent.hpp"

#ifndef NDEBUG
#define dout std::cout
#else
#pragma GCC diagnostic ignored "-Wunused-value"
#define dout 0 && std::cout
#endif
#pragma GCC diagnostic ignored "-Wsign-compare"

#ifndef PARTICLE
#define PARTICLE
extern int new_forest_counter;
extern int delete_forest_counter;


/*! 
 * \brief Derived class from Forest.
 * A particle is a special case of a Forest, with special members
 */
class ForestState : public Forest{      
    friend class ParticleContainer;
    friend class CountModel;
    
    // All members and methods are private
    private:
        //
        // Constructors and Destructors
        //        
        //ForestState():Forest(){ }
        ForestState(Model* model, 
                    RandomGenerator* random_generator);    /*!< \brief ForestState constructer, used when initialize ForestState from absolutely the first time */    
        ForestState(ForestState *current_state); /*!< \brief ForestState constructer, used when copy particle from a given particle */
        ~ForestState();
    
        //
        // Methods
        //           
        void init(double weight=1.0, 
                  double weight_updated_at_site=0.0, 
                  ForestState * previous_state = NULL); /*!< Initialize ForestState member particle_weight_ and site_where_weight_was_updated_ */

        // Update weight
        void include_haplotypes_at_tips(vector <bool> haplotypes_at_tips); /*!< \brief Update data to the particle */        
        double calculate_likelihood(bool withdata, bool finite = false); /*!< \brief Calculate the likelihood of the genealogy */
        valarray<double> cal_marginal_likelihood_finite(Node * node); /*!< Calculate the marginal likelihood of each node */
        valarray<double> cal_marginal_likelihood_infinite(Node * node); /*!< Calculate the marginal likelihood of each node */
        
        // Record events
        void record_all_event(TimeInterval const &ti);
        void record_Coalevent(size_t pop_i,
                          //double start_time, 
                          //double end_time, 
                          double opportunity, 
                          eventCode event_code);
        
        void record_Recombevent(size_t pop_i,
                          //double start_time, 
                          //double end_time, 
                          double opportunity, 
                          eventCode event_code);

        void record_Migrevent(size_t pop_i,
                          size_t mig_pop,
                          //double start_time, 
                          //double end_time, 
                          double opportunity, 
                          eventCode event_code);                          

        void clear_CoaleventContainer();
        void clear_RecombeventContainer();
        void clear_MigreventContainer();
        
        //
        // Setters and getters:
        //
        void setSiteWhereWeightWasUpdated(double site){this->site_where_weight_was_updated_=site;}        
        double site_where_weight_was_updated(){return site_where_weight_was_updated_;}
        void setParticleWeight(double weight){this->particle_weight_ = weight;};
        double weight(){return this->particle_weight_;};                
        void setAncestor ( size_t ancestor ) { this->ancestor_ = ancestor; }
        size_t ancestor(){ return this->ancestor_;}
        
        //
        // Debugging tools
        //
        bool print_Coalevent();
        bool print_Recombevent();
        bool print_Migrevent();

        //
        // Members
        //   
        ForestState *previous_state;                 /*!< \brief Pointer to the previous particle, i.e. the previous genealogy */
        int pointer_counter;                         /*!< \brief Pointer counter. Record the number of particles that are pointing towards to THIS PARTICLE. */   
        vector < Coalevent*  > CoaleventContainer;   /*!< \brief Coalescent events recorder */
        vector < Recombevent*> RecombeventContainer; /*!< \brief Recombination events recorder */
        vector < Migrevent*  > MigreventContainer;   /*!< \brief Migration events recorder */
        double site_where_weight_was_updated_;
        double particle_weight_;
        size_t ancestor_;
    };
#endif
