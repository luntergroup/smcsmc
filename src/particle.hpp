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
#include<iterator>
#include<valarray>
#include"vcf.hpp"
#include"pfparam.hpp"

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
    public:

        /*!
         * Constructors and Destructors
         */      
        //ForestState():Forest(){ }
        ForestState(Model* model, 
                    RandomGenerator* random_generator);    /*!< ForestState constructer, used when initialize ForestState from absolutely the first time */    
        ForestState(ForestState *current_state); /*!< \brief ForestState constructer, used when copy particle from a given particle */
        ~ForestState();
    
        /*!
         * Methods
         */ 
        void init(double weight=1.0, 
                  double weight_updated_at_site=0.0, 
                  ForestState * previous_state = NULL); /*!< Initialize ForestState member particle_weight_ and site_where_weight_was_updated_ */

        void clear_CoaleventContainer();
        void clear_RecombeventContainer();
        void clear_MigreventContainer();
        
        void record_all_event(TimeInterval const &ti);
        void record_Coalevent(size_t pop_i,
                          double start_time, 
                          double end_time, 
                          double opportunity, 
                          eventCode event_code);
        
        void record_Recombevent(size_t pop_i,
                          double start_time, 
                          double end_time, 
                          double opportunity, 
                          eventCode event_code);

        void record_Migrevent(size_t pop_i,
                          size_t mig_pop,
                          double start_time, 
                          double end_time, 
                          double opportunity, 
                          eventCode event_code);                          
        
        void include_haplotypes_at_tips(vector <bool> haplotypes_at_tips); /*!< Update data to the particle */        
        valarray<double> cal_marginal_likelihood(Node * node); /*!< Calculate the marginal likelihood of each node */
        double calculate_likelihood(bool withdata); /*!< Calculate the likelihood of the genealogy */
        
        /*!
         *  Setters and getters:
         */ 
        void setSiteWhereWeightWasUpdated(double site){this->site_where_weight_was_updated_=site;}        
        double site_where_weight_was_updated(){return site_where_weight_was_updated_;}

        void setParticleWeight(double weight){this->particle_weight_ = weight;};
        double weight(){return particle_weight_;};        

        /*!
         * Debugging tools
         */ 
        bool print_Coalevent();
        bool print_Recombevent();
        //bool print_Coalevent_out();

        /*!
         * Members
         */ 
        ForestState *previous_state; /*!< Pointer to the previous particle, i.e. the previous genealogy */
        int pointer_counter; /*!< Pointer counter. Record the number of particles that are pointing towards to THIS PARTICLE. */   
        
        vector < Coalevent*  > CoaleventContainer;   /*!< Coalescent events recorder */
        vector < Recombevent*> RecombeventContainer; /*!< Recombination events recorder */
        vector < Migrevent*  > MigreventContainer;   /*!< Migration events recorder */
        
    private:

        /*!
         * Members
         */ 
        double site_where_weight_was_updated_;
        double particle_weight_;
};


class ParticleContainer{
    public:

        /*!
         * Constructors and Destructors
         */      
        ParticleContainer();   
        ParticleContainer(Model* model, 
                  RandomGenerator* rg, 
                  size_t Num_of_states, 
                  vector <bool> data_for_init_states,  
                  bool withdata,
                  double initial_position);// this is used to create the particle initial states
        ~ParticleContainer(); //Destuctor needs to free memory of each particles that the pointers are pointing to...

        /*!
         * Methods
         */ 
        void update_state_to_data(Vcf * VCFfile, 
                                  Model * model, 
                                  valarray<double> & weight_cum_sum);
        void extend_ARGs(double mutation_at,
                         double mutation_rate, 
                         bool withdata);
        void normalize_probability();    
        void update_state_weights_at_A_single_site(double mutation_at,
                                                   double mutation_rate, 
                                                   bool withdata,
                                                   vector <bool> haplotypes_at_tips);
        void push(ForestState * particle, double weight=1); /*!< If particle is new, initialize the weight as 1 */
        void update_cum_sum_array_find_ESS(std::valarray<double> & weight_cum_sum);
        void ESS_resampling(valarray<double> weight_cum_sum, valarray<int> &sample_count, int mutation_at, double ESSthreshold, int num_state);
        void resample(valarray<int> & sample_count);
        void shifting(int number_of_particles);
        void trivial_resampling(size_t N, std::valarray<int> & sample_count);
        void systemetic_resampling(std::valarray<double> cum_sum, std::valarray<int>& sample_count, int sample_size);
        
        int count_total_number_of_nodes();        
        void clear();
        void clean_old_states(double xstart);
        
        bool appendingStuffToFile(double x_end, pfARG::param pfparam);
        
        /*!
         *  Setters and getters:
         */ 
        double ESS() const {return this->ESS_;};
        void set_ESS(double ess){this->ESS_ = ess;};
        
        RandomGenerator* random_generator() const { return this->random_generator_; }
        void set_random_generator(RandomGenerator *rg) {
                this->random_generator_ = rg; }
        
        double current_printing_base() const { return this->current_printing_base_;}
        void set_current_printing_base (double base) { this->current_printing_base_ = base;}
        
        //double max_weight() const {return this->max_weight_;}
        //void set_max_weight(double weight) {this->max_weight_ = weight;}

        /*!
         * Debugging tools
         */ 
        void print();
        bool check_state_orders();

        /*!
         * Members
         */ 
        
        vector <ForestState*> particles;
        
    private:

        /*!
         * Members
         */ 
        double ESS_;
        RandomGenerator* random_generator_;
        double current_printing_base_;
        //double max_weight_;
};
#endif
