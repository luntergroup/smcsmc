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

#include"particle.hpp"
#include"pfparam.hpp"


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
        void ESS_resampling(valarray<double> weight_cum_sum, valarray<int> &sample_count, int mutation_at, double ESSthreshold, int num_state);        
        bool appendingStuffToFile(double x_end, PfParam &pfparam);
        void clean_old_states(double xstart);
        void clear();
        //int count_total_number_of_nodes();                

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
         * Methods
         */ 
        void extend_ARGs(double mutation_at,
                         double mutation_rate, 
                         bool withdata);
        void normalize_probability();    
        void update_state_weights_at_A_single_site(double mutation_at,
                                                   double mutation_rate, 
                                                   bool withdata,
                                                   vector <bool> haplotypes_at_tips);
        void push(ForestState * particle, double weight=1); /*!< If particle is new, initialize the weight as 1 */

        // Resampling
        void resample(valarray<int> & sample_count);
        void shifting(int number_of_particles);
        void trivial_resampling(size_t N, std::valarray<int> & sample_count);
        void systemetic_resampling(std::valarray<double> cum_sum, std::valarray<int>& sample_count, int sample_size);
        void update_cum_sum_array_find_ESS(std::valarray<double> & weight_cum_sum);
        
        
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


        /*!
         * Members
         */ 
        double ESS_;
        RandomGenerator* random_generator_;
        double current_printing_base_;
    };
