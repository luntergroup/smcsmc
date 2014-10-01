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


#ifndef PARTICLECONTAINER
#define PARTICLECONTAINER

/*! \brief ForestState Container, particle filters*/
class ParticleContainer{
    friend class CountModel;

    public:

        //
        // Constructors and Destructors
        //        
        ParticleContainer(Model* model, 
                  //size_t random_seed,
                  MersenneTwister *rg,
                  size_t Num_of_states, 
                  double initial_position,
                  bool heat_bool );// this is used to create the particle initial states
        ~ParticleContainer(); //Destuctor needs to free memory of each particles that the pointers are pointing to...

        //
        // Methods
        //   
        void update_state_to_data(VariantReader * VCFfile, 
                                  Model * model, 
                                  valarray<double> & weight_cum_sum);
                                  //bool finite_bool = false);
        void extend_ARGs(double mutation_at,
                             double mutation_rate, 
                             bool withdata);
        void set_particles_with_random_weight();
        void ESS_resampling(valarray<double> weight_cum_sum, valarray<int> &sample_count, int mutation_at, double ESSthreshold, int num_state);        
        bool appendingStuffToFile(double x_end, PfParam &pfparam);
        void cumulate_recomb_opportunity_at_seq_end( double seqend );
        void normalize_probability();    

        void clear();

        void print_particle_probabilities();

        //
        // Debugging tools
        //
        void print();
        bool check_state_orders();
        
    private:
        
        //
        // Methods
        //   
       
        void update_state_weights_at_A_single_site(double mutation_at,
                                                   double mutation_rate, 
                                                   bool withdata,
                                                   vector <int> &haplotypes_at_tips);
                                                   //bool finite_bool);
        void push(ForestState * particle, double weight=1); /*!< If particle is new, initialize the weight as 1 */        

        // Resampling
        void resample(valarray<int> & sample_count);
        void duplicate_particles ( valarray<int> & sample_count );
        
        void resample_for_check(valarray<int> & sample_count);
        
        void shifting(int number_of_particles);
        void trivial_resampling( std::valarray<int> & sample_count, size_t num_state );

        void systemetic_resampling(std::valarray<double> cum_sum, std::valarray<int>& sample_count, int sample_size);
        void update_cum_sum_array_find_ESS(std::valarray<double> & weight_cum_sum);
        
        
        
        
        
        
        
        //
        // Setters and getters:
        //
        double ESS() const {return this->ESS_;};
        void set_ESS(double ess){this->ESS_ = ess;};
        RandomGenerator* random_generator() const { return this->random_generator_; }
        //void set_random_generator(RandomGenerator *rg) { this->random_generator_ = rg; }        
        double current_printing_base() const { return this->current_printing_base_;}
        void set_current_printing_base (double base) { this->current_printing_base_ = base;}

        //
        // Members
        //
        vector <ForestState*> particles;
        double ESS_;
        RandomGenerator* random_generator_; // This is for particle filter only, 
        double current_printing_base_;
        bool heat_bool_ ;
    };

#endif
