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


#include "particleContainer.hpp"


#ifndef COUNT
#define COUNT


/*! \brief Derived class of Model, used for inference.
 */         
class CountModel: public Model{

    public:    
        //
        // Constructors and Destructors
        //    
        CountModel():Model(){};     
        CountModel(const Model& model, double lag = 0 ):Model( model ){ this->const_lag_ = lag; };
        ~CountModel(){};
        
        //
        // Methods
        //   
        void init();
        double extract_and_update_count( ParticleContainer &Endparticles , double current_base, bool end_data = false);
        void reset_model_Ne( Model * model, bool online = true, bool print = true);
        void count_events_in_one_interval( ParticleContainer &Endparticles, size_t time_i, size_t pop_j, double x_start, double x_end);
        
        // DEBUG
        void print_recomb_counts();
        void print_pop_size();
        void print_change_time();
        void print_Time_count_pop();

        //
        // Members
        //   
        vector < double > previous_base;
        vector < double > lags;
        vector < vector<double> > inferred_mig_rate; // This should be a 3-D vector 
        
        double inferred_recomb_rate;
        
    private:

        //
        // Methods
        //   
        void init_coal_and_recomb();
        void init_migr();
        void init_lags();
        void compute_recomb_rate();
        void compute_mig_rate();
        void check_model_updated_Ne(Model * model);

        void count_events_in_one_interval_alt( ParticleContainer &Endparticles, size_t time_i, size_t pop_j, double x_start, double x_end);

        //void check_CountModel_Ne();

        //
        // Members
        //   
        /*! The dimension of total_coal_count, total_weighted_coal_opportunity, total_recomb_count, total_weighted_recomb_opportunity are number_of_time_interval * number_of_population
         */ 
        vector < vector<double> >   total_coal_count;
        vector < vector<double> >   total_weighted_coal_opportunity;        
        vector < vector<double> >   total_recomb_count;
        vector < vector<double> >   total_weighted_recomb_opportunity;

        /*! The dimension of total_mig_count, total_weighted_mig_opportunity are number_of_population * number_of_population
         */         
        vector < vector<double> >   total_mig_count;
        vector < vector<double> >   total_weighted_mig_opportunity;
        
        double recomb_count_;
        double const_lag_;
    };    
#endif
