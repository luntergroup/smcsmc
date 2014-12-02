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
        //CountModel():Model(){};     
        CountModel(const Model& model, double lag = 0, double minimal_lag_update_ratio = 0.10 ) : Model( model ){ 
			this->const_lag_ = lag;
			this->const_minimal_lag_update_ratio_ = minimal_lag_update_ratio; 
		};
        ~CountModel(){};
        
        //
        // Methods
        //   
        void init();
        void extract_and_update_count( ParticleContainer &Endparticles , double current_base, bool end_data = false);
        void reset_model_parameters(double current_base, Model * model, bool online = true, bool force_update = false, bool print = true);
        void log_counts( PfParam& param );
        
        // DEBUG
        void print_recomb_count();        
    private:

        //
        // Methods
        //   
        // Initializeation
        void init_coal_and_recomb();
        void init_migr();
        void init_lags();
        
        // Reset parameters
        void reset_recomb_rate ( Model *model );
        void reset_Ne ( Model *model );
        void reset_mig_rate ( Model *model );
        //void reset_single_mig_rate ( Model *model );
        void initialize_mig_rate ( vector <vector<double>*> & rates_list );


        void update_coalescent_count( deque < Coalevent *> & CoaleventContainer_i, double weight, size_t x_end, vector<double>& total_coal_count, vector<double>& total_coal_opportunity ) ;
        void update_recombination_count( deque < Recombevent *> & RecombeventContainer_i, double weight, size_t x_start, size_t x_end, vector<double>& total_recomb_count, vector<double>& total_recomb_opportunity ) ;
        void update_migration_count( deque < Migrevent *> & MigreventContainer_i, double weight, size_t x_end, size_t epoch_idx );

        void compute_recomb_rate();
        void compute_mig_rate();

        void resize_Starevent ( deque < Coalevent *> & CoaleventContainer_i , int index) ;
        void resize_Starevent ( deque < Recombevent *> & RecombeventContainer_i , int index) ;
        void resize_Migrevent ( deque < Migrevent *> & MigreventContainer_i , int index) ;

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
        /*! The dimension of total_mig_count is number_of_epochs * number_of_population (from) * number_of_population (to).  For total_weighted_mig_opportunity only the 'from' population is important
         */         
        vector < vector < vector<double> > >  total_mig_count;
        vector < vector<double> >             total_weighted_mig_opportunity;

        vector < double > counted_to;
        vector < double > lags;
        vector < vector < vector<double> > > inferred_mig_rate; // This should be a 3-D vector 
        
        double recomb_count_;
        double recomb_opportunity_;
        double inferred_recomb_rate;

        double const_lag_;               
        double const_minimal_lag_update_ratio_; 
        double update_param_threshold_;
        double update_param_interval_;

        // DEBUG
        
        void print_pop_size();
        void print_change_time();
        void print_coal_count();
        bool print_mig_rate ( vector <vector<double>*> & rates_list );
        void check_model_updated_mig(Model * model);
        void check_model_updated_Ne(Model * model);    
    };    
#endif
