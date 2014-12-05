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

#ifndef BIG_TO_SMALL_RATIO
#define BIG_TO_SMALL_RATIO 256
//#define BIG_TO_SMALL_RATIO 1
#endif

#ifndef CUM_TO_BIG_RATIO
#define CUM_TO_BIG_RATIO 5000
//#define CUM_TO_BIG_RATIO 1
#endif

class Two_doubles {
	public:
		Two_doubles( double init = 0 ){
			big_  = 0;
			small_ = init;
			//big_added_counter_ = 0;
		};
		~Two_doubles(){};
		
		void add_to_small ( double added ){
			//cout.precision(15);
			//cout << "adding "<< added << " to small " << small_ <<" wheas big is "<< big_ <<endl;
			small_ += added;
			//assert ( big_ > small_ );
			if ( small_ * BIG_TO_SMALL_RATIO > big_ && big_ != 0 ){
				add_small_to_big();
			}
			assert ( big_ > small_ );
		}
		
		void add_to_big ( double added ) {
			//cout.precision(15);
			//cout << "adding "<< added << " to big " << big_ <<" wheas small is "<< small_ <<endl;
			big_ += added;
			//big_added_counter_++;
			assert ( big_ > small_ );
			if ( big_ > small_ * BIG_TO_SMALL_RATIO * CUM_TO_BIG_RATIO ){
				//cout << " ok, bump up big!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
				add_big_to_cumsum( );
				assert ( cumsum_ > big_ );
			}			
		}
		
		void add ( double added ){
			if ( big_ == 0 ){
				if ( added > BIG_TO_SMALL_RATIO * small_ ){
					add_to_big ( added );
					return;
				}
				
				if ( added * BIG_TO_SMALL_RATIO < small_ ){
					switch_big_and_small();					
				}
				
				this->add_to_small ( added );	
				return;
			}
			else if ( added * BIG_TO_SMALL_RATIO < big_ ){
				this->add_to_small ( added );
				return;
			}
			else{
				add_to_big ( added );
				return;
			}
		}

		void compute_final_answer(){
			//cout.precision(15);
			//cout << "adding small "<< small_ << " to big " << big_ <<endl;
			this->add_small_to_big( );
			//cout << "adding big "<< big_ << " to cumsum " << cumsum_ <<endl;
			this->add_big_to_cumsum();	
		}
		
		double final_answer () {
			return this->cumsum_;
		}
		
		void add_small_to_big( ){
			//cout << "adding small "<< small_ << " to big " << big_ <<endl;
			add_to_big ( small_ );
			small_ = 0;
		}
		
		void add_big_to_cumsum () {
			//cout.precision(15);
			//cout << "adding big "<< big_ << " to cumsum " << cumsum_ <<endl;
			//big_added_counter_ = 0;
			cumsum_ += big_;
			big_ = 0;
		}
		
		void switch_big_and_small ( ){
			cout << " switching big and small !!!!!!!!!!!!!!!!!!!!!!"<<endl;
			double tmp = big_;
			big_ = small_;
			small_ = tmp;
		}
	
	private:
		double cumsum_;
		double big_, small_;
		//size_t big_added_counter_;
};

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


        void update_coalescent_count( deque < Coalevent *> & CoaleventContainer_i, double weight, double x_end, vector<Two_doubles>& total_coal_count, vector<Two_doubles>& total_coal_opportunity ) ;
        void update_recombination_count( deque < Recombevent *> & RecombeventContainer_i, double weight, double x_start, double x_end, vector<Two_doubles>& total_recomb_count, vector<Two_doubles>& total_recomb_opportunity ) ;
        void update_migration_count( deque < Migrevent *> & MigreventContainer_i, double weight, double x_end, size_t epoch_idx );

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
        vector < vector<Two_doubles> >   total_coal_count;
        vector < vector<Two_doubles> >   total_weighted_coal_opportunity;        
        vector < vector<Two_doubles> >   total_recomb_count;
        vector < vector<Two_doubles> >   total_weighted_recomb_opportunity;
        /*! The dimension of total_mig_count is number_of_epochs * number_of_population (from) * number_of_population (to).  For total_weighted_mig_opportunity only the 'from' population is important
         */         
        vector < vector < vector<Two_doubles> > >  total_mig_count;
        vector < vector<Two_doubles> >             total_weighted_mig_opportunity;

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
