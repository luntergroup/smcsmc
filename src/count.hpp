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


#include "particleContainer.hpp"
#include <ostream>

#ifndef COUNT
#define COUNT

#ifndef BIG_TO_SMALL_RATIO
#define BIG_TO_SMALL_RATIO 67108864 // 2^26
#endif

/*! \brief Derived class of Model, used for inference.
 */
class CountModel: public Model {

    public:
        //
        // Constructors and Destructors
        //
        CountModel(const Model& model, double lag = 0, double minimal_lag_update_ratio = 0.10 ) : Model( model ){
            this->const_lag_ = lag;
            this->const_minimal_lag_update_ratio_ = minimal_lag_update_ratio;
            this->sample_size_ = model.sample_size();
        };
        ~CountModel(){};

        //
        // Methods
        //
        void init();
        void extract_and_update_count( ParticleContainer &Endparticles , double current_base, bool end_data = false);
        void reset_model_parameters(double current_base, Model * model, bool useCap = false, double cap = 200000, bool online = true, bool force_update = false, bool print = true);
        void log_counts( PfParam& param );
        void reset_lag( std::vector<double> survival, double lag_fraction = 1 );
        void update_resample_count() { total_resample_count++; }

        // DEBUG
        void print_recomb_count();
        vector<double> check_lags() const {return lags;}
        void dump_local_recomb_logs( ostream& stream, double locus_length, int iteration );

    private:
        // Initialisation
        void init_coal_and_recomb();
        void init_migr();
        void init_lags();
        void init_local_recomb();

        // Reset parameters
        void reset_recomb_rate( Model *model );
        void reset_Ne( Model *model, bool useCap = false, double cap = 200000 );
        void reset_mig_rate( Model *model );
        void clear_2d_vector( vector <vector<double>> & rates_list );

        void update_all_counts( EvolutionaryEvent** event_ptr, double weight, vector<double>& update_to, size_t epoch_idx );
        void update_all_counts_single_evolevent( EvolutionaryEvent* event, double weight, vector<double>& update_to, size_t first_epoch_to_update );
        void record_local_recomb_events( double x_start, double x_end, double weight, double opportunity, double event_base, Descendants_t descendants );
        void update_delayed_weight_count( double opportunity ) { total_delayed_weight_count += opportunity; }
        void update_delayed_weight_opportunity( double opportunity ) { total_delayed_weight_opportunity += opportunity; }
        //
        // Members
        //
        /*! The dimension of total_coal_count, total_weighted_coal_opportunity, total_recomb_count, total_weighted_recomb_opportunity is
         *      number_of_epochs * number_of_population
         *  The dimension of total_mig_count is
         *      number_of_epochs * number_of_population (from) * number_of_population (to).
         *  For total_weighted_mig_opportunity only the 'from' population is important
         */
        vector < vector <double> >   total_coal_count;
        vector < vector <double> >   total_coal_opportunity;
        vector < vector <double> >   total_coal_weight;
        vector < vector <double> >   total_recomb_count;
        vector < vector <double> >   total_recomb_opportunity;
        vector < vector <double> >   total_recomb_weight;
        vector < vector < vector<double> > >  total_mig_count;
        vector < vector <double> >            total_mig_opportunity;
        vector < vector <double> >            total_mig_weight;
        vector < double >                local_recomb_opportunity;
        vector < vector < double > >     local_recomb_counts;
        double                           total_delayed_weight_opportunity;
        double                           total_delayed_weight_count;
        double                           total_resample_count;

        vector < double > counted_to;
        vector < double > lags;

        double const_lag_;
        double const_minimal_lag_update_ratio_;
        double update_param_threshold_;
        double update_param_interval_;

        double local_recording_interval_ = 100;
        int sample_size_;
    
        // DEBUG
        void print_pop_size();
        void print_change_time();
        void print_coal_count();
        bool print_mig_rate ( vector <vector<double>> & rates_list );
        void check_model_updated_mig(Model * model);
        void check_model_updated_Ne(Model * model);
    };
#endif
