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

#include"count.hpp"


void CountModel::init() {

    this->init_coal_and_recomb();
    this->init_migr();
    this->init_lags();

    this->update_param_interval_  = 5e6; // ONLINE EM
    this->update_param_threshold_ = 1e7; // ONLINE EM

}


void CountModel::reset_model_parameters(double current_base, Model * model, bool online, bool force_update, bool print){

    bool update = false;
    if ( online && current_base > this->update_param_threshold_ ) {
        update = true;
        this->update_param_threshold_ += this->update_param_interval_;
    }

    if ( update || force_update  ) {
        cout<<" MODEL IS RESET at base " << current_base <<endl;
        this->reset_recomb_rate ( model );
        this->reset_Ne ( model );
        this->reset_mig_rate ( model );
        model->finalize();
    } else {
        if (print) {
            this->print_coal_count();
        }
    }
}


void CountModel::log_counts( PfParam& param ) {
    // log coalescent counts
    for (size_t epoch_idx = 0; epoch_idx < change_times_.size(); epoch_idx++ ) {
        for (size_t pop_idx = 0; pop_idx < this->population_number(); pop_idx++ ) {
            param.append_to_count_file( epoch_idx, "Coal", pop_idx, -1, this->total_coal_opportunity[epoch_idx][pop_idx].final_answer(),
                                                                        this->total_coal_count[epoch_idx][pop_idx].final_answer(),
                                                                        this->total_coal_weight[epoch_idx][pop_idx].final_answer());
        }
    }
    // log recombination counts
    for (size_t epoch_idx = 0; epoch_idx < change_times_.size(); epoch_idx++ ) {
        param.append_to_count_file( epoch_idx, "Recomb", -1, -1, this->total_recomb_opportunity[epoch_idx][0].final_answer(),
                                                                 this->total_recomb_count[epoch_idx][0].final_answer(),
                                                                 this->total_recomb_weight[epoch_idx][0].final_answer());
    }
    // log migration counts
    for (size_t epoch_idx = 0; epoch_idx < change_times_.size(); epoch_idx++ ) {
        for (size_t from_pop_idx = 0; from_pop_idx < this->population_number(); from_pop_idx++ ) {
            for (size_t to_pop_idx = 0; to_pop_idx < this->population_number(); to_pop_idx++ ) {
                if (from_pop_idx != to_pop_idx) {
                    param.append_to_count_file( epoch_idx, "Migr", from_pop_idx, to_pop_idx, this->total_mig_opportunity[epoch_idx][from_pop_idx].final_answer(),
                                                                                             this->total_mig_count[epoch_idx][from_pop_idx][to_pop_idx].final_answer(),
                                                                                             this->total_mig_weight[epoch_idx][from_pop_idx].final_answer());
                }
            }
        }
    }
}


void CountModel::init_coal_and_recomb() {

    total_coal_count.clear();
    total_coal_opportunity.clear();
    total_coal_weight.clear();
    total_recomb_count.clear();
    total_recomb_opportunity.clear();
    total_recomb_weight.clear();

    resetTime();
    for (size_t epoch_idx = 0 ; epoch_idx < change_times_.size(); epoch_idx++) {
        // move to the next epoch
        if (epoch_idx > 0)
            increaseTime();
        // populate coalescent and recombination event counters
        vector <Two_doubles> tmp_count(this->population_number(), Two_doubles(0));
        total_coal_count.push_back(tmp_count);
        total_recomb_count.push_back(tmp_count);
        // enter initial value
        for (size_t pop_i = 0 ; pop_i < this->population_number(); pop_i++ ) {
            total_coal_count[epoch_idx][pop_i] = 1 / ( 2 * this->population_size( pop_i ) );
            /*! Note that this uses recombination rate at position -1 */
            total_recomb_count[epoch_idx][pop_i] = this->recombination_rate();
        }
        // populate and enter initial value for opportunity
        vector <Two_doubles> tmp_opportunity(this->population_number(), Two_doubles(1));
        total_coal_opportunity.push_back(tmp_opportunity);
        total_recomb_opportunity.push_back(tmp_opportunity);
        // and same for weight counters
        total_coal_weight.push_back(tmp_opportunity);
        total_recomb_weight.push_back(tmp_opportunity);
    }
}


void CountModel::init_migr() {

    total_mig_count.clear();
    total_mig_opportunity.clear();
    total_mig_weight.clear();

    // set initial counts/rates for all epochs
    resetTime();
    for (size_t epoch_idx = 0 ; epoch_idx < change_times_.size(); epoch_idx++) {

        // move to the next epoch
        if (epoch_idx > 0)
            increaseTime();

        // populate and set up initial values for the event count, opportunity, and inferred rate vectors, for one epoch
        vector < vector < Two_doubles > > tmp_count_Time_i;               // Event counts for migrations pop_i -> pop_j
        vector < Two_doubles >            tmp_opp_Time_i;                 // Opportunity  for migrations from pop_i
        vector < vector < double > >      tmp_rate_Time_i_double;         // Rates        for migrations pop_i -> pop_j
        for (size_t pop_i = 0 ; pop_i < this->population_number(); pop_i++ ){
            tmp_count_Time_i.      push_back( vector<Two_doubles>( this->population_number(), Two_doubles( 0.0 ) ) );
            tmp_opp_Time_i.        push_back( Two_doubles(1) );
            tmp_rate_Time_i_double.push_back( vector<double>( this->population_number(), 0 ) );
            for (size_t pop_j = 0 ; pop_j < this->population_number(); pop_j++) {
                tmp_count_Time_i[ pop_i ][ pop_j ]       = this->migration_rate( pop_i, pop_j );
                tmp_rate_Time_i_double[ pop_i ][ pop_j ] = this->migration_rate( pop_i, pop_j );
            }
        }
        total_mig_count.      push_back(tmp_count_Time_i);
        total_mig_opportunity.push_back(tmp_opp_Time_i);
        total_mig_weight.     push_back(tmp_opp_Time_i);
    }
}


void CountModel::init_lags(){

    this->counted_to.clear();
    this->lags.clear();
    for (size_t epoch_idx = 0 ; epoch_idx < change_times_.size(); epoch_idx++){
        this->counted_to.push_back( (double)0 );
        double top_t = epoch_idx == (change_times_.size() -1) ? change_times_[change_times_.size()-1] : change_times_[epoch_idx+1];
        double lag_i = this->const_lag_ > 0 ? this->const_lag_ : double(4) / (this->recombination_rate() * top_t) ;
        this->lags.push_back( lag_i );
    }
}

// Do we want to allow for different lags in diff populations?
void CountModel::reset_lag ( std::vector<double> new_lags ){
    for (size_t epoch_idx = 0; epoch_idx < change_times_.size(); epoch_idx++) {
        lags.at(epoch_idx) = new_lags.at(epoch_idx);
    }
}

void CountModel::reset_Ne ( Model *model ){

    for (size_t epoch_idx = 0; epoch_idx < change_times_.size(); epoch_idx++) {
        for (size_t pop_j = 0 ; pop_j < this->population_number(); pop_j++ ) {
            double coal_opp    = total_coal_opportunity[epoch_idx][pop_j].final_answer();
            double coal_count  = total_coal_count[epoch_idx][pop_j].final_answer();
            double coal_weight = total_coal_weight[epoch_idx][pop_j].final_answer();
            double coal_rate  = coal_count / coal_opp;
            double pop_size   = 1.0 / (2.0 * coal_rate);
            model->addPopulationSize(change_times_[epoch_idx], pop_j, pop_size ,false, false);
            cout << " Setting size of population " << pop_j << " @ " << setw(8) << change_times_[epoch_idx] << " to "
                 << setw(8) << pop_size
                 << " ( 0.5 * " << total_coal_opportunity[epoch_idx][pop_j].final_answer() << " / " << this->total_coal_count[epoch_idx][pop_j].final_answer() << "; post-lag ESS "
                 << 1.0 / (coal_weight / coal_opp) << " )" << endl;
        }
    }
    this->check_model_updated_Ne( model );
}


void CountModel::reset_recomb_rate ( Model *model ){

    double recomb_opportunity = 0;
    double recomb_count = 0;
    double recomb_weight = 0;

    for ( size_t epoch_idx = 0; epoch_idx < change_times_.size(); epoch_idx++ ) {
        recomb_opportunity += total_recomb_opportunity[epoch_idx][0].final_answer();
        recomb_count       += total_recomb_count[epoch_idx][0].final_answer();
        recomb_weight      += total_recomb_weight[epoch_idx][0].final_answer();
    }
    double inferred_recomb_rate = recomb_count / recomb_opportunity;

    model->setRecombinationRate( inferred_recomb_rate, false, false, 0);
    cout << " Setting recombination rate to " << model->recombination_rate(0)
         << " ( " << recomb_count<<" / " << recomb_opportunity << "; post-lag ESS "
         << 1.0 / (recomb_weight / recomb_opportunity) << " )" << endl;
}


void CountModel::clear_2d_vector ( vector <vector<double>> & rates_list ) {
    for (size_t i = 0; i < rates_list.size(); i++ ) {
        if (!rates_list[i].empty()) {
            for (size_t j = 0; j < rates_list[i].size() ; j++) {
                rates_list[i][j] = 0;
            }
        }
    }
}


void CountModel::reset_mig_rate ( Model *model ) {

    if (!has_migration()) return;

    clear_2d_vector( model->mig_rates_list_ );
    clear_2d_vector( model->total_mig_rates_list_ );

    vector < vector < vector<double> > > inferred_mig_rate;
    for (size_t epoch_idx = 0; epoch_idx<change_times_.size(); epoch_idx++) {
        for (size_t pop_i = 0 ; pop_i < this->population_number(); pop_i++ ) {
            for (size_t pop_j = 0 ; pop_j < this->population_number(); pop_j++ ) {
                if ( pop_i == pop_j) continue;
                double migration_rate =
                    total_mig_count[epoch_idx][pop_i][pop_j].final_answer() /
                    total_mig_opportunity[epoch_idx][pop_i].final_answer();

                model->addMigrationRate(change_times_[epoch_idx], pop_i, pop_j, migration_rate, false, false);
            }
        }
    }

    this->check_model_updated_mig (model);

    assert( print_mig_rate (model->mig_rates_list_) );
    assert( print_mig_rate (model->total_mig_rates_list_) );
}
