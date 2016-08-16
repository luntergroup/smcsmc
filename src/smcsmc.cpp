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

#include "general.hpp"
#include "count.hpp"
#include "arena.hpp"
#include "debug/usage.hpp"
#include "help.hpp"
#include "model_summary.hpp"
#include <set>

#ifdef GPERF
#include <gperftools/heap-profiler.h>
#endif

/*!
 * Global variables for debugging the number of times that ForestState was constructed and removed.
 */
int new_forest_counter    = 0; // DEBUG
int delete_forest_counter = 0; // DEBUG
int recombination_counter = 0; // DEBUG
int recombination_event_called = 0;
double recomb_opp = 0; // DEBUG
int node_created = 0;
int node_deleted = 0;

/*!
 * Global variable for the memory arena
 */
class Arena* Arena::globalArena;


void pfARG_core(PfParam &pfARG_para,
                CountModel *countNe,
                bool print_update_count);


int main(int argc, char *argv[]){

    bool print_update_count = false; // DEBUG

    if ( argc == 1 ) {
        Help_header();
    }

    #ifdef GPERF
    //std::cout << "Starting heap profiling..." << std::endl;
    //HeapProfilerStart("smcsmc_gprofile");
    #endif    
    
    try {
        /*! Extract pfARG parameters */
        PfParam pfARG_para( argc, argv );

        /*! Initialize memory arena */
        Arena* arena = new Arena( pfARG_para.model.getNumEpochs(), 1000000 );

        /*!  INITIALIZE CountModel */
        CountModel *countNe = new CountModel( pfARG_para.model , pfARG_para.lag);
        pfARG_para.appending_Ne_file( true ); // Append initial values to History file

        /*! EM step */
        for (int i = 0; i <= pfARG_para.EM_steps; i++) {

            cout << "EM step " << i << endl;
            pfARG_core( pfARG_para,
                        countNe,
                        print_update_count);
            cout << "End of EM step " << i << endl;
        }

        pfARG_para.appending_Ne_file( );

        int exit_success = pfARG_para.log( );

        /*! Clean up */
        delete countNe;
        delete arena;
        cout << "Actual recombination "<<recombination_counter<<endl;// DEBUG
	#ifdef GPERF
	//HeapProfilerStop();
	//std::cout << "Stopped heap profiling." << std::endl;
	#endif
        return exit_success;
    }
    catch (const exception &e) {
        std::cerr << "Error: " << e.what() << std::endl;
	std::cout << "Error: " << e.what() << std::endl;
        Help_option();
	#ifdef GPERF
	//HeapProfilerStop();
	//std::cout << "Stopped heap profiling." << std::endl;
	#endif
        return EXIT_FAILURE;
    }
}


bool is_height_in_tree( const Forest& forest, double height ) {
    for (auto it = forest.getNodes()->iterator(); it.good(); ++it) {
        if ( abs( it.height() - height ) <= .000001 && it.local() ){
            return true;
        }
    }
    return false;
}

vector<double> calculate_median_survival_distances( Model model, int Num_events = 200 ) {

    // we want to change model attributes inside this function only
    model.biased_sampling = false;

    //// Simulating trees in order to calibrate lag and bias ratios
    // lag calibration
    const int num_epochs = model.change_times().size();
    //const int num_pops = model.population_number();
    //cout << "Number of populations: " << num_pops << endl;
    //cout << "std::vector<std::vector<double> > single_mig_probs_list_; " << model.single_mig_probs_list_ << endl;
    //cout << "single_mig_probs_list_.size() " << model.single_mig_probs_list_.size() << endl;
    //for (size_t i = 0; i<model.change_times().size(); i++){
        //cout << "single_mig_probs_list_.at(" << i << ") " << model.single_mig_probs_list_.at(i) << endl;
    //}
    ////vector < vector < vector<double> > > survival_distances( num_pops ); // survival_distances[pop_idx][epoch_idx][coal_sample_index]
    //// to size each vector to number of epochs, need to find which one involves the absorbtion of that pop...
    vector < vector <double> > survival_distances( num_epochs );  // survival_distances[epoch_idx][coal_sample_index]
    int num_epochs_not_done = num_epochs;
    MersenneTwister randomgenerator( true, 1 );
    while ( num_epochs_not_done > 0 ) {
      Forest arg( &model, &randomgenerator );
      arg.buildInitialTree();
      std::set<double> original_node_heights;
      for ( auto it = arg.getNodes()->iterator(); it.good(); ++it) {
          if (!(*it)->in_sample()) {
              original_node_heights.insert( (*it)->height() );
          }
      }
      while( !original_node_heights.empty() && arg.next_base() < (model.loci_length() * 0.6) ) {
        arg.sampleNextGenealogy( false );
        arg.sampleRecSeqPosition( false );
        for( auto it = original_node_heights.begin(); it != original_node_heights.end(); it++){
          if( !is_height_in_tree( arg, *it ) ) {
            //// add genomic position to the histogram fo survival for the appropriate epoch
            // loop over change times to get appropriate index
            size_t epoch_idx = 0;
            while (epoch_idx+1 < model.change_times().size() && model.change_times().at(epoch_idx+1) <= *it) {
                 epoch_idx++;
            }
            survival_distances[epoch_idx].push_back( arg.current_base() );
            if (survival_distances[epoch_idx].size() == Num_events) {
                --num_epochs_not_done;
            }
            //// remove height from the set
            original_node_heights.erase( *it );
            break;
          }
        }
      }
      // Need to break if we have troublesome epochs due to large estimated population sizes
      bool exists_highly_explored_epoch = false;
      bool exists_insufficiently_explored_epoch = false;
      for( size_t epoch_idx = 0; epoch_idx < num_epochs; epoch_idx++) {
        if( survival_distances[epoch_idx].size() > 1000000 ) {
          exists_highly_explored_epoch = true;
        }
        if( survival_distances[epoch_idx].size() < 3 ) {
          exists_insufficiently_explored_epoch = true;
        }
      }
      if ( exists_highly_explored_epoch && !exists_insufficiently_explored_epoch ) {
        break;
      }
    }
    vector <double> median_survival;
    cout << "ARGs simulated under model:" << endl;
    for( size_t epoch_idx = 0; epoch_idx < survival_distances.size(); epoch_idx++ ) {
        std::sort( survival_distances[epoch_idx].begin(), survival_distances[epoch_idx].end() );
        int median_idx = (survival_distances[epoch_idx].size()-1) / 2;
        if (median_idx == -1) {
            throw std::length_error("No survival_distance instances in epoch");
        }
        median_survival.push_back( survival_distances[epoch_idx][ median_idx ] );
        cout << " Epoch " << epoch_idx << " contains " << survival_distances[epoch_idx].size() << " coalescent events" << endl;
        cout << "   These have a median survival distance of " << survival_distances[epoch_idx][ median_idx ]
             << ", a minimum of " << *survival_distances[epoch_idx].begin()
             << ", and a maximum of " << *survival_distances[epoch_idx].rbegin() << endl;
    }
    return median_survival;
}


void pfARG_core(PfParam &pfARG_para,
                CountModel *countNe,
                bool print_update_count ) {
    recombination_counter = 0; // DEBUG
    recomb_opp = 0;// DEBUG

    int who = RUSAGE_SELF;     // PROFILING
    struct rusage usage;       // PROFILING
    struct rusage *p = &usage; // PROFILING

    Model *model         = &pfARG_para.model;
    MersenneTwister *rg  = pfARG_para.rg;
    size_t Nparticles    = pfARG_para.N ;
    Segment *Segfile     = pfARG_para.Segfile;
    double mutation_rate = model->mutation_rate();
    cout << "model has bias heights " << model->bias_heights() << endl;
    cout << "model has bias strengths " << model->bias_strengths() << endl;

    // bias ratio calibration
    int Num_bias_ratio_simulation_trees = 10000;
    ModelSummary model_summary = ModelSummary(model, pfARG_para.top_t());
    for(size_t tree_idx = 0 ; tree_idx < Num_bias_ratio_simulation_trees ; tree_idx++){
      //cout << "Adding tree " << tree_idx << endl;
      model_summary.addTree();
      //cout << "current tree B " << model_summary.current_tree_B() << endl;
    }
    model_summary.finalize();
    cout << "Information from pre-sequence-analysis tree simulation:" << endl;
    cout << "    model_summary.times_: " << model_summary.times_ << endl;
    cout << "    avg B: " << model_summary.avg_B() << endl;
    cout << "    avg B below: " << model_summary.avg_B_below() << endl;
    cout << "    avg B within: " << model_summary.avg_B_within() << endl;
    cout << "    avg B below bh: " << model_summary.avg_B_below_bh() << endl;
    cout << "    avg lineage count: " << model_summary.avg_lineage_count() << endl;
    cout << "    single lineage count: " << model_summary.single_lineage_count() << endl;
    cout << "    tree count: " << model_summary.tree_count_ << endl;
    
    if(model->biased_sampling) {
        model->clearBiasRatios();
        for( size_t idx=0; idx < model->bias_strengths().size(); idx++ ){
            model->addToBiasRatios( model_summary.getBiasRatio(idx) );
        }
        cout << "    Bias ratios set to: " << model->bias_ratios() << endl;
        assert( model->bias_ratios().size() == model->bias_heights().size()-1 );
    }
    //// Done with calibration



    dout<< "############# starting seg file at base " << Segfile->segment_start()<<endl;
    Segfile->read_new_line();
    /*! Initial particles */
    //double initial_position = 0;
    ParticleContainer current_states(model, rg, pfARG_para.record_event_in_epoch,
                                     Nparticles,
                                     Segfile->segment_start(),
                                     pfARG_para.heat_bool,
                                     pfARG_para.Segfile->empty_file(),
                                     pfARG_para.Segfile->allelic_state_at_Segment_end);
    dout<<"######### finished initial particle building"<<endl;

    valarray<int> sample_count( Nparticles ); // if sample_count is in the while loop, this is the initializing step...
    //#ifdef _SCRM
    //current_states.print_particle_newick();
    //#endif
    /*! Initialize prior Ne */
    countNe->init();
    vector<double> median_survival = calculate_median_survival_distances( *model );
    if ( pfARG_para.calibrate_lag ){
      countNe->reset_lag( median_survival, pfARG_para.lag_fraction );
    }
    model->lags_to_application_delays( countNe->check_lags() );
    for( size_t i=0; i<model->application_delays.size(); i++ ){
	    cout << " Application delay for epoch " << i << " set to " << model->application_delays.at(i) << endl;
	}
    cout << "    Lags set to: " << countNe->check_lags() << endl;

    //cout << "This model has " << model->population_number() << " populations" << endl;
    //model->resetTime();
    //for( size_t idx=0; idx < model->change_times().size(); idx++){
        //cout << "hasFixedTimeEvent is " << model->hasFixedTimeEvent( model->change_times()[idx] )
             //<< " for epoch " << idx << " at time " << model->change_times()[idx] << endl;
        //if( idx < model->change_times().size() - 1 ) {
            //model->increaseTime();
        //} else if( idx == model->change_times().size() - 1 ) {
            //model->resetTime();
        //} else {
            //cout << "PROBLEM ITERATING OVER TIMES" << endl;
        //}
    //}

    /*! Go through seg data */
    bool force_update = false;
    do{
        getrusage(who,p);            // PROFILING
        process(p, Segfile->segment_start()); // PROFILING

        /*! A particle path, where x-----o is a ForestState.
         *  The particle weight is the weight at the ForestState 6
         *  Even though the particle has extended to state 6,
         *  which means the Weight has included upto Segfile->site()
         *
         * The ForestState we are intested in is upon until backbase,
         * which is state 2 in this case.
         *
         * Update count step should include all the coalescent events of
         * State 0, 1, and 2
         *
         * backbase is the sequence position that the current mutation "Segfile->site()" go back lag
         *
         * previous_backbase is the sequence position that the previous mutation "Segfile->site()" go back lag
         *
         *  \verbatim

           remove all the state prior to the minimum of
            current_printing_base and
                previous_backbase
                .                      backbase                     Segfile->site()
                .                      .                            .
                .                      .     3                      .
                .                      .     x---o              6   .
                .                  2   .     |   |              x-------o
                .                  x---------o   |              |   .
                .                  |   .         |              |   .
             0  .                  |   .         x---o          |   .
             x---------o           |   .         4   |          |   .
                .      |           |   .             x----------o   .
                .      |           |   .             5              .
                .      x-----------o   .                            .
                .      1               .-------------lag------------.
                .                      .                            .
                Count::update_e_count( .                            ParticleContainer::update_state_to_data(
           x_start = previous_backbase .                              mutation data comes in here
           x_end = backbase            .
                                       .
                                       ParticleContainer::appendingStuffToFile(
                                           x_end = backbase,
          \endverbatim
         */

        dout <<"current base is at "<<Segfile->segment_start()<<" and it is ";
        dout << (Segfile->end_data()? "YES":"NOT");
        dout << " the end of data "<<endl;

        valarray<double> weight_cum_sum((Nparticles+1)); //Initialize the weight cumulated sums

        /*!     Sample the next genealogy, before the new data entry is updated to the particles
         *      In this case, we will be update till Segfile->site()
         */
        current_states.update_state_to_data( mutation_rate, (double)model->loci_length(), Segfile, weight_cum_sum);

        /*! WRITE TMRCA AND BL TO FILE, This is used when generating the heatmap */
        current_states.appendingStuffToFile( min(Segfile->segment_end(), (double)model->loci_length()), pfARG_para);

        /*! UPDATE CUM COUNT AND OPPORTUNITIES ACCORDING TO THE PARTICLE WEIGHT */
        countNe->extract_and_update_count( current_states , min(Segfile->segment_end(), (double)model->loci_length()) );

        /*! Reset population sizes in the model */
        countNe->reset_model_parameters( min(Segfile->segment_end(), (double)model->loci_length()), model, pfARG_para.useCap, pfARG_para.Ne_cap, pfARG_para.online_bool, force_update = false, false );

        if ( pfARG_para.ESS() == 1 ) {
            dout << " random weights" <<endl;
            current_states.set_particles_with_random_weight();
        }
        /*! ESS resampling. Filtering step*/
        current_states.ESS_resampling(weight_cum_sum, sample_count, min(Segfile->segment_end(), (double)model->loci_length()), pfARG_para.ESSthreshold, Nparticles);

        if ( Segfile->segment_end() >= (double)model->loci_length() ) {
            cout << "\r" << " Particle filtering step" << setw(4) << 100 << "% completed." << endl;
            if ( Segfile->segment_end() > (double)model->loci_length() ) {
                cout << "  Segment data is beyond loci length" << endl;
            }
            Segfile->set_end_data (true);
        }

        Segfile->read_new_line(); // Read new line from the seg file

    } while( !Segfile->end_data() );

    cout << "Got to end of sequence" << endl;

    double sequence_end = pfARG_para.default_loci_length; // Set the sequence_end to the end of the sequence
    // EXDEND THE ARG TO THE END OF THE SEQUENCE AS MISSING DATA ...
    //bool with_data_to_the_end = true;
    //current_states.extend_ARGs( sequence_end, model->mutation_rate(),  with_data_to_the_end); // DEBUG
    // In case the rest of the sequence is too long, this needs some "ghost" snp ... invariants
    // should include coalescent events as well ...

    //current_states.cumulate_recomb_opportunity_at_seq_end( sequence_end ); // This is to make up the recomb opportunities till the end of the sequence.
    dout <<endl << "### PROGRESS: end of the sequence at "<< sequence_end << endl;

    // This is mandatory, as the previous resampling step will set particle probabilities to ones.
    current_states.normalize_probability();

    countNe->extract_and_update_count( current_states , sequence_end, true ); // Segfile->end_data()
    cout << " Inference step completed." << endl;
    countNe->reset_model_parameters(sequence_end, model, pfARG_para.useCap, pfARG_para.Ne_cap, true, force_update = true, true); // This is mandatory for EM steps

    bool append_to_history_file = true;
    pfARG_para.appending_Ne_file( append_to_history_file );
    countNe->log_counts( pfARG_para );

    /*! WRITE TMRCA AND BL TO FILE, This is used when generating the heatmap */
    current_states.appendingStuffToFile( sequence_end, pfARG_para);
    
    current_states.print_ln_normalization_factor();

    //#ifdef _SCRM
        //current_states.print_particle_newick();
    //#endif
    current_states.clear(); // This line is sufficient to clear the memory.
    Segfile->reset_data_to_first_entry();

    dout << "Actual recombination "<<recombination_counter<<endl;// DEBUG
    dout << "Actual recombination recorder called "<<recombination_event_called<<endl;// DEBUG
    dout << "Actual node_created "<<node_created<<endl;// DEBUG
    dout << "Actual node_deleted "<<node_deleted<<endl;// DEBUG
    dout << "Actual recomb op is "<<recomb_opp <<endl;

} // End of void pfARG_core( ... )
