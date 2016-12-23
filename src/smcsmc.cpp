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

#include "general.hpp"
#include "count.hpp"
#include "arena.hpp"
#include "usage.hpp"
#include "model_summary.hpp"
#include <set>


/*!
 * Global variable for the memory arena
 */
class Arena* Arena::globalArena;


void pfARG_core(PfParam &pfARG_para,
                CountModel *countNe,
                bool print_update_count);


int main(int argc, char *argv[]){

    bool print_update_count = false; // DEBUG

    try {
        /*! Extract pfARG parameters */
        PfParam pfARG_para;
        (void)pfARG_para.parse( argc, argv );

        if (pfARG_para.version()){
            pfARG_para.printVersion(&std::cout);
            return EXIT_SUCCESS;
        }

        if (pfARG_para.help()){
            pfARG_para.printHelp();
            return EXIT_SUCCESS;
        }

        /*! Initialize memory arena */
        Arena* arena = new Arena( pfARG_para.model.getNumEpochs(), 1000000 );

        /*!  INITIALIZE CountModel */
        CountModel *countNe = new CountModel( pfARG_para.model , pfARG_para.lag);

        /*! EM step */
        pfARG_para.outFileHeader();
        for (int i = 0; i <= pfARG_para.EM_steps; i++) {

            cout << "EM step " << i << endl;
            clog << "EM step " << i << endl;
            pfARG_core( pfARG_para,
                        countNe,
                        print_update_count);
            pfARG_para.increaseEMcounter();
            cout << "End of EM step " << i << endl;
            clog << "End of EM step " << i << endl;
        }

        int exit_success = pfARG_para.log( );

        /*! Clean up */
        delete countNe;
        delete arena;
        return exit_success;
    }
    catch (const exception &e) {
        std::cerr << "Error: " << e.what() << std::endl;
        std::cout << "Error: " << e.what() << std::endl;
        return EXIT_FAILURE;
    }
}


bool is_height_in_tree( const Forest& forest, double height ) {
    for (auto it = forest.getNodes()->iterator(); it.good(); ++it) {
        if ( (abs( it.height() - height ) < 1e-6) && ((*it)->countChildren( true ) > 0 ) ){
            return true;
        }
    }
    return false;
}

vector<double> calculate_median_survival_distances( Model model, int min_num_events_per_epoch = 200 ) {

    // we want to change model attributes inside this function only
    model.biased_sampling = false;

    // Simulating trees in order to calibrate lag and bias ratios
    const int num_epochs = model.change_times().size();
    vector < vector <double> > survival_distances( num_epochs );  // survival_distances[epoch_idx][coal_sample_index]
    int num_epochs_not_done = num_epochs;
    MersenneTwister randomgenerator( true, 1 );
    int iterations = 0;
    const int max_iterations = 10000000;
    const int max_safe_sample_size = 10;
    size_t epoch_idx;
    
    while ( num_epochs_not_done > 0 ) {

      if (++iterations > max_iterations) break;
      Forest arg( &model, &randomgenerator );
      bool has_node_in_incomplete_epoch = false;
      arg.buildInitialTree();
      std::set<double> original_node_heights;
      for ( auto it = arg.getNodes()->iterator(); it.good(); ++it) {
          if (!(*it)->in_sample()) {
              double height = (*it)->height();
              original_node_heights.insert( height );
              for (epoch_idx = 0;
                   epoch_idx+1 < model.change_times().size() && model.change_times()[ epoch_idx+1 ] <= height;
                   ++epoch_idx) {}
              if (survival_distances[epoch_idx].size() < min_num_events_per_epoch)
                  has_node_in_incomplete_epoch = true;
          }
      }
      // if this tree cannot contribute to an epoch that needs more counts, then try the next one
      if (!has_node_in_incomplete_epoch) continue;

      // evolve the tree, and see how long nodes at particular heights survive
      while( !original_node_heights.empty() && arg.next_base() < (model.loci_length() * 0.6) ) {
        arg.sampleNextGenealogy( false );
        arg.sampleRecSeqPosition( false );
        for( auto it = original_node_heights.begin(); it != original_node_heights.end(); it++){
          if( !is_height_in_tree( arg, *it ) ) {
            // add genomic position to the histogram fo survival for the appropriate epoch
            for (epoch_idx = 0;
                 epoch_idx+1 < model.change_times().size() && model.change_times()[ epoch_idx+1 ] <= *it;
                 ++epoch_idx) {}
            survival_distances[epoch_idx].push_back( arg.current_base() );
            if (survival_distances[epoch_idx].size() == min_num_events_per_epoch) {
                --num_epochs_not_done;
            }
            // remove height from the set
            original_node_heights.erase( *it );
            break;
          }
        }
      }
    }
    clog << "Median survival: " << iterations << " iterations; per-epoch counts";
    for (size_t idx = 0; idx < survival_distances.size(); idx++ ) clog << " " << survival_distances[idx].size();
    clog << endl;

    // now deal with the possiblity that no or insufficient events have been sampled for an epoch
    vector <double> median_survival_distances;
    double earliest_survival_distance = -1;
    for(epoch_idx = 0; epoch_idx < survival_distances.size(); epoch_idx++ ) {
        std::sort( survival_distances[epoch_idx].begin(), survival_distances[epoch_idx].end() );
        int median_idx = (survival_distances[epoch_idx].size()-1) / 2;
        if (median_idx < max_safe_sample_size) {
            // insufficiently many samples for a borderline-reliable estimate
            median_survival_distances.push_back( -1 );
        } else {
            median_survival_distances.push_back( survival_distances[epoch_idx][ median_idx ] );
            if (earliest_survival_distance < 0)
                earliest_survival_distance = median_survival_distances.back();
        }
    }
    // supply reasonable defaults for epochs with insufficient samples
    for( epoch_idx = 0; epoch_idx < survival_distances.size(); epoch_idx++ ) {
        if (median_survival_distances[epoch_idx] < 0) {
            if (epoch_idx > 0) {
                median_survival_distances[epoch_idx] = median_survival_distances[epoch_idx - 1];
            } else {
                median_survival_distances[epoch_idx] = earliest_survival_distance;
            }
        }
        clog << " Epoch " << epoch_idx << " contains " << survival_distances[epoch_idx].size() << " coalescent events" << endl;
        clog << "   Survival distance set to " << median_survival_distances[epoch_idx] << endl;
    }
    return median_survival_distances;
}


void calibrate_bias_ratios(Model* model, double top_t_value) {

    if (!model->biased_sampling) return;
    
    clog << "model has bias heights " << model->bias_heights() << endl;
    clog << "model has bias strengths " << model->bias_strengths() << endl;

    int Num_bias_ratio_simulation_trees = 10000;
    ModelSummary model_summary = ModelSummary(model, top_t_value);
    for(size_t tree_idx = 0 ; tree_idx < Num_bias_ratio_simulation_trees ; tree_idx++){
      model_summary.addTree();
    }
    model_summary.finalize();
    
    clog << "Information from pre-sequence-analysis tree simulation:" << endl;
    clog << "    model_summary.times_: " << model_summary.times_ << endl;
    clog << "    avg B: " << model_summary.avg_B() << endl;
    clog << "    avg B below: " << model_summary.avg_B_below() << endl;
    clog << "    avg B within: " << model_summary.avg_B_within() << endl;
    clog << "    avg B below bh: " << model_summary.avg_B_below_bh() << endl;
    clog << "    avg lineage count: " << model_summary.avg_lineage_count() << endl;
    clog << "    single lineage count: " << model_summary.single_lineage_count() << endl;
    clog << "    tree count: " << model_summary.tree_count_ << endl;

    model->clearBiasRatios();
    for( size_t idx=0; idx < model->bias_strengths().size(); idx++ ){
        model->addToBiasRatios( model_summary.getBiasRatio(idx) );
    }
    clog << "    Bias ratios set to: " << model->bias_ratios() << endl;
    assert( model->bias_ratios().size() == model->bias_heights().size()-1 );
}


void pfARG_core(PfParam &pfARG_para,
                CountModel *countNe,
                bool print_update_count ) {

    Model *model         = &pfARG_para.model;
    MersenneTwister *rg  = pfARG_para.rg;
    size_t Nparticles    = pfARG_para.N ;
    Segment *Segfile     = pfARG_para.Segfile;
    double mutation_rate = model->mutation_rate();

    calibrate_bias_ratios( model, pfARG_para.top_t() );

    dout<< "############# starting seg file at base " << Segfile->segment_start()<<endl;    
    Segfile->read_new_line();
    /*! Initial particles */
    ParticleContainer current_states(model, rg, pfARG_para.record_event_in_epoch,
                                     Nparticles,
                                     Segfile->segment_start(),
                                     pfARG_para.Segfile->empty_file(),
                                     pfARG_para.Segfile->allelic_state_at_Segment_end);
    dout<<"######### finished initial particle building"<<endl;

    valarray<int> sample_count( Nparticles ); // if sample_count is in the while loop, this is the initializing step...

    /*! Initialize prior Ne */
    countNe->init();
    vector<double> median_survival = calculate_median_survival_distances( *model );

    // lags_to_application_delays sets app delay to half the argument, we want the app delay to be half the survival
    model->lags_to_application_delays( median_survival );
    for( size_t i=0; i<model->application_delays.size(); i++ ){
	    clog << " Application delay for epoch " << i << " set to " << model->application_delays.at(i) << endl;
    }
    // Now set the lag correctly according to the command
    if ( pfARG_para.calibrate_lag ){
      countNe->reset_lag( median_survival, pfARG_para.lag_fraction );
    }
    clog << "    Lags set to: " << countNe->check_lags() << endl;


    /*! Go through seg data */
    bool force_update = false;
    do{
        dout <<"current base is at "<<Segfile->segment_start()<<" and it is "
            << (Segfile->end_data()? "":"NOT") << " the end of data "<<endl;

        valarray<double> weight_cum_sum((Nparticles+1)); //Initialize the weight cumulated sums

        /*!     Sample the next genealogy, before the new data entry is updated to the particles
         *      In this case, we will be update till Segfile->site()
         */
        current_states.update_state_to_data( mutation_rate, (double)model->loci_length(), Segfile, weight_cum_sum);

        /*! Add posterior event counts to global counters */
        countNe->extract_and_update_count( current_states , min(Segfile->segment_end(), (double)model->loci_length()) );

        /*! Reset population sizes in the model */
        countNe->reset_model_parameters( min(Segfile->segment_end(), (double)model->loci_length()), model,
                                         pfARG_para.useCap, pfARG_para.Ne_cap, pfARG_para.online_bool,
                                         force_update = false, false );

        if ( pfARG_para.ESS() == 1 ) {
            dout << " random weights" <<endl;
            current_states.set_particles_with_random_weight();
        }
        /*! ESS resampling. Filtering step*/
        current_states.ESS_resampling(weight_cum_sum, sample_count, min(Segfile->segment_end(),
                                                                        (double)model->loci_length()),
                                      pfARG_para, Nparticles);

        if ( Segfile->segment_end() >= (double)model->loci_length() ) {
            cout << "\r" << " Particle filtering step" << setw(4) << 100 << "% completed." << endl;
            if ( Segfile->segment_end() > (double)model->loci_length() ) {
                clog << "  Segment data is beyond loci length" << endl;
            }
            Segfile->set_end_data (true);
        }

        Segfile->read_new_line(); // Read new line from the seg file

    } while( !Segfile->end_data() );

    clog << "Got to end of sequence" << endl;

    double sequence_end = pfARG_para.default_loci_length; // Set the sequence_end to the end of the sequence

    dout << endl << "### PROGRESS: end of the sequence at "<< sequence_end << endl;

    // This is mandatory, as the previous resampling step will set particle probabilities to ones.
    current_states.normalize_probability();

    countNe->extract_and_update_count( current_states , sequence_end, true );
    clog << " Inference step completed." << endl;

    ofstream recombFile(pfARG_para.get_recombination_map_filename(), ios::out | ios::app | ios::binary);
    countNe->dump_local_recomb_logs( recombFile, model->loci_length(), pfARG_para.EMcounter() );
    recombFile.close();
    
    countNe->reset_model_parameters(sequence_end, model, pfARG_para.useCap, pfARG_para.Ne_cap,
                                    true, force_update = true, true); // This is mandatory for EM steps

    countNe->log_counts( pfARG_para );

    clog << " Estimates log likelihood: " << current_states.ln_normalization_factor() << endl;
    current_states.clear();
    Segfile->reset_data_to_first_entry();

}
