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
#include "model_summary.hpp"

#include <set>
#include <fstream>
#include <iostream>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include <fenv.h>

// Global variable for the memory arena
class Arena* Arena::globalArena;


void pfARG_core(PfParam &pfARG_para,
                CountModel *countNe,
                bool print_update_count);


int main(int argc, char *argv[]){

    bool print_update_count = false; // DEBUG

    // catch floating point exceptions, to aid debugging, 
    // but not on OSX as it uses SSE for all FP math
#ifndef __APPLE__
    feenableexcept(FE_DIVBYZERO|FE_INVALID|FE_OVERFLOW);
#endif
    
    try {
        /*! Extract pfARG parameters and start timer */
        PfParam pfARG_para;
	pfARG_para.startTimer();
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

            clog << "EM step " << i << endl;
            pfARG_core( pfARG_para,
                        countNe,
                        print_update_count);
            pfARG_para.increaseEMcounter();
            clog << "End of EM step " << i << endl;
	    /*! Print EM step timer */
	    clog << "Elapsed time: " << setprecision(3) << pfARG_para.stopAndPrintTimer() << endl;
        }

        int exit_success = pfARG_para.log( );

        /*! Clean up */
        delete countNe;
        delete arena;
        return exit_success;
    }
    catch (const exception &e) {
        std::cerr << "Error: " << e.what() << std::endl;
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


double parent_height_ignoring_migrations( const Node* node ) {
    while (!node->is_root()) {
        node = node->parent();
        if ((!node->is_migrating()) && 
            (node->countChildren( /* only_local = */ true ) == 2 )) {
            break;
        }
    }
    return node->height();
}


TerminalBranchLengthQuantiles calculate_terminal_branch_length_quantiles( Model& model ) {

    const int max_iterations = 1000000;
    MersenneTwister randomgenerator( true, 1 );
    TerminalBranchLengthQuantiles tblq;
    double total_length = 0.0;
    tblq.quantiles = {0.001,0.003,0.01,0.03,0.1,0.5,0.95};
    tblq.lengths.resize( model.sample_size() );
    vector<vector<double>> tbls( model.sample_size() );

    // (for auxiliary particle filter:)
    cout << "Calculating terminal branch length quantiles..." << endl;
    
    for (int iterations = 0; iterations < max_iterations; iterations++) {

        Forest tree( &model, &randomgenerator );
        tree.buildInitialTree( false );
        total_length += tree.getLocalTreeLength();
        for ( auto it = tree.getNodes()->iterator(); it.good(); ++it) {
            if ((*it)->in_sample()) {
                double height = parent_height_ignoring_migrations( (*it) );
                int sample = (*it)->label() - 1;
                tbls[ sample ].push_back( height );
            }
        }
    }
    for (int sample = 0; sample < model.sample_size(); ++sample) {
        cout << "Lineage " << sample << "; Terminal branch length quantiles:";
        std::sort( tbls[ sample ].begin(), tbls[ sample ].end() );
        for (int q = 0; q < tblq.quantiles.size(); ++q) {
            tblq.lengths[ sample ].push_back( tbls[ sample ][ (int)(tblq.quantiles[q] * max_iterations) ] );
            cout << " [" << tblq.quantiles[q] << ":] " << tblq.lengths[ sample ].back();
        }
        cout << endl;
    }
    tblq.mean_total_branch_length = total_length / max_iterations;
    
    return tblq;
}


vector<double> calculate_median_survival_distances( Model model, int min_num_events_per_epoch = 200 ) {

    // we want to change model attributes inside this function only
    model.biased_sampling = false;
    const int num_epochs = model.change_times().size();
    const int max_iterations = 1000000;
    const long long max_num_trees = (long long)50 * max_iterations;
    const int max_safe_sample_size = 10;

    // Simulating trees in order to calibrate lag and bias ratios
    int iterations = 0;
    long long num_trees = 0;
    int num_epochs_not_done = num_epochs;
    size_t epoch_idx;
    vector < vector <double> > survival_distances( num_epochs );  // [epoch_idx][coal_sample_index]
    MersenneTwister randomgenerator( true, 1 );
    
    while ( num_epochs_not_done > 0 ) {

      if (++iterations > max_iterations) break;
      if (num_trees > max_num_trees) {
          throw std::runtime_error("Using an excessive number of trees to calculate survival distances."
                                   "  Check if model parameters are sensible.");
      }
      Forest arg( &model, &randomgenerator );
      bool has_node_in_incomplete_epoch = false;
      arg.buildInitialTree( false );
      std::set<double> original_node_heights;
      for ( auto it = arg.getNodes()->iterator(); it.good(); ++it) {
          if ((!(*it)->in_sample()) &&
              (!(*it)->is_migrating())) {
              double height = (*it)->height();
              original_node_heights.insert( height );
              for (epoch_idx = 0;
                   epoch_idx+1 < num_epochs && model.change_times()[ epoch_idx+1 ] <= height;
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
        ++num_trees;
        for( auto it = original_node_heights.begin(); it != original_node_heights.end(); it++){
          if( !is_height_in_tree( arg, *it ) ) {
            // add genomic position to the histogram fo survival for the appropriate epoch
            for (epoch_idx = 0;
                 epoch_idx+1 < num_epochs && model.change_times()[ epoch_idx+1 ] <= *it;
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
        clog << " Epoch " << epoch_idx << " (" << survival_distances[epoch_idx].size()
             << " events): survival distance " << median_survival_distances[epoch_idx] << endl;
    }
    return median_survival_distances;
}


int max_epoch_to_update( const vector<double>& lags, double distance_to_mutation ) {
    int epoch = 0;
    // There's not real justification for the factor 0.5 below; but large gaps in the data are infrequent,
    // and cause significant memory usage, so stop recording events a bit sooner than you would idealy want to do.
    double scale_factor = 0.5;
    while (epoch < lags.size() && 
           distance_to_mutation < scale_factor * lags[epoch])
        epoch++;
    return epoch - 1;
}


void pfARG_core(PfParam &pfARG_para,
                CountModel *countNe,
                bool print_update_count ) {

    Model *model         = &pfARG_para.model;
    MersenneTwister *rg  = pfARG_para.rg;
    size_t Nparticles    = pfARG_para.N ;
    Segment *Segfile     = pfARG_para.Segfile;

    vector<double> median_survival = calculate_median_survival_distances( *model );
    TerminalBranchLengthQuantiles tblq = calculate_terminal_branch_length_quantiles( *model );
    
    /* load the recombination guide rates into the model */
    pfARG_para.setModelRates();

    dout<< "############# starting seg file at base " << Segfile->segment_start()<<endl;
    Segfile->read_new_line();
    /*! Initial particles */
    ParticleContainer current_states(model, rg, &pfARG_para,
                                     Nparticles,
                                     Segfile->segment_start(),
                                     Segfile->empty_file(),
                                     Segfile->allelic_state_at_Segment_end);
    dout<<"######### finished initial particle building"<<endl;

    /*! Initialize prior Ne */
    countNe->init();

    // lags_to_application_delays sets app delay to half the argument, we want the app delay to be half the survival
    model->lags_to_application_delays( median_survival, pfARG_para.delay );
    if (pfARG_para.delay > 0) {
        for( size_t i=0; i<model->application_delays.size(); i++ ){
	    clog << " Application delay for epoch " << i << " set to " << model->application_delays.at(i) << endl;
        }
    }
    // Now set the lag correctly according to the command
    if ( pfARG_para.calibrate_lag ){
      countNe->reset_lag( median_survival, pfARG_para.lag_fraction );
    }
    const vector<double>& lags = countNe->get_lags();
    clog << "    Lags set to: " << lags << endl;
    clog << " Starting position: " << fixed << setprecision(0) << pfARG_para.start_position << setprecision(6) << scientific << endl;

    /*! Go through seg data */
    bool force_update = false;
    int total_resample_count = 0;
    do{
        dout <<"Current base is " << fixed << setprecision(0) << Segfile->segment_start() << setprecision(6) << scientific << endl;

        /*!     Sample the next genealogy, before the new data entry is updated to the particles
         *      In this case, we will be update till Segfile->site()
         */
        int max_update_epoch = max_epoch_to_update( lags, Segfile->distance_to_mutation() );
        current_states.update_state_to_data( Segfile, pfARG_para.ancestral_aware, tblq, max_update_epoch );

        /*! Add posterior event counts to global counters */
        countNe->extract_and_update_count( current_states , min(Segfile->segment_end(), (double)model->loci_length()) );

        /*! Reset population sizes in the model */
        countNe->reset_model_parameters( min(Segfile->segment_end(), (double)model->loci_length()), model,
                                         pfARG_para.useCap, pfARG_para.Ne_cap, pfARG_para.online_bool,
                                         force_update = false, false );

        if ( pfARG_para.get_ESS_fraction() == 1 ) {
            dout << " random weights" <<endl;
            current_states.set_particles_with_random_weight();
        }

        /*! Resample if the ESS becomes too low */
        double update_to = min( Segfile->segment_end(), (double)model->loci_length() );
        if (current_states.resample( update_to, pfARG_para )) {
            total_resample_count++;
            countNe->update_resample_count();
        }

        if ( Segfile->segment_end() >= (double)model->loci_length() ) {
            cout << "\r Particle filtering step 100% completed." << endl;
            Segfile->set_end_data (true);
        }

        Segfile->read_new_line(); // Read new line from the seg file

    } while( !Segfile->end_data() );

    clog << "Got to end of sequence; resampled " << total_resample_count << " times" << endl;

    current_states.print_recent_recombination_histogram();    //DEBUG
    
    double sequence_end = pfARG_para.default_loci_length; // Set the sequence_end to the end of the sequence

    dout << endl << "### PROGRESS: end of the sequence at "<< sequence_end << endl;

    // This is mandatory, as the previous resampling step will set particle probabilities to ones.
    current_states.normalize_probability();

    countNe->extract_and_update_count( current_states , sequence_end, true );
    clog << " Inference step completed." << endl;

    {
        // put in a scope to close the files automatically after dumping recomb logs
        ofstream recombFile(pfARG_para.get_recombination_map_filename(), ios::out | ios::app | ios::binary);
        boost::iostreams::filtering_ostream out;
        out.push( boost::iostreams::gzip_compressor() );
        out.push( recombFile );
        countNe->dump_local_recomb_logs( out, model->loci_length(), pfARG_para.EMcounter(), pfARG_para.start_position );
    }
    
    countNe->reset_model_parameters(sequence_end, model, pfARG_para.useCap, pfARG_para.Ne_cap,
                                    true, force_update = true, true); // This is mandatory for EM steps

    countNe->log_counts( pfARG_para );

    // write log likelihood to out file
    pfARG_para.appendToOutFile( pfARG_para.EMcounter(), -1, 0, 1e+99, "LogL", -1, -1, 1.0, current_states.ln_normalization_factor(), 1.0 );
    clog << " Estimated log likelihood: " << current_states.ln_normalization_factor() << endl;

    // write out trees
    current_states.resample( (double)model->loci_length(), pfARG_para, NULL, 1 );
    current_states.printTrees( pfARG_para );
    
    current_states.clear();
    Segfile->reset_data_to_first_entry();

}
