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
        return exit_success;
    }
    catch (const exception &e) {
        std::cerr << "Error: " << e.what() << std::endl;
        Help_option();
        return EXIT_FAILURE;
    }
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



    //// Simulating trees in order to calibrate lag and bias ratios
    ModelSummary model_summary = ModelSummary(model, pfARG_para.top_t());
    int Num_trees = 10000;

    for(size_t tree_idx = 0 ; tree_idx < Num_trees ; tree_idx++){
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
      model->setBiasRatioUpper( model_summary.getBiasRatioUpper() );
      model->setBiasRatioLower( model_summary.getBiasRatioLower() );
      cout << "    Bias ratio upper set to: " << model->bias_ratio_upper() << endl;
      cout << "    Bias ratio lower set to: " << model->bias_ratio_lower() << endl;
      //cout << "To insure correct num of rec events we need rho_u*B_u+rho_l*B_l=rho*B" << endl;
      //cout << "br_upper*B_upper " << model->bias_ratio_upper()*(model_summary.avg_B()-model_summary.avg_B_below_bh()) << endl;
      //cout << "br_lower*B_lower " << model->bias_ratio_lower()*model_summary.avg_B_below_bh() << endl;
      //cout << "B " << model_summary.avg_B() << endl;
    }
    // we reset the lags later as they are currently initialized later
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
    countNe->reset_lag(model_summary.getLags());
    cout << "    Lags set to: " << countNe->check_lags() << endl;

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
        countNe->reset_model_parameters( min(Segfile->segment_end(), (double)model->loci_length()), model, pfARG_para.online_bool, force_update = false, false);

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
    countNe->reset_model_parameters(sequence_end, model, true, force_update = true, true); // This is mandatory for EM steps

    bool append_to_history_file = true;
    pfARG_para.appending_Ne_file( append_to_history_file );
    countNe->log_counts( pfARG_para );

    /*! WRITE TMRCA AND BL TO FILE, This is used when generating the heatmap */
    current_states.appendingStuffToFile( sequence_end, pfARG_para);

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
