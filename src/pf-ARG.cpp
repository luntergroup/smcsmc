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
#include"vcf.hpp"
//#include"usage.hpp"
#include"count.hpp"
#include"pfparam.hpp"
using namespace pfARG;


/*!
 * Global variables for debugging the number of times that ForestState was constructed and removed.
 */ 
int new_forest_counter = 0;
int delete_forest_counter = 0;

void pfARG_core(Model *model, 
                MersenneTwister *rg,
                pfARG::param pfARG_para,
                Vcf *VCFfile,
                CountModel *countNe,
                bool print_update_count);

int main(int argc, char *argv[]){

    //
    // Settings
    //
    //bool show_progress = false;
    //INITIALIZE USAGE
    //int who = RUSAGE_SELF;
    //struct rusage usage;
    //struct rusage *p=&usage;
    
    /*! 
     * Default values
     */ 

    bool print_update_count = false;

    /*! Extract pfARG parameters */
    pfARG::param pfARG_para( argc, argv );
    
    if ( argc==1 ){
        //pfARG_para.print_help();
    } try {//else, proceed
        
        /*! 
         * INITIALIZE
         */

        /*! Initialize vcf file, and data up to the first data entry says "PASS"   */
        Vcf * VCFfile =  new Vcf(pfARG_para.vcf_NAME, pfARG_para.buff_length);
        pfARG_para.finalize ( VCFfile );

        /*! convert scrm_input string to argv */
        enum { kMaxArgs = 264 };
        int scrm_argc = 0;
        char *scrm_argv[kMaxArgs];        
        char * p2 = strtok((char *)pfARG_para.scrm_input.c_str(), " ");
        while (p2 && scrm_argc < kMaxArgs) {
            scrm_argv[scrm_argc++] = p2;
            p2 = strtok(0, " ");
            }
        
        /*! Extract scrm parameters */        
        Param * scrm_para = new Param(scrm_argc, scrm_argv, false);        
                
        /*! Initialize scrm model */
        Model * model = new Model();        
        scrm_para->parse( *model );

        /*! Initialize mersenneTwister seed */
        MersenneTwister *rg = new MersenneTwister(scrm_para->random_seed);
                                
        /*! Initialize updated weighted coalescent count and BL  */ 
        CountModel *countNe = new CountModel(*model);
        //countNe->check_CountModel_Ne();               
        
        /*!
         * PfARG CORE 
         */ 
        for (int I = 0; I <= pfARG_para.EM_steps; I++){
            cout << "Now starting EM_step " << I << endl;
            pfARG_core(model, rg, pfARG_para, VCFfile, countNe, print_update_count);
            cout << "End of EM_step " << I << endl;
        }
        
        //countNe->reset_model_Ne( model, true, true); // This is mandatory for appending the correct value out ...
        pfARG_para.appending_Ne_file(model);
        
        //int main_return = pfARG_para.log(model, rg->seed(), runningtime, countNe->inferred_recomb_rate);
        int main_return = pfARG_para.log(model, rg->seed(), countNe->inferred_recomb_rate);
        /*! Clean up */
        delete scrm_para;
        delete VCFfile;
        delete model; // 
        delete rg;
        delete countNe;
        
        cout<<"Forest state was created " << new_forest_counter << " times" << endl;
        dout<<"Forest state destructor was called " << delete_forest_counter << " times" << endl;
        
        return main_return;

    } catch (const exception &e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return EXIT_FAILURE;
    }
}


void pfARG_core(Model *model, 
                MersenneTwister *rg,
                pfARG::param pfARG_para,
                Vcf *VCFfile,
                CountModel *countNe,
                bool print_update_count){
    VCFfile->read_new_line();
    cout<< "############# starting vcf file at base " <<VCFfile->site()<<endl;
    /*! Initial particles */ 

    double initial_position = 0;

    ParticleContainer current_states(model, rg, pfARG_para.N, VCFfile->vec_of_sample_alt_bool, VCFfile->withdata(), initial_position );             

    valarray<int> sample_count(pfARG_para.N); // if sample_count is in the while loop, this is the initializing step...

    /*! Initialize prior Ne */
    countNe->init();

    /*! Go through vcf data */
    double previous_x = 0;
    do{
        
        /*! A particle path, where x-----o is a ForestState. 
         *  The particle weight is the weight at the ForestState 6
         *  Even though the particle has extended to state 6, 
         *  which means the Weight has included upto VCFfile->site()
         *  
         * The ForestState we are intested in is upon until backbase,
         * which is state 2 in this case.
         * 
         * Update count step should include all the coalescent events of 
         * State 0, 1, and 2
         * 
         *  \verbatim
           remove all the state prior to the minimum of
            current_printing_base and 
                previous_backbase     
                .                      backbase                     VCFfile->site()
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
         
        // backbase is the sequence position that the current mutation "VCFfile->site()" go back lag    
        double backbase = max(0.0, VCFfile->site() - pfARG_para.lag);
        
        // previous_backbase is the sequence position that the previous mutation "VCFfile->site()" go back lag
        // Therefore, anything before previous_backbase should not have been removed ...
        double previous_backbase = max(0.0, previous_x - pfARG_para.lag); 

        if (VCFfile->withdata() && VCFfile->site() > model->loci_length()){
            cout<<" VCF data is beyond loci length"<<endl;
            VCFfile->force_to_end_data();
            countNe->extract_and_update_count( current_states, previous_backbase, backbase);
            continue;
        }
        
        //cout<<" previous_backbase = " << previous_backbase<<", backbase = "<<backbase<<", current mutation is at base "<< VCFfile->site()<<", previous_x = "<<previous_x<<endl;
        
        valarray<double> weight_cum_sum((pfARG_para.N+1)); //Initialize the weight cumulated sums
            
        /*!     Sample the next genealogy, before the new data entry is updated to the particles 
         *      In this case, we will be update till VCFfile->site() 
         */
        current_states.update_state_to_data(VCFfile, model, weight_cum_sum);
                
        /*! WRITE TMRCA AND BL TO FILE, This is used when generating the heatmap
         */
        current_states.appendingStuffToFile( backbase, pfARG_para);    
        
        /*!     UPDATE CUM COUNT FOR WEIGHT AND BRANCH LENGTH */ 
        countNe->extract_and_update_count( current_states, previous_backbase, backbase);
        // This could return a value for which the earliest base position, things can be removed ...
        
        /*! Reset population sizes in the model */
        countNe->reset_model_Ne( model, pfARG_para.online_bool, false);
        
        /*! ESS resampling */        
        current_states.ESS_resampling(weight_cum_sum, sample_count, VCFfile->site(), pfARG_para.ESSthreshold, pfARG_para.N);
        
        /*! Clean up states that are no longer can be traced back from the present particles 
         *  
         *  As we have already updated the Ne counts, we could clean up any ForestState before backspace, 
         *  However, need to check this with the current_printing_space
         */
        current_states.clean_old_states(previous_backbase); 
        
        /*! update previous_x before move on to the next line of the data */
        previous_x = VCFfile->site();                              
        VCFfile->read_new_line(); // Read new line from the vcf file        
    
        }while(!VCFfile->end_data());
    
    cout << "### PROGRESS: end of the sequence" << endl;
    
    countNe->reset_model_Ne( model, true, true);
    
    pfARG_para.appending_Ne_file(model, true);
    current_states.clear(); // This line is sufficient to clear the memory.
    VCFfile->reset_VCF_to_data();
    
    } // End of void pfARG_core( ... )
