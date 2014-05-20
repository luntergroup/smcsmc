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

#include"count.hpp" 
//#include"debug/usage.hpp"

/*!
 * Global variables for debugging the number of times that ForestState was constructed and removed.
 */ 
int new_forest_counter = 0;
int delete_forest_counter = 0;

void pfARG_core(PfParam &pfARG_para,
                CountModel *countNe,
                bool print_update_count);

int main(int argc, char *argv[]){   
    
    /*! 
     * Default values
     */ 
    bool print_update_count = false;

    /*! Extract pfARG parameters */
    PfParam pfARG_para( argc, argv );
    
    if ( argc==1 ){
        //pfARG_para.print_help();
        } 
    try {//else, proceed
        
        /*! 
         * INITIALIZE CountModel
         */
        CountModel *countNe = new CountModel( *pfARG_para.model , pfARG_para.lag);
        pfARG_para.appending_Ne_file( true ); // Append initial values to History file
        
        /*!
         * EM step
         */                 
        for (int I = 0; I <= pfARG_para.EM_steps; I++){
            cout << "Now starting EM_step " << I << endl;
            pfARG_core( pfARG_para, 
                        countNe, 
                        print_update_count);
            cout << "End of EM_step " << I << endl;
            }

        pfARG_para.appending_Ne_file( );
        
        int exit_success = pfARG_para.log( );
        
        /*! Clean up */
        delete countNe;
        //cout<<"Forest state was created " << new_forest_counter << " times" << endl;
        //dout<<"Forest state destructor was called " << delete_forest_counter << " times" << endl;
        
        return exit_success;
        } 
    catch (const exception &e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return EXIT_FAILURE;
        }
    }


void pfARG_core(PfParam &pfARG_para,
                CountModel *countNe,
                bool print_update_count){
                    
    //int who= RUSAGE_SELF;
    //struct rusage usage;
    //struct rusage *p=&usage;
    //int ret;

    Model *model = pfARG_para.model;
    MersenneTwister *rg = pfARG_para.rg;
    size_t Nparticles = pfARG_para.N ;
    Vcf *VCFfile = pfARG_para.VCFfile;
    
    VCFfile->read_new_line();
    cout<< "############# starting vcf file at base " <<VCFfile->site()<<endl;
    /*! Initial particles */ 

    double initial_position = 0;

    ParticleContainer current_states(model, rg, Nparticles, VCFfile->vec_of_sample_alt_bool, VCFfile->withdata(), initial_position );             
    dout<<"######### finished initial particle building"<<endl;
    valarray<int> sample_count( Nparticles ); // if sample_count is in the while loop, this is the initializing step...

    /*! Initialize prior Ne */
    countNe->init();

    /*! Go through vcf data */
    double previous_x = 0;

    do{
        //ret=getrusage(who,p);
        //process(p, VCFfile->site());

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
         * backbase is the sequence position that the current mutation "VCFfile->site()" go back lag    
         * 
         * previous_backbase is the sequence position that the previous mutation "VCFfile->site()" go back lag
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
         
        
        if (VCFfile->withdata() && VCFfile->site() > model->loci_length()){
            cout<<" VCF data is beyond loci length"<<endl;
            VCFfile->force_to_end_data();
            countNe->extract_and_update_count( current_states , VCFfile->site() );

            continue;
            }
        
        valarray<double> weight_cum_sum((Nparticles+1)); //Initialize the weight cumulated sums
            
        /*!     Sample the next genealogy, before the new data entry is updated to the particles 
         *      In this case, we will be update till VCFfile->site() 
         */
        current_states.update_state_to_data(VCFfile, model, weight_cum_sum, pfARG_para.mid_bool, pfARG_para.finite_bool);
                
        /*! WRITE TMRCA AND BL TO FILE, This is used when generating the heatmap
         */
         
        current_states.appendingStuffToFile( VCFfile->site(), pfARG_para);    

        /*!     UPDATE CUM COUNT FOR WEIGHT AND BRANCH LENGTH */ 
        countNe->extract_and_update_count( current_states , VCFfile->site() );
        // This could return a value for which the earliest base position, things can be removed ...
        
        /*! Reset population sizes in the model */
        countNe->reset_model_parameters( model, pfARG_para.online_bool, false);

        if ( pfARG_para.ESS() == 1 ){
            cout << " random weights" <<endl;
            current_states.set_particles_with_random_weight();    
            }
        
        
        /*! ESS resampling */        
        current_states.ESS_resampling(weight_cum_sum, sample_count, VCFfile->site(), pfARG_para.ESSthreshold, Nparticles);
        
        /*! Clean up states that are no longer can be traced back from the present particles 
         *  
         *  As we have already updated the Ne counts, we could clean up any ForestState before backspace, 
         *  However, need to check this with the current_printing_space
         */
        
        /*! update previous_x before move on to the next line of the data */
        previous_x = VCFfile->site();                              
        VCFfile->read_new_line(); // Read new line from the vcf file        
    
        }while(!VCFfile->end_data());
    
    //ret=getrusage(who,p);
    //process(p, VCFfile->site());

    
    cout <<endl << "### PROGRESS: end of the sequence" << endl;

    countNe->extract_and_update_count( current_states , previous_x, true );
    countNe->reset_model_parameters( model, true, true); // This is mandatory for appending the correct value out ...
    
    pfARG_para.appending_Ne_file( true );

    current_states.clear(); // This line is sufficient to clear the memory.
    VCFfile->reset_VCF_to_data();
    
    } // End of void pfARG_core( ... )
