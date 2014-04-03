#include"vcf_data.hpp"
#include"pf-ARG.hpp"
//#include"scrm/copy_forest.h"

using namespace pfARG;


int main(int argc, char *argv[]){
    if (argc==1 ){
        pfARG::print_help();
    }    //else, proceed
    try {    
        scrm::param scrm_para(argc, argv);
        pfARG::param pfARG_para(argc, argv);
        srand(scrm_para.random_seed);
        if (scrm_para.log_NAME != pfARG_para.log_NAME) {
            remove(scrm_para.log_NAME.c_str());
            scrm_para.log_NAME=pfARG_para.log_NAME;
        };
        
//READ vcf file header, and data up to the first data entry says "PASS"
        vcf_data::vcf_header * VCFheader= new vcf_data::vcf_header(pfARG_para.vcf_NAME);
        size_t vcf_file_length=VCFheader->vcf_file_length;
        size_t previous_vcf_end_at=VCFheader->header_end_pos;
        size_t current_vcf_location = previous_vcf_end_at;
        bool not_the_end= current_vcf_location < vcf_file_length;
        scrm_para.nsites=1000000;
        scrm_para.nsam=2*VCFheader->nsam(); //overwrite the sample size according the data.
        int nsam=scrm_para.nsam; 
        scrm_para.log_param(); // rewrite scrm parameters
        
//READ the first buffer block and update the data entry to the first line !!!        
        int buff_length=pfARG_para.buff_length;    
        vcf_data::vcf_buffer * VCFBUFFER= new vcf_data::vcf_buffer(pfARG_para.vcf_NAME,previous_vcf_end_at, buff_length);
        
        previous_vcf_end_at=VCFBUFFER->end_pos;
        cout << "everything should end at " <<  VCFheader->vcf_file_length +1 <<" and " << previous_vcf_end_at<<endl;
        
        int previous_chrom=0;
        double data_base_at=0;
        
        vcf_data::vcf_line * VCFLINE;
        size_t buff_line_i=0;
        VCFLINE= new vcf_data::vcf_line(VCFBUFFER->buffer_lines[buff_line_i],nsam,previous_chrom,data_base_at );
        current_vcf_location =  current_vcf_location + VCFBUFFER->buffer_lines[buff_line_i].size()+1;
        not_the_end= current_vcf_location < vcf_file_length;
        buff_line_i++;
        
        previous_chrom=VCFLINE->chrom();
        data_base_at=VCFLINE->site();
        
        
        
//Start building the initial particles        
        time_t start_time = time(0);
               
        //Initialize c++ default random number generator 
        srand(scrm_para.random_seed);
        //Initialize mersenneTwister seed
        MersenneTwister *rg = new MersenneTwister(scrm_para.random_seed);
        Model *model = new Model(scrm_para);
        
    //First create a bigger list of initial particles
        int initial_num_particles=100*pfARG_para.N;
        Vparticle * pointer_to_current_particles=new Vparticle(model, rg, initial_num_particles, VCFLINE->vec_of_sample_alt_bool, scrm_para.theta);    
        //Vparticle * pointer_to_current_particles=new Vparticle(model, rg, 100*pfARG_para.N);    //build initial particles sample pool 100 times the number of particles needed, then rescale
        //pointer_to_current_particles->update_particle_weights_at_A_single_site(VCFLINE->vec_of_sample_alt_bool, scrm_para.theta);
        valarray<double> initial_weight_cum_sum((initial_num_particles+1)); //Initialize the weight cumulated sums        
        
        update_cum_sum_array(pointer_to_current_particles,initial_weight_cum_sum,initial_num_particles);    
        
        valarray<int> initial_sample_count(initial_num_particles);    // if sample_count is in the while loop, this is the initializing step...

        systemetic_resampling( initial_weight_cum_sum, initial_sample_count, pfARG_para.N);

        for (size_t i=0;i<initial_sample_count.size();i++){dout<<initial_sample_count[i]<<"  ";}dout<<std::endl;
        
//cout << "pointer_to_current_particles->particles.size()=" <<pointer_to_current_particles->particles.size()<<endl;
        pointer_to_current_particles=new Vparticle(pointer_to_current_particles, initial_sample_count);  
//cout << "pointer_to_current_particles->particles.size()=" <<pointer_to_current_particles->particles.size()<<endl;
//exit(1);
        double next_base=data_base_at;
        pfARG_para.appending_TMRCA(next_base, pointer_to_current_particles);

        double curly_L = pointer_to_current_particles->particles[0]->seq_posi_4_gt_change();// curly_L is the l_i of the first particle
        next_base=next_base+curly_L; //next position of genealogy to sequence

time_t initial_particle_end_time = time(0);

        size_t I=0; // Initialize the iteration
        valarray<int> sample_count(pfARG_para.N);    // if sample_count is in the while loop, this is the initializing step...
        while (not_the_end){
            
            
            valarray<double> weight_cum_sum((pfARG_para.N+1)); //Initialize the weight cumulated sums
            
            dout << "everything should end at " <<  VCFheader->vcf_file_length +1 <<" and currently is at " << current_vcf_location<<endl;
            dout<< " ******************** Update the weight of the particles " << I<<"th iteration ********** " <<endl;
            bool keep_updating_weights=true;
            while(keep_updating_weights ){
                assert(buff_line_i<VCFBUFFER->buffer_lines.size());            
                dout<<VCFBUFFER->buffer_lines.size() << " lines, the " << buff_line_i<<"th line: " <<  VCFBUFFER->buffer_lines[buff_line_i]<<endl;
                VCFLINE= new vcf_data::vcf_line(VCFBUFFER->buffer_lines[buff_line_i],nsam,previous_chrom,data_base_at );
                if (!VCFLINE->skip){                    
                    double buff_base=VCFLINE->site();
                    dout << "next change at " << next_base<<", data site: " << buff_base<<endl;
                    if ( next_base > buff_base){
                        previous_chrom=VCFLINE->chrom();
                        data_base_at=VCFLINE->site();
                        cout << "yes update weights, at " << data_base_at <<", next change is at " << next_base <<endl;
                    //dout<<VCFBUFFER->buffer_lines[i]<<endl;
                    //VCFLINE->print_vcf_line(VCFheader->sample_names);
                    
                    //particles read VCFLINE data
                        pointer_to_current_particles->update_particle_weights_at_A_single_site(VCFLINE->vec_of_sample_alt_bool, scrm_para.theta);
                                
                    }
                    else{
                        dout << " since it's ready to change" << endl;
                        keep_updating_weights=false;
                        break;
                    }
                }
                dout << "next_base is " << next_base<<" currently data is at " << data_base_at<<endl;
                if (keep_updating_weights){
                    current_vcf_location =  current_vcf_location + VCFBUFFER->buffer_lines[buff_line_i].size()+1;
                    buff_line_i++;
                    not_the_end= current_vcf_location < vcf_file_length;
                    //if (buff_line_i==VCFBUFFER->buffer_lines.size()){cout << "should end" << endl;}
                    if (buff_line_i==VCFBUFFER->buffer_lines.size() && not_the_end){ // need to reload and reset the buffer line counter
                        cout << "Reload more buff" << endl;
                        VCFBUFFER= new vcf_data::vcf_buffer(pfARG_para.vcf_NAME,previous_vcf_end_at, buff_length);
                        previous_vcf_end_at=VCFBUFFER->end_pos;
                        buff_line_i=0;
                        //buffer_is_long_enough=false;
                        
                    }
                    //else{
                        //break;
                    //}
                    if (buff_line_i==VCFBUFFER->buffer_lines.size()){break;}
                }
                else{
                    break;
                }
                
            }
            //cout<<endl<<"keep_update_weigth=" << keep_updating_weights<<endl<<"out of updating from data" << endl;
            if (!keep_updating_weights){
                update_cum_sum_array(pointer_to_current_particles,weight_cum_sum,pfARG_para.N);    
                dout<< " ******************** Update the weight of the particles END********** " <<endl;
                dout<<std::endl<<"   ------------------  Iteration # " <<  I << " END  ---------------  " << std::endl; 
                I++;    
                dout<<std::endl<<"   ------------------  Iteration # " <<  I << "   ---------------  " << std::endl;                      
                pointer_to_current_particles->particles[0]->particle_forest->set_current_base(curly_L+pointer_to_current_particles->particles[0]->particle_forest->current_base() ); // curly_L needs to be greater than 1, otherwise, it does not move
                pointer_to_current_particles->particles[0]->particle_forest->sampleNextGenealogy();
                dout << "updated new lambda " << pointer_to_current_particles->particles[0]->particle_forest->local_tree_length()<<std::endl;
                pointer_to_current_particles->particles[0]->Cal_next_position();
                
                //update the weight
                /* code */
                //using a random number first 
                //pointer_to_current_particles->particles[0]->weight=unifRand();
                //dout << "proposed new weight " << pointer_to_current_particles->particles[0]->weight<<endl;
                //Resampling 
                //sample_count=0; // Initialize sample count at the beginning of a resampling procudue, do this in the function systemetic_resampling();
                systemetic_resampling( weight_cum_sum, sample_count,pfARG_para.N);
                //trivial_resampling(pfARG_para.N,sample_count);             
                for (size_t i=0;i<sample_count.size();i++){dout<<sample_count[i]<<"  ";}dout<<std::endl;
                //updating_current_particles(pointer_to_current_particles, sample_count, weight_cum_sum,curly_L, pfARG_para.N);
                 //pointer_to_current_particles=updating_current_particles(pointer_to_current_particles, sample_count, weight_cum_sum,curly_L, pfARG_para.N);
                 pointer_to_current_particles=new Vparticle(pointer_to_current_particles, sample_count);  
                 pfARG_para.appending_TMRCA(next_base, pointer_to_current_particles);
                    //cout << "Printing the particle " <<  pointer_to_current_particles->particles.size()<<endl;
    
        //THIS IS FOR CHECKING
            //dout << "##CHECK AGAINER!!! THEY DO NOT MATCH THE PREVIOUS TABLE, LAMBDA AND WEIGHT !!! ########current particles#####################" << endl;
            //for (size_t i=0; i<pfARG_para.N;i++){
                //pointer_to_current_particles->particles[i]->particle_forest->printTree();
            //}
        //THIS IS FOR CHECKING
            
                curly_L = pointer_to_current_particles->particles[0]->seq_posi_4_gt_change();// curly_L is the l_i of the first particle
                next_base=next_base+curly_L;
                cout << "new next genealogy change at " << next_base<<endl;
                assert(pointer_to_current_particles->particles.size() == pfARG_para.N);
            }
        }
        
        time_t end_time = time(0);
        
        if (pfARG_para.log_bool){  
            std::ofstream log_file;
            log_file.open (pfARG_para.log_NAME.c_str(), std::ios::out | std::ios::app | std::ios::binary); 
            log_file << "Initial particles building took about " << initial_particle_end_time - start_time << " second(s) \n";
            log_file << "Simulation took about " << end_time - initial_particle_end_time << " second(s) \n";
            log_file.close();
            string log_cmd="cat "+pfARG_para.log_NAME;
            system(log_cmd.c_str());
        }
        
        //delete pointer_to_current_particles;
        delete rg;
        delete model;
        delete VCFheader;
    }
    catch (const exception &e)
    {
      std::cerr << "Error: " << e.what() << std::endl;
      return EXIT_FAILURE;
    }
}


void cal_marginal_likelihood_old(Node * node, double mutation_rate){// Genealogy branch lengths are in number of generations, the mutation rate is unit of per site per generation, often in the magnitute of 10 to the power of negative 8.
    
    if (node->first_child() == NULL){
        if (node->mutation_state){
            node->marginal_likelihood[1]=1.0;
            node->marginal_likelihood[0]=0.0;    
        }
        else{
            node->marginal_likelihood[1]=0.0;
            node->marginal_likelihood[0]=1.0;        
        }
        //dout << "Marginal probability at " << node->label() << " is " << node->marginal_likelihood[0]<<"," << node->marginal_likelihood[1]<<endl;

    }
    else{
        double t1=node->height()- node->first_child()->height();
        double ut1=exp(-t1*mutation_rate); // let ut1 be the probability of no mutation on the branch to the first child
assert(ut1>=0 && ut1<=1);

        double t2=node->height()- node->second_child()->height();
        double ut2=exp(-t2*mutation_rate); // let ut2 be the probability of no mutation on the branch to the second child
assert(ut2>=0 && ut2<=1);
        cal_marginal_likelihood(node->first_child(),mutation_rate);
        double y[2];
        y[0]= node->first_child()-> marginal_likelihood[0];
        y[1]= node->first_child()-> marginal_likelihood[1];
        cal_marginal_likelihood(node->second_child(),mutation_rate);
        double z[2];
        z[0]= node->second_child()-> marginal_likelihood[0];
        z[1]= node->second_child()-> marginal_likelihood[1];
        
        //node->marginal_likelihood[0] = y[0] * ut1 * z[0] * ut2 + y[0] * ut1 * z[1] * (1-ut2) + y[1] * (1-ut1) * z[0] * ut2 + y[1] * (1-ut1) * z[1] * (1-ut2);
        //node->marginal_likelihood[1] = y[1] * ut1 * z[1] * ut2 + y[1] * ut1 * z[0] * (1-ut2) + y[0] * (1-ut1) * z[1] * ut2 + y[0] * (1-ut1) * z[0] * (1-ut2);
        node->marginal_likelihood[0] = (y[0]*ut1 + y[1]*(1-ut1)) * (z[0]*ut2 + z[1]*(1-ut2)) ;
        node->marginal_likelihood[1] = (y[1]*ut1 + y[0]*(1-ut1)) * (z[1]*ut2 + z[0]*(1-ut2)) ;
        
        //dout << "t1 = " <<t1<<", and t2= " << t2<<endl;
        //dout << "ut1 = " <<ut1<<", and ut2= " << ut2<<endl;
        //dout << "Marginal probability at " << node << " is " << node->marginal_likelihood[0]<<"," << node->marginal_likelihood[1]<<endl;

    }    
}

