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


#include"pfparam.hpp"
#include"pattern.hpp"


PfParam::PfParam(int argc, char *argv[]): argc_(argc), argv_(argv) {
    this->init(); // Initialize pfARG program parameters
    argc_i=1; // Skipping argv[0], which is pf-ARG calling
    
    while( argc_i < argc ){
        
        string argv_i(argv_[argc_i]);

        // ------------------------------------------------------------------
        // Parameters 
        // ------------------------------------------------------------------        
        if ( argv_i == "-Np" ){ 
            nextArg( argv_i );
            N = readInput<size_t>(argv_[argc_i]);
            }

        else if ( argv_i == "-nsam" ){ 
            nextArg( argv_i );
            default_nsam = readInput<size_t>(argv_[argc_i]);
            }
                        
        
        else if (argv_i == "-ESS"){ 
            nextArg( argv_i );
            ESS = readInput<double>(argv_[argc_i]);
            ESS_default_bool = false;
            }
                
        else if (argv_i == "-EM"){
            this->EM_bool = true;
            nextArg( argv_i );
            EM_steps = readInput<int>(argv_[argc_i]);
            }        
        
        else if (argv_i == "-tmax"){
            nextArg( argv_i );
            top_t = readInput<double>(argv_[argc_i]);
            }
        
        else if (argv_i == "-p"){ 
            nextArg( argv_i );
            pattern = argv_[argc_i];
            }
            
        // ------------------------------------------------------------------
        // Input 
        // ------------------------------------------------------------------
        else if (argv_i == "-vcf"){ 
            nextArg( argv_i );
            vcf_NAME = argv_[argc_i];
            }
        
        else if (argv_i == "-buff"){ 
            nextArg( argv_i );
            buff_length = readInput<int>(argv_[argc_i]);
            }    
        
        // ------------------------------------------------------------------
        // Action 
        // ------------------------------------------------------------------
        else if (argv_i == "-lag"){ 
            nextArg( argv_i );
            lag = readInput<double>(argv_[argc_i]);
            }

        else if (argv_i == "-online"){
            this->online_bool = true;
            }
            
        // ------------------------------------------------------------------
        // Output 
        // ------------------------------------------------------------------
        else if (argv_i == "-o"){ 
            nextArg( argv_i );
            out_NAME_prefix = argv_[argc_i];
            }
        
        else if ( argv_i == "-log"  ){ this->log_bool  = true; }        
        else if ( argv_i == "-heat" ){ this->heat_bool = true; }        
        else if ( argv_i == "-hist" ){ this->hist_bool = true; }
        
        else if (argv_i == "-h" || argv_i == "-help") {
            this->print_option();
            exit(1);
            }
            
        else {
            scrm_input += argv_i + " ";
            }
        
        argc_i++;
        }
        
        this->finalize( );
    }


PfParam::~PfParam(){ 
    //cout<<"~PfParam() is called"<<endl;
    delete this->VCFfile; 
    delete this->model;
    delete this->SCRMparam;
    }
    
    
//void PfParam::get_scrm_argv(){
    //scrm_input += convert_pattern(pattern, top_t);
    //cout << scrm_input << endl;
    //}
    

void PfParam::nextArg(std::string option) {
    ++argc_i;
    if (argc_i >= argc_) {
        throw std::invalid_argument(std::string("Not enough parameters when parsing option ") + option);
        }
    }

/*! 
 * Set default parameters
 */  
void PfParam::init(){
    this->default_nsam        = 2;
    this->default_mut_rate    = 0.00000001;
    this->default_recomb_rate = 0.000000001;
    this->default_loci_length = 5000000;    
    
    this->N                = 100;
    this->buff_length      = 200;
    this->lag              = 0;
    this->out_NAME_prefix  = "pfARG";
    this->ESS              = 0.5;
    this->ESS_default_bool = true;
    this->log_bool         = false;
    this->heat_bool        = false;
    this->hist_bool        = false;
    this->online_bool      = false;
    this->window           = 100; 
    this->EM_steps         = 0;
    this->EM_bool          = false;
    
    this->VCFfile          = NULL;
    this->SCRMparam        = NULL;
    this->model = new Model();
    this->scrm_input       = "";
    this->top_t            = 2;
}


void PfParam::insert_mutation_rate_in_scrm_input ( ) {
    size_t found = scrm_input.find("-t");
    if ( found == std::string::npos ) { // if "-t" option is not find ... 
        this->scrm_input = "-t " + to_string ( this->default_mut_rate * this->default_loci_length * 4 * 10000 ) + " " + this->scrm_input;
            //this->scrm_input = to_string ( nsam ) + " 1 " + this->scrm_input; 
        }
    //else { // if "-t" is find ... 
        
        //}
    
    }

    
void PfParam::insert_recomb_rate_and_seqlen_in_scrm_input (  ){
    size_t found = scrm_input.find("-r");
    if ( found == std::string::npos ) { // if "-r" option is not find ... 
        this->scrm_input = "-r " + to_string ( this->default_recomb_rate * this->default_loci_length * 4 * 10000 )  // number of recombination events in 4N0, as the scaling N0 in screen is 10000
                                 + " " + to_string ((size_t)this->default_loci_length) + " " + this->scrm_input;            // sequence length
        }
    else {
        // skipping the number of recombination first
        size_t pos_start = scrm_input.find(" ", found+2, 1);
        size_t pos_end = scrm_input.find(" ", pos_start+1, 1);        
        // extract the loci length
        pos_start = scrm_input.find(" ", pos_start+2, 1);
        pos_end = scrm_input.find(" ", pos_start+1, 1);        
        this->default_loci_length = std::strtod ( (char*)scrm_input.substr(pos_start, pos_end - pos_start).c_str(), NULL);
        }

    if (!VCFfile->withdata()){
            VCFfile->set_even_interval( this->default_loci_length / 10 );
        }            

    }

    
void PfParam::insert_sample_size_in_scrm_input (  ){
    size_t nsam = 2*VCFfile->nsam(); // Extract number of samples from VCF file
    if (!VCFfile->withdata()){
        nsam = this->default_nsam;
        }    
    this->scrm_input = to_string ( nsam ) + " 1 " + this->scrm_input; 
    }


void PfParam::finalize_scrm_input (  ){
    // These options insert parameters to the beginning of the current scrm_input
    this->insert_recomb_rate_and_seqlen_in_scrm_input ( );
    this->insert_mutation_rate_in_scrm_input ( );
    this->insert_sample_size_in_scrm_input ( );       
    
    this->scrm_input = "scrm " + this->scrm_input + convert_pattern(pattern, top_t);    
    cout << scrm_input <<endl;
    this->convert_scrm_input ();
    }


void PfParam::convert_scrm_input (){
    enum { kMaxArgs = 264 };
    int scrm_argc = 0;
    char *scrm_argv[kMaxArgs];        
    char * p2 = strtok((char *)this->scrm_input.c_str(), " ");
    while (p2 && scrm_argc < kMaxArgs) {
        scrm_argv[scrm_argc++] = p2;
        p2 = strtok(0, " ");
        }
    
    ///*! Extract scrm parameters */ 
    this->SCRMparam = new Param(scrm_argc, scrm_argv, false);
    this->SCRMparam->parse( *this->model );
    }


void PfParam::finalize(  ){
     /*! Initialize vcf file, and data up to the first data entry says "PASS"   */
    this->VCFfile =  new Vcf(this->vcf_NAME, this->buff_length);

    this->finalize_scrm_input ( );
    //this->lag = this->default_loci_length / 20; // TESTING, just to try use full lagging...
    this->ESSthreshold = this->N * this->ESS;        
 
    this->TMRCA_NAME   = out_NAME_prefix + "TMRCA";
    this->WEIGHT_NAME  = out_NAME_prefix + "WEIGHT";
    this->BL_NAME      = out_NAME_prefix + "BL";
    this->Ne_NAME      = out_NAME_prefix + "Ne";
    this->log_NAME     = out_NAME_prefix + ".log";
    this->HIST_NAME    = out_NAME_prefix + "HIST";
    this->SURVIVOR_NAME= out_NAME_prefix + "SURVIVOR";
    
    remove( this->TMRCA_NAME.c_str() );
    remove( this->WEIGHT_NAME.c_str());
    remove( this->BL_NAME.c_str()    );
    remove( this->Ne_NAME.c_str()    );
    remove( this->log_NAME.c_str()   );
    remove( this->HIST_NAME.c_str()  );
    remove( this->SURVIVOR_NAME.c_str());
    }


//int PfParam::log(Model *model, size_t random_seed, pfTime * runningtime, double inferred_recomb_rate){
//int PfParam::log(Model *model, size_t random_seed, double inferred_recomb_rate){
int PfParam::log( double inferred_recomb_rate ){
    if (log_bool){  
        this->log_param(inferred_recomb_rate);
        //log_end(runningtime);
        string log_cmd="cat " + log_NAME;
        return system(log_cmd.c_str());  
    } else {
        return 0;    
    }
}


//void PfParam::log_param(Model *model, size_t random_seed, double inferred_recomb_rate){
void PfParam::log_param( double inferred_recomb_rate){
    ofstream log_file;
    string emptyfile("EMPTY FILE");
    string vcf_file = ( vcf_NAME.size() > 0 ) ?  vcf_NAME : emptyfile;
    
    log_file.open (log_NAME.c_str(), ios::out | ios::app | ios::binary); 
    
    log_file << "pf-ARG parameters: \n";
    if (this->heat_bool){
    log_file << "TMRCA saved in file: "  << TMRCA_NAME  << "\n";
    log_file << "WEIGHT saved in file: " << WEIGHT_NAME << "\n";
    log_file << "BL saved in file: "     << BL_NAME     << "\n";
    }
    if (this->hist_bool){
    log_file << "HIST saved in file: "   << HIST_NAME   << "\n";
    }
    log_file << "Ne saved in file: "     << Ne_NAME     << "\n";        
    log_file << "VCF data file: "        << vcf_file    <<"\n";
    log_file << setw(15) <<     " EM steps =" << setw(10) << EM_steps                    << "\n";
    log_file << setw(15) <<          " lag =" << setw(10) << lag                         << "\n";
    log_file << setw(15) <<             "N =" << setw(10) << N                           << "\n";
    log_file << setw(15) <<           "ESS =" << setw(10) << ESS; 
    if (ESS_default_bool){ log_file << " (by default)";}                        log_file << "\n";
    log_file << setw(15) <<        "buffer =" << setw(10) << buff_length                 << "\n";
    
    log_file<<"scrm model parameters: \n";
    log_file << setw(15) <<   "Sample size =" << setw(10) << this->model->sample_size()        << "\n";
    log_file << setw(15) << "mutation rate =" << setw(10) << this->model->mutation_rate()      << "\n";
    log_file << setw(15) <<   "recomb rate =" << setw(10) << this->model->recombination_rate() << "\n";
    log_file << setw(15) <<"inferred recomb rate = " << setw(10) << inferred_recomb_rate       << "\n";
    log_file << setw(15) <<    "Seq length =" << setw(10) << this->model->loci_length()        << "\n";
    log_file << setw(15) <<"Extract window =" << setw(10) << this->model->exact_window_length()<< "\n";
    log_file << setw(15) <<         "seed : " << setw(10) << this->SCRMparam->random_seed      << "\n";

    this->model->resetTime();
    log_file<<setw(15)<<"Pop size (at Generation):\n";
    for (size_t i = 0; i < this->model->change_times_.size()-1; i++){
        log_file<<setw(3)<<"(" << setw(8) << this->model->getCurrentTime() <<" )"; 
        for (size_t pop_j = 0 ; pop_j < this->model->population_number() ; pop_j++){
            log_file << " | " << setw(10)<<this->model->population_size(pop_j);
        }
        
        log_file<< "\n";
        this->model->increaseTime();
    }
    log_file<<setw(3)<<"(" << setw(8) << this->model->getCurrentTime() <<" )" ;
    for (size_t pop_j = 0 ; pop_j < this->model->population_number() ; pop_j++){
        log_file << " | " << setw(10)<<this->model->population_size(pop_j);
    }
    log_file<< "\n";
    
    this->model->resetTime();
    log_file << setw(15) << "Migration rate :" << "\n";
    for (size_t pop_i = 0 ; pop_i < this->model->population_number() ; pop_i++){
        log_file << setw(15) << " ";
        for (size_t pop_j = 0 ; pop_j < this->model->population_number() ; pop_j++){
            log_file << setw(10) << this->model->migration_rate(pop_i, pop_j)  ;
        }
        log_file << "\n";
    }
    
    log_file.close();
}        

//void PfParam::log_end(pfTime * running_time){
    //ofstream log_file ( log_NAME.c_str(), ios::out | ios::app | ios::binary); 
    //log_file << "Initial particles building took about " << running_time->timing_[0] << " second(s) \n";
    ////log_file << "Simulation took about " << end_time - initial_state_end_time << " second(s) \n";
    //log_file << "    Update weight took " << running_time->timing_[2]<<" second(s)\n";
    //log_file << "    Resampling took " << running_time->timing_[1]<<" second(s)\n";
    //log_file.close();
    //}


void PfParam::appending_Ne_file(Model *model, bool hist){
    string file_name = hist ? HIST_NAME : Ne_NAME ;
    ofstream Ne_file( file_name.c_str(), ios::out | ios::app | ios::binary);   
    model->resetTime();
    for (size_t i = 0; i < model->change_times_.size()-1; i++){
        Ne_file << model->getCurrentTime() / model->default_pop_size / 4 << "\t" << model->population_size() / model->default_pop_size << "\n" ;
        model->increaseTime();
        }
    Ne_file << model->getCurrentTime() / model->default_pop_size / 4 << "\t" << model->population_size() / model->default_pop_size << "\n" ;
    Ne_file.close();
    return;
    }        


void PfParam::print_help(){
    cout << endl;
    cout << endl;
    cout << "*****************************************************************" << endl;
    cout << "*                          pf-ARG                               *" << endl;
    cout << "*              Author:  Sha Zhu, Gerton Lunter                  *" << endl;
    cout << "*****************************************************************" << endl;
    cout << endl<<endl;
    cout << "Too few command line arguments" << endl;
    cout << "    Options:" << endl;
    this->print_option();
    this->print_example();
    exit(1);
    }


void PfParam::print_option(){
    cout <<"   Options:" << endl;
    cout << setw(10)<<"-Np"     << setw(5) << "INT" << "  --  " << "Number of particles [ 1000 ]." << endl;
    cout << setw(10)<<"-ESS"    << setw(5) << "FLT" << "  --  " << "Proportion of the effective sample size [ 0.6 ]." << endl;
    cout << setw(10)<<"-p"      << setw(5) << "STR" << "  --  " << "Pattern of time segment [ \"3*1+2*3+4\" ]." <<endl;
    cout << setw(10)<<"-tmax"   << setw(5) << "FLT" << "  --  " << "Maximal time, in unit of 4N0 [ 3 ]." <<endl;
    cout << setw(10)<<"-EM"     << setw(5) << "INT" << "  --  " << "EM steps [ 20 ]." << endl;
    cout << setw(10)<<"-lag"    << setw(5) << "FLT" << "  --  " << "Lagging step [ 1000 ]." << endl;
    cout << setw(10)<<"-vcf"    << setw(5) << "STR" << "  --  " << "Data file in vcf format [ Chrom1.vcf ]." << endl;
    //cout << setw(20)<<"-buff BUFFSIZE" << "  --  " << "User define the size of buffer for the vcf file BUFFSIZE." << endl;
    cout << setw(10)<<"-o"      << setw(5) << "STR" << "  --  " << "Prefix for output files" << endl;
    cout << setw(10)<<"-online" << setw(5) << " "   << "  --  " << "Perform online EM" << endl;
    cout << setw(10)<<"-log"    << setw(5) << " "   << "  --  " << "Generate *.log file" << endl;
    cout << setw(10)<<"-heat"   << setw(5) << " "   << "  --  " << "Generate *TMRCA and *WEIGHT for heatmap" << endl;
    //cout << setw(20)<<"-TMRCA" << "  --  " << "User define the output file for the time to the most recent common ancester, if it is not speciefied, the file name is TMRCA by default." << endl;   
    }


void PfParam::print_example(){
    cout << "Example:" << endl;
    cout << "pf-ARG 10 -nsam 3" << endl;
    cout << "./pf-ARG -Np 5 -t 0.002 -r 400 -npop 20000 -vcf eg_vcf.vcf -buff 4" << endl;
    cout << "./pf-ARG -Np 5 -t 0.002 -r 400 -npop 20000 -vcf eg_vcf.vcf -log -TMRCA myTMRCA" << endl;
    cout << "./pf-ARG -Np 5 -t 0.002 -r 400 -npop 20000 -vcf eg_vcf.vcf" << endl;
    cout << "./pf-ARG -Np 6 -t 0.0002 -r 30 -npop 10000 -seed 1314 -vcf eg_vcf.vcf" << endl;
    cout << "./pf-ARG -Np 7 -t 0.002 -log -r 400 -vcf eg_vcf.vcf " << endl;
    cout << "./pf-ARG -Np 8 -t 0.002 -r 400 -log LOGFILE -vcf eg_vcf.vcf" << endl;
    }

