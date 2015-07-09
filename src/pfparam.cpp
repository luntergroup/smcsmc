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


#include "pfparam.hpp"
#include "pattern.hpp"
#include "rescue.hpp"
#include "help.hpp"

PfParam::PfParam(int argc, char *argv[]): argc_(argc), argv_(argv) {
    this->init(); // Initialize pfARG program parameters
    argc_i=1; // Skipping argv[0], which is smcsmc calling
    
    while( argc_i < argc ){
        
        string argv_i(argv_[argc_i]);

        // ------------------------------------------------------------------
        // Parameters 
        // ------------------------------------------------------------------        
        if      ( argv_i == "-Np"   ){ this->N = readNextInput<size_t>(); }
        else if ( argv_i == "-nsam" ){ this->default_nsam = readNextInput<size_t>(); }
        else if ( argv_i == "-ESS"  ){ this->ESS_ = readNextInput<double>();
                                       this->ESS_default_bool = false; }
        else if ( argv_i == "-EM"   ){ this->EM_steps = readNextInput<int>();
                                       this->EM_bool = true; }
        else if ( argv_i == "-xr" || argv_i == "-xc" ) 
	{
            int last_epoch;
            int first_epoch = readRange(last_epoch);        // obtain 1-based closed interval
	    first_epoch--;                                  // turn into 0-based half-open
            for (int i=0; i<last_epoch; i++) {
                if (record_event_in_epoch.size() <= i) {
		    // extend vector, and set default: record both recomb and coal/migr events
		    record_event_in_epoch.push_back( PfParam::RECORD_COALMIGR_EVENT | PfParam::RECORD_RECOMB_EVENT );
		}
		if (i >= first_epoch) {
		    // reset bit signifying recording of either recombination or coal/migr events
		    // " &= " is and-update (cf. +=, sum-update); " ~ " is bitwise not
		    record_event_in_epoch[i] &= ~( (argv_i == "-xc") ? PfParam::RECORD_COALMIGR_EVENT : PfParam::RECORD_RECOMB_EVENT );
		}
	    }
	}
        else if ( argv_i == "-tmax" ){ this->top_t_ = readNextInput<double>(); }        
        else if ( argv_i == "-p"    ){ this->nextArg();
                                       this->pattern = argv_[argc_i]; }
            
        // ------------------------------------------------------------------
        // Input files
        // ------------------------------------------------------------------
        else if ( argv_i == "-seg"  ){ this->nextArg(); 
                                       this->input_SegmentDataFileName = argv_[argc_i]; 
                                       }

        //else if ( argv_i == "-buff" ){ this->buff_length = this->readNextInput<int>(); }    
        //else if ( argv_i == "-ghost"){ this->ghost = readNextInput<int>(); }          
        
        // ------------------------------------------------------------------
        // Action 
        // ------------------------------------------------------------------
        else if ( argv_i == "-lag"    ){ this->lag = this->readNextInput<double>(); }
        //else if ( argv_i == "-filter" ){ this->filter_window_ = this->readNextInput<int>(); }
        //else if ( argv_i == "-missing"){ this->missing_data_threshold_ = this->readNextInput<int>(); }
        else if ( argv_i == "-online" ){ this->online_bool = true; }
        else if ( argv_i == "-rescue" ){ this->rescue_bool = true; }
            
        // ------------------------------------------------------------------
        // Output 
        // ------------------------------------------------------------------
        else if ( argv_i == "-o"     ){ this->nextArg();
                                        this->out_NAME_prefix = argv_[argc_i]; }        
        else if ( argv_i == "-log"   ){ this->log_bool  = true; }
        else if ( argv_i == "-heat"  ){ this->heat_bool = true; }
        //else if ( argv_i == "-hist" ){ this->hist_bool = true; }
        //else if ( argv_i == "-finite"  ){ this->finite_bool  = true; }        
        
        else if (argv_i == "-h" || argv_i == "-help") {
            Help_header();            
        }

        else if (argv_i == "-v") {
            Help_version(this->compileTime, this->smcsmcVersion, this->scrmVersion);
            exit(0);
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
    delete this->Segfile;
    delete this->model;
    delete this->SCRMparam;
    this->rg->clearFastFunc();
    delete this->rg;
    }


/*! 
 * Set default parameters
 */  
void PfParam::init(){
    this->default_nsam        = 2;
    this->default_mut_rate    = 1e-8;
    this->default_recomb_rate = 1e-9;
    this->default_loci_length = 2e7;
    this->default_num_mut = this->default_mut_rate*40000*this->default_loci_length;
    //this->ghost = 10;

    #ifdef COMPILEDATE
        compileTime = COMPILEDATE;
    #else
        compileTime = "";
    #endif

    #ifdef SMCSMCVERSION
        smcsmcVersion = SMCSMCVERSION;
    #else
        smcsmcVersion = "";
    #endif

    #ifdef SCRMVERSION
        scrmVersion = SCRMVERSION;
    #else
        scrmVersion = "";
    #endif

    this->original_recombination_rate_ = 0;
    this->N                = 100;
    //this->buff_length      = 200;
    this->lag              = 0;
    //this->lag              = 5000000;
    this->out_NAME_prefix  = "pfARG";
    this->ESS_              = 0.5;
    this->ESS_default_bool = true;
    this->log_bool         = true; // Enable log by default
    this->heat_bool        = false;
    //this->hist_bool        = false;
    this->online_bool      = false;
    //this->finite_bool      = false;
    this->heat_seq_window  = 400; 
    //this->heat_seq_window  = 100; 
    this->EM_steps         = 0;
    this->EM_bool          = false;
    
    //this->FileType         = EMPTY;
    this->Segfile          = NULL;
    this->SCRMparam        = NULL;
    this->model = new Model();
    this->rg               = NULL;  
    this->scrm_input       = "";
    this->top_t_            = 2;
    //this->filter_window_   = 0;
    //this->filter_window_   = 2;
    //this->missing_data_threshold_ = INT_MAX;
    this->rescue_bool = false;
    }


void PfParam::insert_mutation_rate_in_scrm_input ( ) {
    size_t found = scrm_input.find("-t");
    if ( found == std::string::npos ) { // if "-t" option is not find ... 
        this->default_num_mut = this->default_mut_rate * this->default_loci_length * 4 * 10000;
        this->scrm_input = "-t " + to_string ( this->default_num_mut ) + " " + this->scrm_input;
            //this->scrm_input = to_string ( nsam ) + " 1 " + this->scrm_input; 
        }    
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

    //if ( this->FileType == EMPTY ){
        //VCFfile->ghost_num_mut = this->ghost;
        //VCFfile->set_even_interval( this->default_loci_length / VCFfile->ghost_num_mut );
        //}            

    }

    
void PfParam::insert_sample_size_in_scrm_input (  ){
    //size_t nsam = ( this->FileType == EMPTY ) ? this->default_nsam: 2*Segfile->nsam(); // Extract number of samples from Segment file
    this->scrm_input = to_string ( this->default_nsam ) + " 1 " + this->scrm_input; 
    }


void PfParam::finalize_scrm_input (  ){
    // These options insert parameters to the beginning of the current scrm_input
    this->insert_recomb_rate_and_seqlen_in_scrm_input ( );
    this->insert_mutation_rate_in_scrm_input ( );
    this->insert_sample_size_in_scrm_input ( );       
    if (pattern.size() > 0) {
		Pattern tmp_pattern( pattern, top_t_ );
		this->scrm_input = "scrm " + this->scrm_input + tmp_pattern.pattern_str;
	} else {
		this->scrm_input = "scrm " + this->scrm_input;
	}
    cout << scrm_input << endl;
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
    // By default, use SMC' model, check if exact window length has been set or not, if it not (i.e. exact window length is infinity, which is denoted by -1) set it to 0.
    if ( this->model->exact_window_length() == -1 ){
        this->model->set_exact_window_length ( 0 );
        }
    this->rg = new MersenneTwister(this->SCRMparam->random_seed());  /*! Initialize mersenneTwister seed */
    this->original_recombination_rate_ = model->recombination_rate();
    }


void PfParam::finalize(  ){
    
    this->ESSthreshold = this->N * this->ESS();
    this->TMRCA_NAME   = out_NAME_prefix + "TMRCA";
    this->WEIGHT_NAME  = out_NAME_prefix + "WEIGHT";
    //this->BL_NAME      = out_NAME_prefix + "BL";
    this->Ne_NAME      = out_NAME_prefix + "Ne";
    this->Count_NAME   = out_NAME_prefix + "Count";
    this->log_NAME     = out_NAME_prefix + ".log";
    this->HIST_NAME    = out_NAME_prefix + "HIST";
    this->SURVIVOR_NAME= out_NAME_prefix + "SURVIVOR";
    
    remove( this->TMRCA_NAME.c_str() );
    remove( this->WEIGHT_NAME.c_str());
    //remove( this->BL_NAME.c_str()    );
    remove( this->Ne_NAME.c_str()    );
    remove( this->Count_NAME.c_str() );
    remove( this->log_NAME.c_str()   );
    remove( this->SURVIVOR_NAME.c_str()); 
    if ( this->rescue_bool ){ // By default, no rescue
        cout << " Rescue from " << this->HIST_NAME.c_str() << endl;
        //RescueHist rescueHist( this->HIST_NAME );
        //this->scrm_input = rescueHist.rescured_param_string;
        cout << this->scrm_input <<endl;
        }
    else { 
        remove( this->HIST_NAME.c_str() );
        }
    
    this->finalize_scrm_input ( );
    
    // if necessary, extend the vector specifying what epochs to collect events for,
    // and check it hasn't been made too large (which wouldn't strictly be a problem,
    // but clearly the user specified something that's silly and should hear that.)
    while (record_event_in_epoch.size() < this->model->change_times_.size()) {
        record_event_in_epoch.push_back( PfParam::RECORD_COALMIGR_EVENT | PfParam::RECORD_RECOMB_EVENT );
    }
    if (record_event_in_epoch.size() > this->model->change_times_.size()) {
        throw std::invalid_argument(std::string("Problem: epochs specified in -xr/-xc options out of range"));
    }

     /*! Initialize seg file, and data up to the first data entry says "PASS"   */
    this->Segfile = new Segment( this->input_SegmentDataFileName, this->default_nsam, (double)this->model->loci_length(), this->default_num_mut );
    //this->VCFfile->filter_window_ = this->filter_window_;
    //this->VCFfile->missing_data_threshold_ = this->missing_data_threshold_;
}
        

int PfParam::log( ){
    if (log_bool){  
        this->log_param( );
        string log_cmd="cat " + log_NAME;
        return system(log_cmd.c_str());  
        } 
    else {
        return 0;    
        }
    }


void PfParam::log_param( ){
    ofstream log_file;
    log_file.open (log_NAME.c_str(), ios::out | ios::app | ios::binary); 

    log_file << "###########################\n";
    log_file << "#        smcsmc log       #\n";
    log_file << "###########################\n";
    Help_version(this->compileTime, this->smcsmcVersion, this->scrmVersion, log_file);
    log_file << "smcsmc parameters: \n";
    if (this->heat_bool){
        log_file << "TMRCA saved in file: "  << TMRCA_NAME  << "\n";
        log_file << "WEIGHT saved in file: " << WEIGHT_NAME << "\n";
        //log_file << "BL saved in file: "     << BL_NAME     << "\n";
        }
    //if (this->hist_bool){
        log_file << "HIST saved in file: "   << HIST_NAME   << "\n";
        //}
    
    log_file << "Ne saved in file: "     << Ne_NAME     << "\n";        
    log_file << "Segment Data file: " ;
    log_file << (( this->input_SegmentDataFileName.size() == 0 )? "empty" : input_SegmentDataFileName.c_str() ) << "\n";
    log_file << setw(15) <<     " EM steps =" << setw(10) << EM_steps                    << "\n";
    if (lag > 0){
        log_file << setw(15) << "Constant lag =" << setw(10) << lag                      << "\n";
        }
    if (online_bool){
        log_file << setw(15) << "Online update = TRUE\n";
    }
    log_file << setw(15) <<             "N =" << setw(10) << N                           << "\n";
    log_file << setw(15) <<           "ESS =" << setw(10) << ESS_; 
    if (ESS_default_bool){ log_file << " (by default)";}                        log_file << "\n";
    //log_file << setw(15) <<        "buffer =" << setw(10) << buff_length                 << "\n";
    
    log_file<<"scrm model parameters: \n";
    log_file << setw(17) <<"Extract window =" << setw(10) << this->model->exact_window_length()<< "\n";
    log_file << setw(17) <<   "Random seed =" << setw(10) << this->SCRMparam->random_seed()    << "\n";

    log_file << setw(17) <<   "Sample size =" << setw(10) << this->model->sample_size()        << "\n";
    log_file << setw(17) <<    "Seq length =" << setw(10) << this->model->loci_length()        << "\n";
    log_file << setw(17) << "mutation rate =" << setw(10) << this->model->mutation_rate()      << "\n";
    log_file << setw(17) <<   "recomb rate =" << setw(10) << this->original_recombination_rate_ << "\n";
    //log_file << setw(17) <<   "recomb rate =" << setw(10) << this->model->recombination_rate() << "\n";
    log_file << setw(17) <<"inferred recomb rate = " << setw(10) << this->model->recombination_rate()   << "\n";

    this->model->resetTime();
    log_file<<setw(17)<<"Pop size (at Generation):\n";
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
    
    //this->model->resetTime();
    //log_file << "Migration rate :" << "\n";
    //for (size_t pop_i = 0 ; pop_i < this->model->population_number() ; pop_i++){
        //log_file << setw(3) << " ";
        //for (size_t pop_j = 0 ; pop_j < this->model->population_number() ; pop_j++){
            //log_file << setw(14) << this->model->migration_rate(pop_i, pop_j)  ;
            //}
        //log_file << "\n";
        //}
        
    //log_file << "Inferred Migration rate :" << "\n";
    //for (size_t pop_i = 0 ; pop_i < this->model->population_number() ; pop_i++){
        //log_file << setw(3) << " ";
        //for (size_t pop_j = 0 ; pop_j < this->model->population_number() ; pop_j++){
            //log_file << setw(14) << migrate[pop_i][pop_j]  ;
            //}
        //log_file << "\n";
        //}    
        
    log_file.close();
    }


void PfParam::append_to_count_file( size_t epoch, string label, int from_pop, int to_pop, double opportunity, double count, double weight) {
    string file_name = Count_NAME;
    ofstream count_file( file_name.c_str(), ios::out | ios::app | ios::binary );
    int field_length_1 = 7;
    int field_length_2 = 15;
    count_file << setw(field_length_1) << epoch << " "
               << setw(field_length_1) << label << " "
               << setw(field_length_1) << from_pop << " "
               << setw(field_length_1) << to_pop << " "
               << setw(field_length_2) << opportunity << " "
               << setw(field_length_2) << count << " "
               << setw(field_length_2) << count/(opportunity+1e-10) << " "
               << setw(field_length_2) << 1.0 / (weight/opportunity+1e-10)
               << endl;
    count_file.close();
}


void PfParam::appending_Ne_file( bool hist ){
    string file_name = hist ? HIST_NAME : Ne_NAME ;
    ofstream Ne_file( file_name.c_str(), ios::out | ios::app | ios::binary);   
    if (hist){
        Ne_file << "=========\n"; 
        }
    Ne_file << "RE\t" << this->model->recombination_rate() << "\n";

    this->model->resetTime();
    for (size_t i = 0; i < this->model->change_times_.size()-1; i++){
    Ne_file << "ME\t" << this->model->getCurrentTime() / this->model->default_pop_size / 4  ; 
        for (size_t pop_i = 0 ; pop_i < this->model->population_number() ; pop_i++){
            for (size_t pop_j = 0 ; pop_j < this->model->population_number() ; pop_j++){
                Ne_file <<  "\t" << this->model->migration_rate(pop_i, pop_j)  ;
                }
            if ( pop_i < (this->model->population_number()-1)){ 
                Ne_file << "\t|";
                }
            }
            Ne_file << "\n";

            this->model->increaseTime();
        }
    Ne_file << "ME\t" << this->model->getCurrentTime() / this->model->default_pop_size / 4  ; 
    for (size_t pop_i = 0 ; pop_i < this->model->population_number() ; pop_i++){        
        for (size_t pop_j = 0 ; pop_j < this->model->population_number() ; pop_j++){
            Ne_file <<  "\t" << this->model->migration_rate(pop_i, pop_j)  ;
            }
        if ( pop_i < (this->model->population_number()-1)){ 
            Ne_file << "\t|";
            }
        }
        Ne_file << "\n";
    
    this->model->resetTime();
    for (size_t i = 0; i < this->model->change_times_.size()-1; i++){
        Ne_file << "NE\t" << this->model->getCurrentTime() / this->model->default_pop_size / 4  ; 
        for ( size_t pop_j = 0 ; pop_j < this->model->population_number() ; pop_j++ ){
            Ne_file << "\t" << this->model->population_size(pop_j) / this->model->default_pop_size ;
            }
        Ne_file << "\n" ;
        this->model->increaseTime();
        }
    Ne_file << "NE\t" << this->model->getCurrentTime() / this->model->default_pop_size / 4  ; 
    for ( size_t pop_j = 0 ; pop_j < this->model->population_number() ; pop_j++ ){
        Ne_file << "\t" << this->model->population_size(pop_j) / this->model->default_pop_size ;
        }
    Ne_file << "\n" ;
    Ne_file.close();
    return;
    }        



//void PfParam::get_scrm_argv(){
    //scrm_input += convert_pattern(pattern, top_t);
    //cout << scrm_input << endl;
    //}
