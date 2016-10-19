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
#include <ostream>
PfParam::PfParam(){
    this->init();  // Initialize pfARG program parameters
}


void PfParam::reInit(){
    if ( this->Segfile != NULL ){
        delete this->Segfile;
    }
    if ( this->SCRMparam != NULL ){
        delete this->SCRMparam;
    }
    if ( this->rg != NULL ){
        delete this->rg;
    }
    this->init();
}


void PfParam::parse(int argc, char *argv[]) {
    this->argc_ = argc;
    this->argv_ = std::vector<std::string>(argv + 1, argv + argc);
    this->argv_i = argv_.begin();

    if ( argv_.size() == 0 ) {
        this->setHelp(true);
        return;
    }

    this->reInit();

    do {
        // ------------------------------------------------------------------
        // Parameters
        // ------------------------------------------------------------------
        if        ( *argv_i == "-Np"   ){
            this->N = readNextInput<size_t>();
        } else if ( *argv_i == "-nsam" ){
            this->default_nsam = readNextInput<size_t>();
        } else if ( *argv_i == "-ESS"  ){
            this->ESS_ = readNextInput<double>();
            this->ESS_default_bool = false;
        } else if ( *argv_i == "-EM"   ){
            this->EM_steps = readNextInput<int>();
            this->EM_bool = true;
        } else if ( *argv_i == "-xr" || *argv_i == "-xc" ) {
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
                    record_event_in_epoch[i] &= ~( (*argv_i == "-xc") ? PfParam::RECORD_COALMIGR_EVENT : PfParam::RECORD_RECOMB_EVENT );
                }
            }
        } else if ( *argv_i == "-cap"  ){
            this->useCap = true;
            this->Ne_cap = readNextInput<double>();
        } else if ( *argv_i == "-tmax" ){
            this->top_t_ = readNextInput<double>();
        } else if ( *argv_i == "-p"    ){
            this->nextArg();
            this->pattern = argv_[argc_i];
        // ------------------------------------------------------------------
        // Input files
        // ------------------------------------------------------------------
        } else if ( *argv_i == "-seg"  ){
            this->nextArg();
            this->input_SegmentDataFileName = argv_[argc_i];
        // ------------------------------------------------------------------
        // Action
        // ------------------------------------------------------------------
        } else if ( *argv_i == "-lag"    ){
            this->lag = this->readNextInput<double>();
        } else if ( *argv_i == "-calibrate_lag"){
            this->calibrate_lag = true;
            this->lag_fraction = readNextInput<double>();
        } else if ( *argv_i == "-online" ){
            this->online_bool = true;
        // ------------------------------------------------------------------
        // Output
        // ------------------------------------------------------------------
        } else if ( *argv_i == "-o"     ){
            this->nextArg();
            this->out_NAME_prefix = argv_[argc_i];
        } else if ( *argv_i == "-log"   ){
            this->log_bool  = true;
        } else if ( *argv_i == "-heat"  ){
            this->heat_bool = true;
        } else if (*argv_i == "-h" || *argv_i == "-help") {
            //Help_header();
            this->setHelp(true);
        } else if (*argv_i == "-v" || *argv_i == "-version") {
            this->setVersion(true);
        } else {
            scrm_input += *argv_i + " ";
        }
    } while ( ++argv_i != argv_.end());

    if ( this->help() || version() ){
        return;
    }

    this->finalize( );
}


PfParam::~PfParam(){
    delete this->Segfile;
    delete this->SCRMparam;
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
    this->lag              = 0;
    this->out_NAME_prefix  = "smcsmc";
    this->ESS_             = 0.5;
    this->ESS_default_bool = true;
    this->log_bool         = true; // Enable log by default
    this->heat_bool        = false;
    this->online_bool      = false;
    this->heat_seq_window  = 400;
    this->EM_steps         = 0;
    this->EM_bool          = false;

    this->Segfile          = NULL;
    this->SCRMparam        = NULL;
    this->rg               = NULL;
    this->scrm_input       = "";
    this->top_t_           = 2;
    this->EMcounter_       = 0;
    this->argc_            = 0;
    this->useCap           = false;
    this->Ne_cap           = 200000;
    this->setHelp(false);
    this->setVersion(false);
}


void PfParam::insert_mutation_rate_in_scrm_input ( ) {

    size_t found = scrm_input.find("-t");
    if ( found == std::string::npos ) {
    // if "-t" option is not found ...
        this->default_num_mut = this->default_mut_rate * this->default_loci_length * 4 * 10000;
        this->scrm_input = "-t " + to_string ( this->default_num_mut ) + " " + this->scrm_input;
    }
}


void PfParam::insert_recomb_rate_and_seqlen_in_scrm_input (  ){

    size_t found = scrm_input.find("-r");
    if ( found == std::string::npos ) { // if "-r" option is not found ...
    // number of recombination events in 4N0, as the scaling N0 in scrm is 10000
        this->scrm_input = "-r " + to_string ( this->default_recomb_rate * this->default_loci_length * 4 * 10000 )
                           + " " + to_string ((size_t)this->default_loci_length) + " " + this->scrm_input;
    } else {
        // skipping the number of recombination first
        size_t pos_start = scrm_input.find(" ", found+2, 1);
        size_t pos_end = scrm_input.find(" ", pos_start+1, 1);
        // extract the loci length
        pos_start = scrm_input.find(" ", pos_start+2, 1);
        pos_end = scrm_input.find(" ", pos_start+1, 1);
        this->default_loci_length = std::strtod ( (char*)scrm_input.substr(pos_start, pos_end - pos_start).c_str(), NULL);
    }
}


void PfParam::insert_sample_size_in_scrm_input (  ) {
    this->scrm_input = to_string ( this->default_nsam ) + " 1 " + this->scrm_input;
}


void PfParam::finalize_scrm_input (  ) {

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
    clog << scrm_input << endl;
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
    //this->SCRMparam->parse( *this->model );
    this->model = this->SCRMparam->parse();
    this->model.has_window_seq_ = true;
    this->rg = new MersenneTwister(this->SCRMparam->seed_is_set(), this->SCRMparam->random_seed());  /*! Initialize mersenneTwister seed */
    this->original_recombination_rate_ = model.recombination_rate();
}


void PfParam::finalize(  ){

    this->ESSthreshold = this->N * this->ESS();
    this->TMRCA_NAME   = out_NAME_prefix + "TMRCA";
    this->WEIGHT_NAME  = out_NAME_prefix + "WEIGHT";
    this->outFileName   = out_NAME_prefix + ".out";
    this->log_NAME     = out_NAME_prefix + ".log";
    //this->SURVIVOR_NAME= out_NAME_prefix + "SURVIVOR";
    this->Resample_NAME= out_NAME_prefix + "Resample";

    // remove any existing files with these names
    remove( this->TMRCA_NAME.c_str() );
    remove( this->WEIGHT_NAME.c_str());
    remove( this->outFileName.c_str() );
    remove( this->log_NAME.c_str()   );
    //remove( this->SURVIVOR_NAME.c_str());
    remove( this->Resample_NAME.c_str());

    this->finalize_scrm_input ( );

    // if necessary, extend the vector specifying what epochs to collect events for,
    // and check it hasn't been made too large (which wouldn't strictly be a problem,
    // but clearly the user specified something that's silly and should hear that.)
    while (record_event_in_epoch.size() < this->model.change_times_.size()) {
        record_event_in_epoch.push_back( PfParam::RECORD_COALMIGR_EVENT | PfParam::RECORD_RECOMB_EVENT );
    }
    if (record_event_in_epoch.size() > this->model.change_times_.size()) {
        throw std::invalid_argument(std::string("Problem: epochs specified in -xr/-xc options out of range"));
    }

     /*! Initialize seg file, and data up to the first data entry says "PASS"   */
    this->Segfile = new Segment( this->input_SegmentDataFileName, this->default_nsam, (double)this->model.loci_length(), this->default_num_mut );
}


int PfParam::log( ) {

    if (log_bool){
        this->log_param( );
        string log_cmd="cat " + log_NAME;
        return system(log_cmd.c_str());
    } else {
        return 0;
    }
}


void PfParam::log_param( ){
    ofstream log_file;
    log_file.open (log_NAME.c_str(), ios::out | ios::app | ios::binary);

    log_file << "###########################\n";
    log_file << "#        smcsmc log       #\n";
    log_file << "###########################\n";
    this->printVersion(&log_file);
    log_file << "smcsmc parameters: \n";
    if (this->heat_bool){
        log_file << "TMRCA saved in file: "  << TMRCA_NAME  << "\n";
        log_file << "WEIGHT saved in file: " << WEIGHT_NAME << "\n";
        //log_file << "BL saved in file: "     << BL_NAME     << "\n";
        }
    //if (this->hist_bool){
        //log_file << "Resample saved in file: " << Resample_NAME << "\n";
        //}

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
    log_file << setw(17) <<"Extract window =" << setw(10) << this->model.window_length_seq()<< "\n";
    //log_file << setw(17) <<   "Random seed =" << setw(10) << this->SCRMparam->random_seed()    << "\n";

    log_file << setw(17) <<   "Sample size =" << setw(10) << this->model.sample_size()        << "\n";
    log_file << setw(17) <<    "Seq length =" << setw(10) << this->model.loci_length()        << "\n";
    log_file << setw(17) << "mutation rate =" << setw(10) << this->model.mutation_rate()      << "\n";
    log_file << setw(17) <<   "recomb rate =" << setw(10) << this->original_recombination_rate_ << "\n";
    log_file << setw(17) <<"inferred recomb rate = " << setw(10) << this->model.recombination_rate()   << "\n";

    this->model.resetTime();
    log_file<<setw(17)<<"Pop size (at Generation):\n";
    for (size_t i = 0; i < this->model.change_times_.size()-1; i++){
        log_file<<setw(3)<<"(" << setw(8) << this->model.getCurrentTime() <<" )";
        for (size_t pop_j = 0 ; pop_j < this->model.population_number() ; pop_j++){
            log_file << " | " << setw(10)<<this->model.population_size(pop_j);
            }

        log_file<< "\n";
        this->model.increaseTime();
        }
    log_file<<setw(3)<<"(" << setw(8) << this->model.getCurrentTime() <<" )" ;
    for (size_t pop_j = 0 ; pop_j < this->model.population_number() ; pop_j++){
        log_file << " | " << setw(10)<<this->model.population_size(pop_j);
    }
    log_file<< "\n";

    log_file << "Out file is saved in file: "  << outFileName   << "\n";

    log_file.close();
}


void PfParam::outFileHeader(){
    string file_name = outFileName;
    ofstream count_file( file_name.c_str(), ios::out | ios::app | ios::binary );
    int field_length_1 = 10;
    int field_length_2 = 15;
    count_file << setw(field_length_1) << "EMstep"      << " "
               << setw(field_length_1) << "EpochIndex"  << " "
               << setw(field_length_1) << "EpochBegin"  << " "
               << setw(field_length_1) << "EpochEnd"    << " "
               << setw(field_length_1) << "EventType"   << " "
               << setw(field_length_1) << "FromPop"     << " "
               << setw(field_length_1) << "ToPop"       << " "
               << setw(field_length_2) << "Opportunity" << " "
               << setw(field_length_2) << "Count"       << " "
               << setw(field_length_2) << "Rate"        << " "
               << setw(field_length_2) << "NE"<< " "
               << setw(field_length_2) << "ESS"
               << endl;
    count_file.close();

}


void PfParam::appendToOutFile( size_t EMstep,
                               int epoch,
                               string epochBegin,
                               string epochEnd,
                               string eventType,
                               int from_pop,
                               int to_pop,
                               double opportunity,
                               double count,
                               double weight) {
    string file_name = outFileName;
    ofstream count_file( file_name.c_str(), ios::out | ios::app | ios::binary );
    int field_length_1 = 10;
    int field_length_2 = 15;
    count_file << setw(field_length_1) << EMstep << " "
               << setw(field_length_1) << epoch << " "
               << setw(field_length_1) << epochBegin << " "
               << setw(field_length_1) << epochEnd << " "
               << setw(field_length_1) << eventType << " "
               << setw(field_length_1) << from_pop << " "
               << setw(field_length_1) << to_pop << " "
               << setw(field_length_2) << opportunity << " "
               << setw(field_length_2) << count << " "
               << setw(field_length_2) << count/(opportunity+1e-10) << " " // Rate
               << setw(field_length_2) << ((eventType=="Coal") ? to_string((opportunity+1e-10)/(2.0*count)) : "-1") << " " // Ne
               << setw(field_length_2) << 1.0 / (weight/opportunity+1e-10) // ESS
               << endl;
    count_file.close();
}


void PfParam::append_resample_file( int position, double ESS) const {
    string file_name = Resample_NAME;
    ofstream resample_file( file_name.c_str(), ios::out | ios::app | ios::binary );
    resample_file << "Resampling at position " << position
                  << " ; ESS is " << ESS
                  << endl;
    resample_file.close();
}


void PfParam::helpOption(){
    //cout << "Too few command line arguments" << endl;
    cout << "Options:" << endl;
    cout << setw(10)<<"-Np"     << setw(5) << "INT" << "  --  " << "Number of particles [ 1000 ]" << endl;
    cout << setw(10)<<"-ESS"    << setw(5) << "FLT" << "  --  " << "Fractional ESS threshold for resampling (1 = use random likelihoods) [ 0.6 ]" << endl;
    cout << setw(10)<<"-p"      << setw(5) << "STR" << "  --  " << "Pattern of time segments [ \"3*1+2*3+4\" ]" <<endl;
    cout << setw(10)<<"-tmax"   << setw(5) << "FLT" << "  --  " << "Maximum time, in unit of 4N0 [ 3 ]" <<endl;
    cout << setw(10)<<"-EM"     << setw(5) << "INT" << "  --  " << "EM iterations [ 20 ]" << endl;
    cout << setw(10)<<"-seg"    << setw(5) << "STR" << "  --  " << "Data file in seg format [ Chrom1.seg ]" << endl;
    cout << setw(10)<<"-o"      << setw(5) << "STR" << "  --  " << "Prefix for output files" << endl;
    cout << setw(10)<<"-online" << setw(5) << " "   << "  --  " << "Perform online EM" << endl;
    cout << setw(10)<<"-xr"     << setw(5) << "INT" << "  --  " << "Epoch or epoch range to exclude from recombination EM (1-based, closed)" << endl;
    cout << setw(10)<<"-xc"     << setw(5) << "INT" << "  --  " << "Epoch or epoch range (e.g. 1-10) to exclude from coalescent/migration EM" << endl;
    cout << setw(10)<<"-log"    << setw(5) << " "   << "  --  " << "Generate *.log file" << endl;
    cout << setw(10)<<"-heat"   << setw(5) << " "   << "  --  " << "Generate *TMRCA and *WEIGHT for heatmap" << endl;
    cout << setw(10)<<"-v"      << setw(5) << " "   << "  --  " << "Display timestamp and git versions" << endl;
};


void PfParam::helpExample(){
    cout << "    Examples:" << endl;
    cout << "smcsmc 10 -nsam 3" << endl;
    cout << "./smcsmc -Np 5 -t 0.002 -r 400 -npop 20000 -seg eg_seg.seg -buff 4" << endl;
    cout << "./smcsmc -Np 5 -t 0.002 -r 400 -npop 20000 -seg eg_seg.seg" << endl;
    cout << "./smcsmc -Np 6 -t 0.0002 -r 30 -npop 10000 -seed 1314 -seg eg_seg.seg" << endl;
    cout << "./smcsmc -Np 7 -t 0.002 -log -r 400 -seg eg_seg.seg " << endl;
};


void PfParam::printHelp(){
    cout << "smcsmc" << endl;
    cout << "  version: " << VERSION << endl;
    cout << "  authored by Sha (Joe) Zhu and Gerton Lunter " <<endl;
    this->helpOption();
    this->helpExample();
};


void PfParam::printVersion(std::ostream *output){
    (*output) << "Program was compiled on: " << this->compileTime << endl;
    (*output) << "smcsmc version: " << this->smcsmcVersion << endl;
    (*output) << "scrm version:   " << this->scrmVersion   << endl;
}
