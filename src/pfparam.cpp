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


#include "pfparam.hpp"
#include "pattern.hpp"
#include <ostream>
#include <iostream>
#include <chrono>
#include <ctime>

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
            this->ESS_fraction = readNextInput<double>();
            this->ESS_default_bool = false;
            if ( this->ESS_fraction > 1.0 || this->ESS_fraction < 0.0 ){
                throw OutOfRange ("-ESS", *argv_i );
            }
        } else if ( *argv_i == "-arg" ) {
            record_trees = true;
        } else if ( *argv_i == "-EM"   ){
            this->EM_steps = readNextInput<int>();
            this->EM_bool = true;
        } else if ( *argv_i == "-xr" || *argv_i == "-xc" ) {
            string tmpFlag = *argv_i;
            int last_epoch;
            int first_epoch = readRange(last_epoch);        // obtain 0-based closed interval
            last_epoch++;                                   // turn into 0-based half-open
            for (int i=0; i<last_epoch; i++) {
                if (record_event_in_epoch.size() <= i) {
                    // extend vector, and set default: record both recomb and coal/migr events, don't record ARG
                    record_event_in_epoch.push_back( PfParam::RECORD_COALMIGR_EVENT
                                                     | PfParam::RECORD_RECOMB_EVENT );
                }
                if (i >= first_epoch) {
                    // reset bit signifying recording of either recombination or coal/migr events
                    // " &= " is and-update (cf. +=, sum-update); " ~ " is bitwise not
                    if (tmpFlag == "-xc")      record_event_in_epoch[i] &= ~PfParam::RECORD_COALMIGR_EVENT;
                    else if (tmpFlag == "-xr") record_event_in_epoch[i] &= ~PfParam::RECORD_RECOMB_EVENT;
                }
            }
        } else if ( *argv_i == "-cap"  ){
            this->Ne_cap = readNextInput<double>();
            this->useCap = true;
        } else if ( *argv_i == "-tmax" ){
            this->top_t_ = readNextInput<double>();
        } else if ( *argv_i == "-p"    ){
            this->nextArg();
            this->pattern = *argv_i;
        // ------------------------------------------------------------------
        // Input files
        // ------------------------------------------------------------------
        } else if ( *argv_i == "-seg"  ){
            this->nextArg();
            this->input_SegmentDataFileName = *argv_i;
        } else if ( *argv_i == "-guide" ) {
            this->nextArg();
            this->input_RecombinationBiasFileName = *argv_i;
        } else if ( *argv_i == "-startpos" ) {
            this->start_position = readNextInput<double>();
            if (start_position < 1) {
                throw OutOfRange("-startpos", *argv_i);
            }
        // ------------------------------------------------------------------
        // Action
        // ------------------------------------------------------------------
        } else if ( *argv_i == "-lag"    ){
            this->lag = this->readNextInput<double>();
            this->calibrate_lag = false;
        } else if ( *argv_i == "-calibrate_lag"){
            this->calibrate_lag = true;
            this->lag_fraction = readNextInput<double>();
            if ( this->lag_fraction < 0.0 ){
                throw OutOfRange ("-calibrate_lag", *argv_i );
            }
        } else if ( *argv_i == "-delay" ) {
            this->delay = this->readNextInput<double>();
            if (this->delay < 0.0) {
                throw OutOfRange( "-delay", *argv_i );
            }
        } else if ( *argv_i == "-delay_coal" ) {
            this->delay_type = RESAMPLE_DELAY_COAL;
        } else if ( *argv_i == "-delay_migr" ) {
            this->delay_type = RESAMPLE_DELAY_COALMIGR;
        } else if ( *argv_i == "-ancestral_aware" ){
            this->ancestral_aware = true;
        } else if ( *argv_i == "-dephase" ){
            this->dephase = true;
	} else if ( *argv_i == "-apf" ) {
	    this->auxiliary_particle_filter = this->readNextInput<int>();
	    if (auxiliary_particle_filter < 0 || auxiliary_particle_filter > 4) {
		throw OutOfRange( "-apf", *argv_i );
	    }
        // ------------------------------------------------------------------
        // Output
        // ------------------------------------------------------------------
        } else if ( *argv_i == "-o"     ){
            this->nextArg();
            this->out_NAME_prefix = *argv_i;
        } else if ( *argv_i == "-log"   ){
            this->log_bool  = true;
        } else if (*argv_i == "-record_ess" ) {
            this->record_resample_file = true;
        } else if (*argv_i == "-h" || *argv_i == "-help") {
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

    clog << "Command line:-" << endl << argv[0];
    for (int i=1; i<argc; i++) clog << " " << argv[i];
    clog << endl;
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
    this->max_segment_length_factor_ = 2.0;  // allow segments of max length   mslf / (4 Ne rho)
    this->N                = 100;
    this->lag              = 0.0;
    this->calibrate_lag    = true;
    this->lag_fraction     = 2.0;
    this->delay            = 0.5;
    this->delay_type       = RESAMPLE_DELAY_RECOMB;
    this->ancestral_aware  = false;
    this->dephase          = false;
    this->auxiliary_particle_filter = 0;
    this->out_NAME_prefix  = "smcsmc";
    this->ESS_fraction     = 0.5;
    this->ESS_default_bool = true;
    this->log_bool         = true; // Enable log by default
    this->online_bool      = false;
    this->EM_steps         = 0;
    this->EM_bool          = false;

    this->Segfile          = NULL;
    this->start_position   = 1;
    this->SCRMparam        = NULL;
    this->rg               = NULL;
    this->scrm_input       = "";
    this->top_t_           = 2;
    this->EMcounter_       = 0;
    this->argc_            = 0;
    this->useCap           = false;
    this->Ne_cap           = 200000;
    this->record_resample_file = false;
    this->record_trees     = false;
    this->setHelp(false);
    this->setVersion(false);
    this->input_SegmentDataFileName = "";
    this->input_RecombinationBiasFileName = "";
    this->pattern = "";
    this->record_event_in_epoch.clear();
}


void PfParam::insertMutationRateInScrmInput(){
    size_t found = scrm_input.find("-t");
    if ( found == std::string::npos ) {
      throw std::invalid_argument("The option -t is required");
    }
}


void PfParam::insertRecombRateAndSeqlenInScrmInput(){
    size_t found = scrm_input.find("-r");
    if ( found == std::string::npos ) {
      throw std::invalid_argument("The option -r is required");
    }
    // If recombination and sequence length are defined in the input, extract the sequence length
    // skipping the number of recombination first
    size_t pos_start = scrm_input.find(" ", found+2, 1);
    size_t pos_end = scrm_input.find(" ", pos_start+1, 1);
    // extract the loci length
    pos_start = scrm_input.find(" ", pos_start+2, 1);
    pos_end = scrm_input.find(" ", pos_start+1, 1);
    this->default_loci_length = std::strtod ( (char*)scrm_input.substr(pos_start, pos_end - pos_start).c_str(), NULL);
}


void PfParam::insertSampleSizeInScrmInput(){
    this->scrm_input = to_string ( this->default_nsam ) + " 1 " + this->scrm_input;
}


void PfParam::finalize_scrm_input (  ) {
    // These options insert parameters to the beginning of the current scrm_input
    this->insertRecombRateAndSeqlenInScrmInput(); // This must be run before insertMutationRateInScrmInput, need to define the sequence length
    this->insertMutationRateInScrmInput();
    this->insertSampleSizeInScrmInput();
    if (pattern.size() > 0) {
        Pattern tmp_pattern( pattern, top_t_ );
        this->scrm_input = "scrm " + this->scrm_input + tmp_pattern.pattern_str;
    } else {
        this->scrm_input = "scrm " + this->scrm_input;
    }
    this->convert_scrm_input ();
}


void PfParam::convert_scrm_input (){
    enum { kMaxArgs = 1024 };
    int scrm_argc = 0;
    char *scrm_argv[kMaxArgs];
    cout << "Scrm input: " << scrm_input << endl;//DEBUG
    char * p2 = strtok((char *)this->scrm_input.c_str(), " ");
    while (p2 && scrm_argc < kMaxArgs) {
        scrm_argv[scrm_argc++] = p2;
        p2 = strtok(0, " ");
    }
    ///*! Extract scrm parameters */
    this->SCRMparam = new Param(scrm_argc, scrm_argv, false);
    this->model = this->SCRMparam->parse();
    this->model.has_window_seq_ = true;
    this->rg = new MersenneTwister(this->SCRMparam->seed_is_set(), this->SCRMparam->random_seed());  /*! Initialize mersenneTwister seed */
    this->original_recombination_rate_ = model.recombination_rate();
}


void PfParam::finalize(){

    this->ESSthreshold           = this->N * this->ESS_fraction;
    this->outFileName            = out_NAME_prefix + ".out";
    this->log_NAME               = out_NAME_prefix + ".log";
    this->recombination_map_NAME = out_NAME_prefix + ".recomb.gz";
    this->resample_NAME          = out_NAME_prefix + ".resample";
    this->tree_NAME              = out_NAME_prefix + ".trees.gz";

    // currently, using guided sampling is incompatible with auxiliary particle filters (see includeLookaheadLikelihood)
    if ( (input_RecombinationBiasFileName.size() > 0) && (auxiliary_particle_filter > 0) ) {
	throw std::invalid_argument(std::string("Recombination guiding and auxiliary particle filters cannot currently be used together"));
    }

    // remove any existing files with these names
    remove( this->outFileName.c_str() );
    remove( this->log_NAME.c_str() );
    remove( this->recombination_map_NAME.c_str() );
    if (record_resample_file)
        remove( this->resample_NAME.c_str() );

    this->finalize_scrm_input();
    if ( model.change_times().back() >= top_t_ * 40000 ) {
        throw std::invalid_argument(std::string("Problem: -tmax must be larger than bottom of final epoch"));
    }

    // if necessary, extend the vector specifying what epochs to collect events for,
    // and check it hasn't been made too large (which wouldn't strictly be a problem,
    // but clearly the user specified something that's silly and should hear that.)
    while (record_event_in_epoch.size() < this->model.change_times_.size()) {
        record_event_in_epoch.push_back( PfParam::RECORD_COALMIGR_EVENT | PfParam::RECORD_RECOMB_EVENT );
    }
    if (record_trees) {
        for (unsigned int i = 0; i < record_event_in_epoch.size(); i++) {
            record_event_in_epoch[i] |= PfParam::RECORD_TREE_EVENT;
        }
    }
    if (record_event_in_epoch.size() > this->model.change_times_.size()) {
        //throw std::invalid_argument(std::string("Problem: epochs specified in -xr/-xc options out of range"));
        throw OutOfEpochRange(to_string(record_event_in_epoch.size()-1), to_string(this->model.change_times_.size()-1));
    }

     /*! Initialize seg file, and data up to the first data entry says "PASS"   */
    int max_seg_len = (int)(max_segment_length_factor_ / (model.recombination_rate() * 4 * model.default_pop_size()));
    this->Segfile = new Segment( this->input_SegmentDataFileName,
                                 this->default_nsam,
                                 (double)this->model.loci_length(),
                                 this->default_num_mut,
                                 this->start_position,
                                 max_seg_len );

    /* set true rate, and number of samples, in RecombinationBias */
    /* but do not actually set the recombination rates yet, to allow unbiased simulation */
    /* To set the rates, call setModelRates() */
    recomb_bias.set_rate( model.recombination_rate() );
    recomb_bias.set_num_leaves( default_nsam );
    if (input_RecombinationBiasFileName.size() > 0) {
        recomb_bias.parse_recomb_bias_file( input_RecombinationBiasFileName );
    }
}


void PfParam::setModelRates() {
    cout << "Setting model rates" << endl;
    if (input_RecombinationBiasFileName.size() > 0) {
        recomb_bias.set_model_rates( model );
    }        
}



int PfParam::log(){
    if (log_bool){
        ofstream logFile(log_NAME.c_str(), ios::out | ios::app | ios::binary);
        this->writeLog(&logFile );
        logFile.close();
    }
    this->writeLog(&std::cout);
    return 0;
}


void PfParam::writeLog(ostream * writeTo){
    (*writeTo) << "###########################\n";
    (*writeTo) << "#        smcsmc log       #\n";
    (*writeTo) << "###########################\n";
    this->printVersion(writeTo);
    (*writeTo) << "smcsmc parameters: \n";

    (*writeTo) << "Segment Data file: "
               << (( this->input_SegmentDataFileName.size() == 0 )? "empty" : input_SegmentDataFileName.c_str() ) << "\n";
    (*writeTo) << "Recombination bias file: "
               << (input_RecombinationBiasFileName.size() == 0 ? "None" : input_RecombinationBiasFileName.c_str() ) << "\n";
    (*writeTo) << setw(15) <<     " EM steps =" << setw(10) << EM_steps                    << "\n";
    if (lag > 0){
        (*writeTo) << setw(15) << "Constant lag =" << setw(10) << lag                      << "\n";
        }
    if (online_bool){
        (*writeTo) << setw(15) << "Online update = TRUE\n";
    }
    (*writeTo) << setw(15) <<             "N =" << setw(10) << N                           << "\n";
    (*writeTo) << setw(15) <<           "ESS =" << setw(10) << ESS_fraction;
    if (ESS_default_bool){ (*writeTo) << " (by default)";}                        (*writeTo) << "\n";
    //(*writeTo) << setw(15) <<        "buffer =" << setw(10) << buff_length                 << "\n";

    (*writeTo)<<"scrm model parameters: \n";
    (*writeTo) << setw(17) <<"Extract window =" << setw(10) << this->model.window_length_seq()<< "\n";
    (*writeTo) << setw(17) <<"Extract window R =" << setw(10) << this->model.window_length_rec()<< "\n";

    //(*writeTo) << setw(17) <<   "Random seed =" << setw(10) << this->SCRMparam->random_seed()    << "\n";

    (*writeTo) << setw(17) <<   "Sample size =" << setw(10) << this->model.sample_size()        << "\n";
    (*writeTo) << setw(17) <<    "Seq length =" << setw(10) << this->model.loci_length()        << "\n";
    (*writeTo) << setw(17) << "mutation rate =" << setw(10) << this->model.mutation_rate()      << "\n";
    (*writeTo) << setw(17) <<   "recomb rate =" << setw(10) << this->original_recombination_rate_ << "\n";
    (*writeTo) << setw(17) <<"inferred recomb rate = " << setw(10) << this->model.recombination_rate()   << "\n";

    this->model.resetTime();
    (*writeTo)<<setw(17)<<"Pop size (at Generation):\n";
    for (size_t i = 0; i < this->model.change_times_.size()-1; i++){
        (*writeTo)<<setw(3)<<"(" << setw(8) << this->model.getCurrentTime() <<" )";
        for (size_t pop_j = 0 ; pop_j < this->model.population_number() ; pop_j++){
            (*writeTo) << " | " << setw(10)<<this->model.population_size(pop_j);
            }

        (*writeTo)<< "\n";
        this->model.increaseTime();
        }
    (*writeTo)<<setw(3)<<"(" << setw(8) << this->model.getCurrentTime() <<" )" ;
    for (size_t pop_j = 0 ; pop_j < this->model.population_number() ; pop_j++){
        (*writeTo) << " | " << setw(10)<<this->model.population_size(pop_j);
    }
    (*writeTo)<< "\n";

    (*writeTo) << "Out file is saved in file: "  << outFileName   << "\n";
}


void PfParam::outFileHeader(){
    string file_name = outFileName;
    ofstream count_file( file_name.c_str(), ios::binary );
    int field_length_1 = 6;
    int field_length_2 = 14;
    count_file << setw(field_length_1) << "Iter"   << " "
               << setw(field_length_1) << "Epoch"  << " "
               << setw(field_length_2) << "Start"  << " "
               << setw(field_length_2) << "End"    << " "
               << setw(field_length_1) << "Type"   << " "
               << setw(field_length_1) << "From"   << " "
               << setw(field_length_1) << "To"     << " "
               << setw(field_length_2) << "Opp"    << " "
               << setw(field_length_2) << "Count"  << " "
               << setw(field_length_2) << "Rate"   << " "
               << setw(field_length_2) << "Ne"     << " "
               << setw(field_length_2) << "ESS"
               << endl;
    count_file.close();

}


class FormatDouble {
public:
    FormatDouble( double d, double scientific_bound = 0.1, int precision = 2 ) : d(d), scientific_bound(scientific_bound), precision(precision) {}
    double d, scientific_bound;
    int precision;
};

ostream& operator<<(ostream& ostr, const FormatDouble& fd) {
    const int field_length = 14;
    const double maxdouble = exp( (field_length - fd.precision - 1) * log(10.0) );
    if (fd.d < maxdouble && (fd.d > fd.scientific_bound || fd.d == 0.0)) {
        return ostr << setw(field_length) << fixed << setprecision(fd.precision) << fd.d;
    } else {
        return ostr << setw(field_length) << scientific << setprecision(field_length-7) << fd.d;
    }
}


void PfParam::appendToOutFile( size_t EMstep,
                               int epoch,
                               double epochBegin,
                               double epochEnd,
                               string eventType,
                               int from_pop,
                               int to_pop,
                               double opportunity,
                               double count,
                               double weight) {
    string file_name = outFileName;
    ofstream count_file( file_name.c_str(), ios::out | ios::app | ios::binary );
    int field_length_1 = 6;
    count_file << setw(field_length_1) << EMstep << " "
               << setw(field_length_1) << epoch << " "
               << FormatDouble(epochBegin) << " "
               << FormatDouble(epochEnd) << " "
               << setw(field_length_1) << eventType << " "
               << setw(field_length_1) << from_pop << " "
               << setw(field_length_1) << to_pop << " "
               << FormatDouble(opportunity) << " "
               << FormatDouble(count) << " "
               << FormatDouble(count/(opportunity+1e-10)) << " "
               << FormatDouble( (eventType=="Coal") ? (opportunity+1e-10)/(2.0*count) : 0.0 ) << " "
               << FormatDouble( 1.0 / (weight/opportunity+1e-10), 1.0, 3 )
               << endl;
    count_file.close();
}


void PfParam::append_resample_file( int position, double ESS) const {

    if (!record_resample_file) return;
    string file_name = resample_NAME;
    ofstream resample_file( file_name.c_str(), ios::out | ios::app | ios::binary );
    resample_file << position << "\t" << ESS << endl;
    resample_file.close();

}


void PfParam::helpOption(){
    cout << "Options:" << endl;
    cout << setw(15)<<"-Np"            << setw(8) << "INT" << "  --  " << "Number of particles [ " << N << " ]" << endl;
    cout << setw(15)<<"-seg"           << setw(8) << "STR" << "  --  " << "Data file in seg format [ Chrom1.seg ]" << endl;
    cout << setw(15)<<"-o"             << setw(8) << "STR" << "  --  " << "Prefix for output files" << endl;
    cout << setw(15)<<"-EM"            << setw(8) << "INT" << "  --  " << "EM iterations [ 20 ]" << endl;
    cout << setw(15)<<"-startpos"      << setw(8) << "INT" << "  --  " << "First nucleotide position to analyze [ 1 ]" << endl;
    cout << setw(15)<<"-apf"           << setw(8) << "INT" << "  --  " << "Use auxiliary particle filter (0=no, 1=singletons, 2=+doubletons, 3=+first split) [ 0 ]" << endl;
    cout << setw(15)<<"-log"           << setw(8) << " "   << "  --  " << "Generate *.log file" << endl;
    cout << setw(15)<<"-v"             << setw(8) << " "   << "  --  " << "Display timestamp and git versions" << endl;
    cout << endl << "Inference tuning:" << endl;
    cout << setw(15)<<"-dephase"       << setw(8) << " "   << "  --  " << "Dephase heterozygous sites (but use phasing for -apf) [ false ]" << endl;
    cout << setw(15)<<"-calibrate_lag" << setw(8) << "FLT" << "  --  " << "Lag before extracting events (multiple of survival time) [ " << lag_fraction << " ]" << endl;
    cout << setw(15)<<"-tmax"          << setw(8) << "FLT" << "  --  " << "Maximum tree height, in unit of 4N0 [ 3 ]" <<endl;
    cout << setw(15)<<"-ESS"           << setw(8) << "FLT" << "  --  " << "Fractional ESS threshold for resampling [ " << ESS_fraction << " ]" << endl;
    cout << setw(15)<<"-xr"            << setw(8) << "INT" << "  --  " << "Epoch or epoch range to exclude from recombination EM (0-based, closed)" << endl;
    cout << setw(15)<<"-xc"            << setw(8) << "INT" << "  --  " << "Epoch or epoch range (e.g. 0-10) to exclude from coalescent/migration EM" << endl;
    cout << setw(15)<<"-arg"           << setw(8) << " "   << "  --  " << "Record ARGs in .trees file" << endl;
    cout << setw(15)<<"-guide"         << setw(8) << "STR" << "  --  " << "Recombination guide file [ none ]" << endl;
    cout << endl << "Less frequently used options:" << endl;
    cout << setw(15)<<"-record_ess"    << setw(8) << " "   << "  --  " << "Generate *.resample file" << endl;
    cout << setw(15)<<"-delay"         << setw(8) << "FLT" << "  --  " << "How much to delay application of importance weight due to bias (fraction of survival time)[ " << delay << "]" << endl;
    cout << setw(15)<<"-delay_coal"    << setw(8) << " "   << "  --  " << "Application delay depends on time of (first) coalescence event (default: recombination event)" << endl;
    cout << setw(15)<<"-delay_migr"    << setw(8) << " "   << "  --  " << "Application delay depends on time of (first) migration or coalescence event (default: recombination event)" << endl;
    cout << setw(15)<<"-bias_heights"  << setw(8) << "FLT(s)" << "  --  " << "Time boundaries (in generations) between time sections to focus sampling" << endl;
    cout << setw(15)<<"-bias_strengths"<< setw(8) << "FLT(s)" << "  --  " << "Relative sampling focus in time sections; should have one more value than bias_heights" << endl;
    //cout << setw(15)<<"-p"             << setw(8) << "STR" << "  --  " << "Pattern of time segments [ \"3*1+2*3+4\" ]" <<endl;
}


void PfParam::helpExample(){
    cout << "    Examples:" << endl;
    cout << "smcsmc 10 -nsam 3" << endl;
    cout << "./smcsmc -Np 5 -t 0.002 -r 400 -npop 20000 -seg eg_seg.seg -buff 4" << endl;
    cout << "./smcsmc -Np 5 -t 0.002 -r 400 -npop 20000 -seg eg_seg.seg" << endl;
    cout << "./smcsmc -Np 6 -t 0.0002 -r 30 -npop 10000 -seed 1314 -seg eg_seg.seg" << endl;
    cout << "./smcsmc -Np 7 -t 0.002 -log -r 400 -seg eg_seg.seg " << endl;
}


void PfParam::printHelp(){
    cout << "smcsmc --  Sha (Joe) Zhu, Donna Henderson and Gerton Lunter -- Version " << VERSION << endl;
    this->helpOption();
    this->helpExample();
}


void PfParam::printVersion(std::ostream *output){
    (*output) << "Program was compiled on: " << this->compileTime << endl;
    (*output) << "smcsmc version: " << this->smcsmcVersion << endl;
    (*output) << "scrm version:   " << this->scrmVersion   << endl;
}

void PfParam::startTimer(){
    start = std::chrono::system_clock::now();	
}

double PfParam::stopAndPrintTimer(){
    end = std::chrono::system_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end - this->start); 
    return (double) duration.count();
}
