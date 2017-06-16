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

#include "mersenne_twister.h"
#include "param.h"
#include "segdata.hpp"
#include "exception.hpp"

#include <fstream>
#include <iostream>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/algorithm/string/predicate.hpp>

using namespace std;

#ifndef PFARGPfParam
#define PFARGPfParam


struct NotEnoughArg : public InvalidInput{
    NotEnoughArg( string str ):InvalidInput( str ){
        this->reason = "Not enough parameters when parsing option: ";
        throwMsg = this->reason + this->src;
    }
    ~NotEnoughArg() throw() {}
};


struct UnknowArg : public InvalidInput{
  UnknowArg( string str ):InvalidInput( str ){
    this->reason = "Unknow option: ";
    throwMsg = this->reason + this->src;
  }
  ~UnknowArg() throw() {}
};

struct FlagsConflict : public InvalidInput{
  FlagsConflict( string str1, string str2 ):InvalidInput( str1 ){
    this->reason = "Flag: ";
    throwMsg = this->reason + this->src + string(" conflict with flag ") + str2;
  }
  ~FlagsConflict() throw() {}
};


struct OutOfEpochRange : public InvalidInput{
  OutOfEpochRange( string str1, string str2 ):InvalidInput( str1 ){
    this->reason = "Problem: epochs specified in -xr/-xc options out of range: ";
    this->src      = "\033[1;31m" + str1 + string(" is greater than ") + str2 + "\033[0m";
    throwMsg = this->reason + src;
  }
  ~OutOfEpochRange() throw() {}
};


struct OutOfRange : public InvalidInput{
  OutOfRange( string str1, string str2 ):InvalidInput( str1 ){
    this->reason = "Flag \"";
    throwMsg = this->reason + this->src + string(" ") + str2 + string("\" out of range [0, 1].");
  }
  ~OutOfRange() throw() {}
};


struct WrongType : public InvalidInput{
  WrongType( string str ):InvalidInput( str ){
    this->reason = "Wrong type for parsing: ";
    throwMsg = this->reason + this->src;
  }
  ~WrongType() throw() {}
};


class RecombBiasSegment {
public:
    RecombBiasSegment() {};
    void set_data(int locus, int size, double rate, vector<double>& leaf_rel_rates ) {
        _locus = locus;
        _size = size;
        _rate = rate;
        _leaf_rel_rates = leaf_rel_rates;
    }
    int get_locus() const                  { return _locus; }
    int get_size() const                   { return _size; }
    int get_end() const                    { return _locus + _size; }
    int get_num_leaves() const             { return _leaf_rel_rates.size(); }
    double get_rate() const                { return _rate; }
    double get_leaf_rate( int leaf ) const { return _leaf_rel_rates[leaf]; }
    friend std::istream& operator>>(std::istream& str, RecombBiasSegment& rbs);
    friend std::ostream& operator<<(std::ostream& str, const RecombBiasSegment& rbs) {
        str << rbs._locus << "\t" << rbs._size << "\t" << rbs._rate;
        for (double rrate : rbs._leaf_rel_rates) str << "\t" << rrate;
        return str;
    }
protected:
    int _locus, _size;
    double _rate;
    vector<double> _leaf_rel_rates;
};


inline std::istream& operator>>(std::istream& str, RecombBiasSegment& rbs) {
    std::string line, elt;
    try {
        if (std::getline(str,line)) {
            if (std::count(line.begin(), line.end(), ' ') > 0) {
                cerr << "Found spaces in recombination record; columns must be tab-separated" << endl;
                cerr << "Record: '" << line << "'" << endl;
                throw InvalidInput("Found spaces in recombination record");
            }
            std::stringstream iss(line);
            if ( !std::getline(iss, elt, '\t' ) ) throw InvalidInput("Parse error on 1st element");
            rbs._locus = std::stoi( elt );
            if ( !std::getline(iss, elt, '\t' ) ) throw InvalidInput("Parse error on 2nd element");
            rbs._size = std::stoi( elt );
            if ( !std::getline(iss, elt, '\t' ) ) throw InvalidInput("Parse error on 3rd element");
            rbs._rate = std::stod( elt );
            rbs._leaf_rel_rates.clear();
            while ( std::getline(iss, elt, '\t' ) ) { rbs._leaf_rel_rates.push_back( std::stod( elt ) ); }
        }
    } catch (const InvalidInput& e) {
        throw e;
    } catch (...) {
        throw InvalidInput("Problem reading or parsing recombination guide file");
    }
    return str;
}


class RecombinationBias {
public:
    RecombinationBias() : _true_rate(-1), _num_leaves(-1), _using_posterior_biased_sampling(false) {};
    void set_rate( double rate ) {
        assert( _bias_segments.size() == 0 );
        _true_rate = rate;
    }
    double get_true_rate() const {
        return _true_rate;
    }
    bool using_posterior_biased_sampling() const {
        return _using_posterior_biased_sampling;
    }
    void set_num_leaves( int num_leaves ) {
        assert( _bias_segments.size() == 0 );
        _num_leaves = num_leaves;
    }
    void parse_recomb_bias_file( string filename ) {
        // recombination guide files are always 0-based!
        boost::iostreams::filtering_istream in_file;
        ifstream in_raw;
        if (boost::algorithm::ends_with( filename, ".gz")) {
            in_raw = ifstream(filename.c_str(), std::ifstream::in | std::ios_base::binary );
            in_file.push(boost::iostreams::gzip_decompressor());
        } else {
            in_raw = ifstream(filename.c_str(), std::ifstream::in );
        }
        if (!in_raw.is_open()) {
            cout << "Problem opening file " << filename << endl;
            throw InvalidInput("Recombination guide file could not be opened.");
        }
        in_file.push(in_raw);
        string header;
        std::getline(in_file, header);
        if (header.substr(0,5) != "locus") throw InvalidInput("Expected header line (with columns 'locus', 'size', 'recomb_rate', '1', ...)"
                                                              " in recombination guide file");
        RecombBiasSegment tmp;
        while (in_file >> tmp) {
            if ( tmp.get_num_leaves() != _num_leaves ) {
                cerr << "Problem on record at position " << tmp.get_locus() << " with " << tmp.get_num_leaves() 
                     << " leaf columns; expected " << _num_leaves << endl;
                cerr << tmp << endl;
                throw InvalidInput("Did not find expected number of leaf columns");
            }
            if ( tmp.get_locus() != ((_bias_segments.size() == 0) ? 0 : _bias_segments.back().get_end() )) {
                cerr << "Problem on record at position " << tmp.get_locus() << endl;
                throw InvalidInput("Did not get expected locus position (records should start at 0, and leave no gaps)");
            }
            _bias_segments.push_back( tmp );
        }
    }
    const RecombBiasSegment& get_recomb_bias_segment( int idx ) const {
        return _bias_segments[idx];
    }
    const void set_model_rates( Model& model ) {
        assert(_true_rate > 0);
        for (const RecombBiasSegment& rbs : _bias_segments) {
            double mutation_rate = model.mutation_rate();
            if (rbs.get_locus() < model.loci_length()) {
                // need to set the mutation rate as well, otherwise it get default value -1
                model.setRecombinationRate( rbs.get_rate(), false, false, rbs.get_locus() );
                model.setMutationRate( mutation_rate, false, false, rbs.get_locus() );
            }
        }
        _using_posterior_biased_sampling = true;
    }
private:
    double _true_rate;
    int _num_leaves;
    vector<RecombBiasSegment> _bias_segments;
    bool _using_posterior_biased_sampling;
};


/*!
 * \brief smcsmc parameters
 */
class PfParam{
  #ifdef UNITTEST
  friend class TestPfParam;
  #endif
  friend class ParticleContainer;

  public:
    // Constructors and Destructors
    PfParam();
    ~PfParam();

    //
    // Methods
    //
    void parse(int argc, char *argv[]);
    void printHelp();
    bool help() const { return help_; }
    void printVersion(std::ostream *output);
    bool version() const { return version_; }

    int  log( );
    void appending_Ne_file( bool hist = false );
    void outFileHeader();
    void appendToOutFile( size_t EMstep,
                          int epoch,
                          double epochBegin,
                          double epochEnd,
                          string eventType,
                          int from_pop,
                          int to_pop,
                          double opportunity,
                          double count,
                          double weight);

    void append_resample_file( int position, double ESS ) const;
    const char* get_recombination_map_filename() const { return recombination_map_NAME.c_str(); }

    // ------------------------------------------------------------------
    // PfParameters
    // ------------------------------------------------------------------
    size_t N;            /*!< \brief Number of particles */
    int    EM_steps;     /*!< \brief Number of EM iterations */
    double ESSthreshold; /*!< \brief Effective sample size, scaled by the number of particles = ESS * N , 0 < ESS < 1 */
    double max_segment_length_factor_;  /*!< \brief Allow segments of max length   mslf / (4 Ne rho) */
    bool   useCap;
    double Ne_cap;

    // ------------------------------------------------------------------
    // Action
    // ------------------------------------------------------------------
    static const int RECORD_RECOMB_EVENT = 1;
    static const int RECORD_COALMIGR_EVENT = 2;
    enum ResampleDelayType { RESAMPLE_DELAY_RECOMB, RESAMPLE_DELAY_COAL, RESAMPLE_DELAY_COALMIGR };

    double lag;
    bool calibrate_lag;
    double lag_fraction;
    double delay;
    ResampleDelayType delay_type;
    bool ancestral_aware;
    bool online_bool;
    double start_position;

    Segment * Segfile;
    Param *SCRMparam;
    Model model;
    MersenneTwister *rg ;
    vector<int> record_event_in_epoch;
    RecombinationBias recomb_bias;

    double default_loci_length;

    double get_ESS_fraction() const { return this-> ESS_fraction; }
    
    double top_t() const { return this->top_t_;}

    void increaseEMcounter() { this->EMcounter_++; }
    size_t EMcounter() const { return EMcounter_; }
    void setModelRates();

  private:

    //
    // Methods
    //
    void init();
    void reInit();
    void insertMutationRateInScrmInput();
    void insertRecombRateAndSeqlenInScrmInput();
    void insertSampleSizeInScrmInput();
    void finalize_scrm_input ( );
    void finalize ( );
    void convert_scrm_input();
    void writeLog(ostream * writeTo);


    void nextArg(){
        string tmpFlag = *argv_i;
        ++argv_i;
        if (argv_i == argv_.end() || (*argv_i)[0] == '-' ) {
            throw NotEnoughArg (tmpFlag);
        }
    }


    template<class T>
    T convert(const std::string &arg) {
        T value;
        std::stringstream ss(arg);
        ss >> value;
        if (ss.fail()) {
        throw WrongType(arg);
        }
        return value;
    }


    template<class T>
    T readNextInput() {
        this->nextArg();
        return convert<T>(*argv_i);
    }


    int readRange(int& second) {
        this->nextArg();
        char c;
        double input;
        std::stringstream ss(*argv_i);
        ss >> input;
        second = input;
        if (ss.fail() || input < 1) throw std::invalid_argument(std::string("Failed to parse int (e.g. '1') or int range (e.g. '1-10'): ") + *argv_i);
        if (!ss.get(c)) return input;
        if (c != '-') throw std::invalid_argument(std::string("Parsing failure: expected '-' after int in range expression: ") + *argv_i);
        ss >> second;
        if (ss.fail() || second <= input || ss.get(c)) throw std::invalid_argument(std::string("Failed to parse int or int range: ") + *argv_i);
        return input;
    }

    // Members
    int argc_;
    std::vector<std::string> argv_;
    std::vector<std::string>::iterator argv_i;

    // ------------------------------------------------------------------
    // PfParameters
    // ------------------------------------------------------------------
    bool   ESS_default_bool;
    string scrm_input;
    bool   EM_bool;
    size_t EMcounter_;
    double original_recombination_rate_;

    // ------------------------------------------------------------------
    // Default values
    // ------------------------------------------------------------------
    size_t default_nsam;
    double default_mut_rate;
    double default_recomb_rate;
    double default_num_mut;

    // ------------------------------------------------------------------
    // Input
    // ------------------------------------------------------------------
    string pattern;     /*! population segement pattern */
    double top_t_;
    double ESS_fraction;   // scaled between zero and one

    // Segdata related
    string input_SegmentDataFileName;
    string input_RecombinationBiasFileName;

    // ------------------------------------------------------------------
    // Output Related
    // ------------------------------------------------------------------
    bool log_bool;
    bool record_resample_file;
    string out_NAME_prefix;
    string outFileName;
    string recombination_map_NAME;
    string log_NAME;
    string resample_NAME;

    // ------------------------------------------------------------------
    // Version Info
    // ------------------------------------------------------------------
    string smcsmcVersion;
    string scrmVersion;
    string compileTime;

    // Help related
    bool help_;
    void setHelp(const bool help) { this->help_ = help; }
    bool version_;
    void setVersion(const bool version) { this->version_ = version; }

    void helpOption();
    void helpExample();
};

#endif
