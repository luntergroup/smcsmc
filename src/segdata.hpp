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


#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <cassert>      // assert
#include <stdlib.h>     /* strtol, strtod */
#include <climits>      // INT_MAX
#include <cmath>        // ceil
#include "exception.hpp"

using namespace std;

#ifndef SEGDATA
#define SEGDATA


struct InvalidSeg : public InvalidInput{
    InvalidSeg( string str ):InvalidInput( str ){
    }
    virtual ~InvalidSeg() throw() {}
};


struct InvalidInputFile : public InvalidSeg{
    InvalidInputFile( string str ):InvalidSeg( str ){
        this->reason = "Invalid input file: ";
        throwMsg = this->reason + this->src;
    }
    ~InvalidInputFile() throw() {}
};


struct WrongNumberOfEntry : public InvalidSeg{
    WrongNumberOfEntry( string str ):InvalidSeg( str ){
        this->reason = "Number of variant site is wrong: ";
        throwMsg = this->reason + this->src;
    }
    ~WrongNumberOfEntry() throw() {}
};


struct InvalidSegmentStartPosition : public InvalidSeg{
    InvalidSegmentStartPosition( string str1, string str2 ):InvalidSeg( str1 ){
        this->reason = "Segment start position at:";
        throwMsg = this->reason + this->src + string(" expect ") + str2;
    }
    ~InvalidSegmentStartPosition() throw() {}
};


struct NoDataError : public InvalidSeg{
    NoDataError( string str, int start, int end ):InvalidSeg( str ){
        this->reason = "No data found in file ";
        throwMsg = this->reason + str + " between positions " +to_string(start) + " and " + to_string(end);
    }
    ~NoDataError() throw() {}
};


enum Segment_State {SEGMENT_INVARIANT, SEGMENT_MISSING, SEGMENT_INVARIANT_PARTIAL};

class Segment {
    Segment( string file_name , size_t nsam, double seqlen, double num_of_mut,
             long long data_start = 1, double max_segment_length = 1e99 );
    ~Segment(){};

    friend class PfParam;
    friend class ParticleContainer;
    #ifdef UNITTEST
    friend class TestPfParam;
    friend class TestSegment;
    #endif

    const string file_name_;
    const size_t nsam_;
    const long long data_start_;
    const double seqlen_;
    const double max_segment_length_;
    double segment_start_;
    double segment_length_;
    Segment_State segment_state_;
    size_t chrom_;
    bool empty_file_;

    // buffering
    class SegDatum {
    public:
        long long segment_start;
        long long segment_length;
        Segment_State segment_state;
        vector<int> allele_state;
        SegDatum( double ss, double sl, Segment_State s, const vector<int>& as ) :
            segment_start(ss), segment_length(sl), segment_state(s), allele_state(as) {}
    };
    vector< SegDatum > buffer;
    size_t current_buf_index_;
    vector<bool> phased;

    // Methods
    void prepare();
    void set_lookahead();
    vector<int> extract_field_VARIANT( const string field );
    void calculate_num_of_expected_mutations ( size_t nsam, double theta );
    Segment_State segment_state() const { return this->segment_state_; }

    double num_of_expected_mutations_;
    bool end_data_;
public:
    // The allele information
    vector <int> allelic_state_at_Segment_end;

    // Lookahead for Auxiliary Particle Filter implementation
    class Doubleton {
    public:
        int seq_idx_1, seq_idx_2;
        double first_evidence_distance, last_evidence_distance;
	bool unphased_1, unphased_2;
        bool incompatible;
        Doubleton( int i1, int i2, double d1, double d2, bool uph_1, bool uph_2 ) :
	  seq_idx_1(i1), seq_idx_2(i2), first_evidence_distance(d1), last_evidence_distance(d2), unphased_1(uph_1), unphased_2(uph_2), incompatible(false) {}
    };

    // singleton data:
    vector<double> first_singleton_distance;
    vector<double> relative_mutation_rate;       // to account for partially missing data
    // doubleton data:
    vector<Doubleton> doubleton;
    // split data:
    double first_split_distance;
    vector<int> allelic_state_at_first_split;
    int mutation_count_at_first_split;

    // Public methods
    bool empty_file () const      { return this->empty_file_; }
    double segment_start() const  { return segment_start_; }
    double segment_length() const { return segment_length_; }
    double segment_end() const    { return segment_start_ + segment_length_; }

    void read_new_line();
    void reset_data_to_first_entry();

    bool end_data() const {return end_data_; }
    void set_end_data ( bool condition ) { this->end_data_ = condition ; }
};

#endif
