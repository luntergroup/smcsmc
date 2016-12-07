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
#include <cassert>       // assert
#include <stdlib.h>     /* strtol, strtod */
#include <climits> // INT_MAX
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
             double data_start = 1, double max_segment_length = 1e99 );
    ~Segment(){};

    friend class PfParam;
    friend class ParticleContainer;
    #ifdef UNITTEST
    friend class TestPfParam;
    friend class TestSegment;
    #endif

    const string file_name_;
    const size_t nsam_;
    const double data_start_;
    const double seqlen_;
    const double max_segment_length_;
    double segment_start_;
    double segment_length_;
    Segment_State segment_state_;
    bool genetic_break_;
    size_t chrom_;
    bool empty_file_;

    // less important stuff
    string tmp_str;
    vector <string> buffer_lines;
    size_t current_line_index_;

    // Methods
    void prepare();
    void extract_field_VARIANT();
    void calculate_num_of_expected_mutations ( size_t nsam, double theta );
    Segment_State segment_state() const { return this->segment_state_; }
    bool genetic_break() const { return this->genetic_break_; }

    double num_of_expected_mutations_;
    bool end_data_;

public:
    vector <int> allelic_state_at_Segment_end; // Since missing variant can be represented here, variant can be ignored?
    bool empty_file () const      { return this->empty_file_; }
    double segment_start() const  { return segment_start_ - data_start_ + 1; }
    double segment_length() const { return segment_length_; }
    double segment_end() const    { return segment_start_ - data_start_ + 1 + segment_length_; }

    void read_new_line();
    void reset_data_to_first_entry();

    bool end_data() const {return end_data_; }
    void set_end_data ( bool condition ) { this->end_data_ = condition ; }
};

#endif
