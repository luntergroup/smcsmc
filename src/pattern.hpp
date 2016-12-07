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

#include <iostream>
#include <cstring> // strlen
#include <vector>
#include <math.h>
#include <string>
#include "exception.hpp"


using namespace std;

struct InvalidPattern : public InvalidInput{
    InvalidPattern( string str ):InvalidInput( str ){
    }
    virtual ~InvalidPattern() throw() {}
    //virtual const char* what () const noexcept {
        //return throwMsg.c_str();
    //}
};


struct PatternDigitsExpected : public InvalidPattern{
    PatternDigitsExpected( string str1 ):InvalidPattern( str1 ){
        this->reason = "While parsing expression: expected digit as first character in number; got ";
        throwMsg = this->reason + this->src ;
    }
    ~PatternDigitsExpected() throw() {}
};


struct PatternTimesExpected : public InvalidPattern{
    PatternTimesExpected( string str1 ):InvalidPattern( str1 ){
        this->reason = "While parsing expression: expected '*' to separate factors; got ";
        throwMsg = this->reason + this->src ;
    }
    ~PatternTimesExpected() throw() {}
};


struct PatternAddsExpected : public InvalidPattern{
    PatternAddsExpected( string str1 ):InvalidPattern( str1 ){
        this->reason = "While parsing expression: expected '+' to separate factors; got ";
        throwMsg = this->reason + this->src ;
    }
    ~PatternAddsExpected() throw() {}
};


class Pattern{
  friend class PfParam;
  friend class TestPattern;

    Pattern();
    Pattern(string pattern, double top_t);
    ~Pattern(){};

    // Methods
    size_t extract_Number();
    size_t extract_SegmentFactors();
    void extract_NumberOfSegment();
    vector <double> extract_Segment();
    vector <double> regroup_Segment(vector <double> old_seg);
    void init();
    void extractFromPattern(string pattern);


    // Members
    string pattern_str;
    vector <size_t> seg_level1_vec_;
    vector <size_t> seg_level2_vec_;
    const char * expr_;
    double top_t_;
    size_t num_seg_;
};


