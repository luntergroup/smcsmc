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

#include <iostream>
#include <cstring> // strlen
#include <vector>
#include <math.h>
#include <string>
#include <stdexcept>


//#ifndef NDEBUG
//#define dout std::cout
//#else
//#pragma GCC diagnostic ignored "-Wunused-value"
//#define dout 0 && std::cout
//#endif

using namespace std;


class Pattern{
    friend class PfParam;
    friend class TestPattern;
    
    Pattern (string pattern, double top_t);
    ~Pattern( ){};
    
    // Methods
    size_t extract_Number( );
    size_t extract_SegmentFactors( ) ;
    void check_pattern ( );
    void extract_NumberOfSegment ( ) ;
    vector <double> extract_Segment ( );
    vector <double> regroup_Segment (vector <double> old_seg);

    void init(){
        this->expr_ = NULL;
        this->num_seg_ = 0;
        }
    
    // Members
    string pattern_str;
    vector <size_t> seg_level1_vec_;
    vector <size_t> seg_level2_vec_;
    const char * expr_;
    double top_t_;
    size_t num_seg_;
    };


