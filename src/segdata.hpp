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


#include<string>
#include<fstream>
#include<iostream>
#include<vector>
#include<cassert>       // assert
#include<stdexcept>      // std::invalid_argument
#include<stdlib.h>     /* strtol, strtod */
#include<climits> // INT_MAX

using namespace std;

#ifndef NDEBUG
#define dout std::cout
#else
//#pragma GCC diagnostic ignored "-Wunused-value"
#pragma clang diagnostic ignored "-Wdeprecated-register"
#define dout 0 && std::cout
#endif

#ifndef SEGDATA
#define SEGDATA

enum Segment_State {SEGMENT_INVARIANT, SEGMENT_MISSING};

class Segment{
    friend class PfParam;
    friend class ParticleContainer;
    #ifdef UNITTEST
    friend class TestParam;
    #endif

    // Important stuff
    int current_site;
    Segment_State current_state;
    bool variant;                   // Variant state (given it is a variant) at the start of the segment, True or False
    vector <int> allelic_state_at_Segment_start; // Since missing variant can be represented here, variant can be ignored?
    bool genetic_break;
    
    // less important stuff
    
    Segment( string file_name );
    ~Segment(){};
    
    // Methods
    void read_new_line();
    void reset_data_to_first_entry();

    ifstream in_file;
    string file_name_;
    string tmp_line;
    string tmp_str;

    vector <string> buffer_lines;
    size_t current_line_index_; 

    
        // Line related
    size_t feild_start;
    size_t field_end;
    int field_index;

};

#endif
