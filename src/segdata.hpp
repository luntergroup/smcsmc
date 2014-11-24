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


#include<string>
#include<fstream>
#include<iostream>
#include<vector>
#include<cassert>       // assert
#include<stdexcept>      // std::invalid_argument
#include<stdlib.h>     /* strtol, strtod */
#include<climits> // INT_MAX

using namespace std;

#ifndef SEGDATA
#define SEGDATA

enum Segment_State {SEGMENT_INVARIANT, SEGMENT_MISSING};

class Segment{
    friend class PfParam;
    friend class ParticleContainer;
    #ifdef UNITTEST
    friend class TestParam;
    #endif


    Segment_State segment_state() const { return this->segment_state_; }
    //bool variant_state() const { return this->variant_state_; }
    bool genetic_break() const { return this->genetic_break_; }
    size_t seqlen_;
    // Important stuff
    size_t segment_start_;
    size_t segment_length_;
    Segment_State segment_state_;
    //bool variant_state_;                   // Variant state (given it is not missing ) at the start of the segment, True or False
    bool genetic_break_;
    size_t chrom_;
    vector <int> allelic_state_at_Segment_start; // Since missing variant can be represented here, variant can be ignored?
    
    bool empty_file;
    size_t nsam_;
    
    // less important stuff

    string file_name_;
    string tmp_line;
    string tmp_str;

    vector <string> buffer_lines;
    size_t current_line_index_; 
    
    // Line related
    size_t feild_start;
    size_t field_end;
    int field_index;
    
    
    Segment( string file_name , size_t nsam, size_t seqlen );
    ~Segment(){};
    
    // Methods
    void init();
    void initialize_read_newLine();
    void extract_field_VARIANT();

    size_t num_of_entries_;
    bool end_data_;
    
public:
    size_t segment_start() const { return this->segment_start_; }
    size_t segment_length() const { return this->segment_length_; }
    size_t segment_end() const { return ( this->segment_start_ + this->segment_length_ ); }

    void read_new_line();
    void reset_data_to_first_entry(){ 
        this->current_line_index_ = 0; 
        this->set_end_data(false);
        this->segment_start_ = 1;
        this->segment_length_ = 0;
        if ( this->empty_file ){
            this->reset_empty_entry();
        }
    };
    
    void reset_empty_entry(){
        this->segment_length_ = this->seqlen_/20 ;
        this->segment_state_ = SEGMENT_MISSING;
        //this->variant_state_ = false;
        this->genetic_break_ = true;         
    }

    bool end_data() const {return end_data_; }
    void set_end_data ( bool condition ) { this->end_data_ = condition ; }
};

#endif
