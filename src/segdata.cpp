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

#include "segdata.hpp"
using namespace std;

Segment::Segment( string file_name, size_t nsam, double seqlen, size_t num_of_mut ){
    
    this->file_name_ = file_name;
    this->nsam_ = nsam;
    this->seqlen_ = seqlen;
    this->segment_start_ = 1;
    this->current_line_index_ = 0;
    this->end_data_ = false;
    this->empty_file = false;
    this->num_of_expected_mutations_ = 0;
    
    if ( this->file_name_.size() == 0 ){
        this->empty_file = true;
        this->calculate_num_of_expected_mutations( nsam, num_of_mut );
        this->reset_empty_entry();
        for ( size_t i = 0; i < nsam_; i++){
            this->allelic_state_at_Segment_start.push_back ( -1 );
        }
        return;
    }
     
    this->init();
}


void Segment::init(){

    ifstream in_file;
    in_file.open( this->file_name_.c_str() );
    if ( in_file.good() ){
        getline ( in_file, this->tmp_line ); 
        //dout << this->tmp_line <<endl;
        while ( this->tmp_line.size() > 0 ){
            if ( this->tmp_line[0] != '#' ){
                buffer_lines.push_back ( this->tmp_line );
            }
            getline ( in_file, this->tmp_line );
        }    
    
    }
    else {
        throw std::invalid_argument ( std::string("Invaild file: ") + this->file_name_ );
        }
    in_file.close();
    
    this->segment_length_ = 0;    
    
    this->read_new_line ();
    
    
    this->segment_start_ = 1;
    this->segment_length_ = 0;
    this->current_line_index_ = 0;
}


void Segment::initialize_read_newLine(){
    this->feild_start = 0;
    this->field_end   = 0;
    this->field_index = 0;
    
    this->segment_start_ += this->segment_length_;    
    // check for genetic break, see if starts from a new chrom 
    }

void Segment::read_new_line(){
    /*! Read Segment data, extract mutation site and haplotype
     */ 

    this->initialize_read_newLine(); 
    
    if ( this->empty_file ){ 
        //this->end_data_ = true; 
        return; }
    
    if ( this->current_line_index_ == this->buffer_lines.size() ) {
        this->end_data_ = true;
        return;
    }

    this->tmp_line = this->buffer_lines[this->current_line_index_];

    while ( field_end < this->tmp_line.size() ){
        field_end = min ( this->tmp_line.find('\t',feild_start), this->tmp_line.find('\n', feild_start) );
        this->tmp_str = this->tmp_line.substr( feild_start, field_end - feild_start );
        if      ( field_index == 0 ) { assert ( strtol( tmp_str.c_str(), NULL, 0) == this->segment_start_ );  }  
        else if ( field_index == 1 ) { this->segment_length_ =  strtol( tmp_str.c_str(), NULL, 0);    }  
        else if ( field_index == 2 ) { this->segment_state_ = ( this->tmp_str == "T" ) ? SEGMENT_INVARIANT : SEGMENT_MISSING; }
        else if ( field_index == 3 ) { this->genetic_break_ = ( this->tmp_str == "T" ) ? true : false; }  
        else if ( field_index == 4 ) { this->chrom_ = strtol( tmp_str.c_str(), NULL, 0); }  
        else if ( field_index == 5 ) { this->extract_field_VARIANT(); }
        
        feild_start = field_end+1;        
        field_index++;
        }
        
    if ( current_line_index_%300 == 0 ){ cout << " Reading the " << current_line_index_ << "th entry"<<endl; }
    this->current_line_index_++;
}
    
    
void Segment::extract_field_VARIANT ( ){
    allelic_state_at_Segment_start.clear();

    assert( nsam_ == this->tmp_str.size() );
    for ( size_t i = 0; i < nsam_; i++){
        int seg_contant = ( this->tmp_str[i] == '.' ) ? -1 : strtol(this->tmp_str.substr(i, 1).c_str(), NULL, 0);
        this->allelic_state_at_Segment_start.push_back ( seg_contant );
        //this->allelic_state_at_Segment_start.push_back ( strtol(this->tmp_str.substr(i, 1).c_str(), NULL, 0) );
        //cout << this->allelic_state_at_Segment_start.back();
    }
    //cout<<endl;
}    


void Segment::calculate_num_of_expected_mutations ( size_t nsam, size_t theta ){
    double sum_of_one_over = 0;
    for ( double i = 1; i < (double)nsam ; i++ ){
        sum_of_one_over += 1 / i;
    }
    sum_of_one_over *= (double)theta;
    this->num_of_expected_mutations_ = (size_t) sum_of_one_over;
}
