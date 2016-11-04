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

#include <iomanip>      // std::setw
#include "segdata.hpp"
using namespace std;

Segment::Segment( string file_name, size_t nsam, double seqlen, double num_of_mut, double data_start )
    : file_name_(file_name), nsam_(nsam), data_start_(data_start), seqlen_(seqlen)
{

    num_of_expected_mutations_ = 0;

    if ( this->file_name_.size() == 0 ){

        empty_file_ = true;
        calculate_num_of_expected_mutations( nsam, num_of_mut );
        for ( size_t i = 0; i < nsam_; i++){
            allelic_state_at_Segment_end.push_back ( -1 );
        }

    } else {

        empty_file_ = false;
        prepare();

    }
    reset_data_to_first_entry();
}


void Segment::prepare(){
    ifstream in_file;
    string tmp_line;
    in_file.open( this->file_name_.c_str() );
    if ( in_file.good() ){
        this->buffer_lines.clear();
        getline ( in_file, tmp_line );
        while ( tmp_line.size() > 0 ){
            if ( tmp_line[0] != '#' ){
                int first_tab_pos = tmp_line.find('\t');
                int second_tab_pos = tmp_line.find('\t',first_tab_pos + 1);
                int segstart = strtol( tmp_line.substr( 0, first_tab_pos ).c_str(), NULL, 0 );
                int seglen   = strtol( tmp_line.substr( first_tab_pos+1, second_tab_pos-first_tab_pos ).c_str(),
                                       NULL, 0 );
                if (segstart >= data_start_ + seqlen_) {
                    // simulate end-of-file to bail out of the loop
                    tmp_line = "";
                    break;
                }
                if (segstart + seglen > data_start_) {
                    buffer_lines.push_back ( tmp_line );
                }
            }
            getline ( in_file, tmp_line );
        }
        in_file.close();
    } else {
        throw InvalidInputFile(this->file_name_);
    }

    if (buffer_lines.size() == 0) {
        throw NoDataError(file_name_, data_start_, data_start_ + seqlen_);
    }
}


void Segment::reset_data_to_first_entry() {

    current_line_index_ = 0;
    set_end_data(false);
    genetic_break_ = false;
    segment_start_ = data_start_;
    segment_length_ = 0;
    if ( empty_file_ ){
        segment_length_ = (size_t)seqlen_ / num_of_expected_mutations_ ;
        segment_state_ = SEGMENT_MISSING;
        genetic_break_ = true;
    }
}


void Segment::read_new_line(){
    /*! Read Segment data, extract mutation site and haplotype
     */

    this->field_start = 0;
    this->field_end   = 0;
    this->field_index = 0;

    this->segment_start_ += this->segment_length_;
    // check for genetic break, see if starts from a new chrom
    // todo

    if ( this->empty_file() ){
        return;
    }

    if ( this->current_line_index_ == this->buffer_lines.size() ) {
        this->end_data_ = true;
        return;
    }

    string tmp_line = this->buffer_lines[this->current_line_index_];
    int new_seg_start = -1;

    while ( field_end < tmp_line.size() ){
        field_end = min ( tmp_line.find('\t',field_start), tmp_line.find('\n', field_start) );
        this->tmp_str = tmp_line.substr( field_start, field_end - field_start );
        if ( field_index == 0 ) {
            if (this->genetic_break_){    // Genetic break! reset segMent_start
                this->segment_start_ = data_start_;
            }
            new_seg_start = strtol( tmp_str.c_str(), NULL, 0 );  // check values when we know the length
        } else if ( field_index == 1 ) {
            this->segment_length_ =  strtol( tmp_str.c_str(), NULL, 0 );
            // check that current segment_start_ is inside the new segment
            if (new_seg_start > segment_start_) {
                throw InvalidSegmentStartPosition(tmp_str, to_string(this->segment_start_));
            }
            if (new_seg_start + segment_length_ < segment_start_) {
                throw InvalidSegmentStartPosition(tmp_str, to_string(this->segment_start_));
            }
            // tweak the length so that the endpoint is where we want it to be
            int new_seg_end = new_seg_start + segment_length_;
            segment_start_ = max( (int)segment_start_, new_seg_start );
            segment_length_ = new_seg_end - segment_start_;
        } else if ( field_index == 2 ) {
            this->segment_state_ = ( this->tmp_str == "T" ) ? SEGMENT_INVARIANT : SEGMENT_MISSING;
        } else if ( field_index == 3 ) {
            this->genetic_break_ = ( this->tmp_str == "T" ) ? true : false;
        } else if ( field_index == 4 ) {
            this->chrom_ = strtol( tmp_str.c_str(), NULL, 0);
        } else if ( field_index == 5 ) {
            this->extract_field_VARIANT();
        }
        field_start = field_end+1;
        field_index++;
    }

    if ( current_line_index_ > 0 && current_line_index_ % 300 == 0 ){
        cout << "\r" << " Particle filtering step" << setw(4) << int((segment_end() * 100) / seqlen_)
             << "% completed." << flush;
    }
    this->current_line_index_++;
}


void Segment::extract_field_VARIANT ( ) {

    allelic_state_at_Segment_end.clear();
    assert( nsam_ == tmp_str.size() );
    if ( nsam_ != tmp_str.size() ){
        throw WrongNumberOfEntry(tmp_str);
    }

    for ( size_t i = 0; i < nsam_; i++ ) {
        int seg_content;
        switch (tmp_str[i]) {
            case '.': seg_content = -1; break;  // missing data
            case '/': seg_content = 2; break;   // unphased heterozygous genotype
            default: seg_content = strtol( tmp_str.substr(i,1).c_str(), NULL, 0 );
        }
        allelic_state_at_Segment_end.push_back ( seg_content );
    }
}


void Segment::calculate_num_of_expected_mutations ( size_t nsam, double theta ){
    double sum_of_one_over = 0.0;
    for ( size_t i = 1; i < nsam ; i++ ){
        sum_of_one_over += 1.0 / i;
    }
    sum_of_one_over *= theta;
    this->num_of_expected_mutations_ = sum_of_one_over;
}
