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

#include <iomanip>      // std::setw
#include "segdata.hpp"
using namespace std;

Segment::Segment( string file_name, size_t nsam, double seqlen, double num_of_mut,
                  double data_start, double max_segment_length )
    : file_name_(file_name), nsam_(nsam), data_start_(data_start),
      seqlen_(seqlen), max_segment_length_(max_segment_length)
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
    long next_start_pos = -1;
    
    if ( in_file.good() ){
        this->buffer.clear();
        getline ( in_file, tmp_line );
        while ( tmp_line.size() > 0 ){
            if ( tmp_line[0] != '#' ){
                //
                // locate start positions of the columns
                //
                vector<int> col_starts;
                int pos = 0;
                while ( tmp_line[pos] && tmp_line[pos] != '\n' ) {
                    col_starts.push_back( pos );
                    for ( ; tmp_line[pos] && tmp_line[pos] != '\n' && tmp_line[pos] != '\t'; ++pos ) {}
                    if (tmp_line[pos] == '\t') ++pos;
                }
                //
                // parse line
                //
                if (col_starts.size() < 3) throw InvalidSeg("Require 3 or 6 columns");
                char* end_ptr;
                vector<int> allele;
                long new_seg_start = strtol( tmp_line.c_str(), &end_ptr, 10 );
                if (*end_ptr != '\t')
                    throw InvalidSegmentStartPosition(tmp_line, to_string(new_seg_start) );
                long new_seg_len   = strtol( tmp_line.c_str() + col_starts[1], &end_ptr, 10 );
                // don't check the end_ptr -- the last length in our segment files have a decimal point
                if ( ( tmp_line[ col_starts[2] ] == 'T' || tmp_line[ col_starts[2] ]== 'F' ) &&
                     tmp_line[ col_starts[2] + 1] == '\t' ) {
                    // column 3 is T or F; ignore the value, but require the next column to also be T/F,
                    // and the column after that the chromosome
                    if ( col_starts.size() != 6 )
                        throw InvalidSeg("Require 6 columns");
                    if ( ( tmp_line[ col_starts[3] ] != 'T' && tmp_line[ col_starts[3] ] != 'F' ) ||
                         tmp_line[ col_starts[3] + 1 ] != '\t' )
                        throw InvalidSeg("Expected T or F in .seg file column 3 and 4");
                    chrom_ = strtol( tmp_line.c_str() + col_starts[4], &end_ptr, 10 );
                    if (*end_ptr != '\t')
                        throw InvalidSeg("Bad chromosome (not an integer) in column 5");
                    allele = extract_field_VARIANT ( string( tmp_line.c_str() + col_starts[5] ) );
                } else {
                    // column 3 is neither T nor F.  Assume alternative format, with just allele following
                    if ( col_starts.size() != 3 )
                        throw InvalidSeg("Require 3 columns");
                    allele = extract_field_VARIANT ( string( tmp_line.c_str() + col_starts[2] ) );
                    chrom_ = 1; // dummy value
                }
                //
                // validate line
                //
                if ( next_start_pos > -1 && next_start_pos != new_seg_start )
                    throw InvalidSeg("Segments are not consecutive");
                next_start_pos = new_seg_start + new_seg_len;
                if ( new_seg_start >= data_start_ + seqlen_) {
                    // don't require this and following data - bail out
                    break;
                }
                //
                // store data
                //
                Segment_State state;
                if (new_seg_start + new_seg_len > data_start_) {
                    // data overlaps segment of interest -- store
                    // if a very long segment is encountered, feed this to the inference machinery in small
                    // pieces; otherwise events don't get evaluated and pruned causing high memory usage
                    do { 
                        if (new_seg_len > max_segment_length_) {
                            new_seg_len = max_segment_length_;
                            // ensure that likelihood for invariance is charged appropriately, but that no
                            // likelihood for a mutation is charged at the end of the segment, as this should
                            // be accounted for only at the actual end of the segment.  However if the segment
                            // is missing, the initial and final sections should be treated the same.
                            state = SEGMENT_INVARIANT_PARTIAL;
                        } else {
                            state = SEGMENT_INVARIANT;
                        }
                        buffer.push_back( SegDatum( new_seg_start,
                                                    new_seg_len,
                                                    state,
                                                    allele ) );
                        new_seg_start += new_seg_len;
                        new_seg_len = next_start_pos - new_seg_start;
                    } while (new_seg_start < next_start_pos);
                }
            }
            getline ( in_file, tmp_line );
        }
        in_file.close();
    } else {
        throw InvalidInputFile(this->file_name_);
    }
    if (buffer.size() == 0) {
        throw NoDataError(file_name_, data_start_, data_start_ + seqlen_);
    }
}


void Segment::reset_data_to_first_entry() {

    current_buf_index_ = 0;
    set_end_data(false);
    segment_start_ = data_start_;
    segment_length_ = 0;
    if ( empty_file_ ){
        segment_length_ = ceil( (size_t)seqlen_ / num_of_expected_mutations_ ) ;
        segment_state_ = SEGMENT_MISSING;
    }
}


void Segment::read_new_line(){
    /*! Read Segment data, extract mutation site and haplotype
     */

    // move to beginning of next segment
    segment_start_ += segment_length_;

    if ( empty_file() ){
        return;
    }

    if ( current_buf_index_ == buffer.size() ) {
        end_data_ = true;
        return;
    }

    /* new values, using segdatum buffer */
    const SegDatum sd = buffer[ current_buf_index_ ];
    long new_seg_end = sd.segment_start + sd.segment_length;
    if (sd.segment_start > segment_start_)
        throw InvalidSeg("Internal error - segment computation problem (start)");
    if (new_seg_end < segment_start_)
        throw InvalidSeg("Internal error - segment computation problem (end)");

    segment_length_ = new_seg_end - segment_start_;
    segment_state_ = sd.segment_state;
    allelic_state_at_Segment_end = sd.allele_state;
    chrom_ = 1; // dummy value

    if ( current_buf_index_ > 0 && current_buf_index_ % 2000 == 0 ){
    cout << "\r" << " Particle filtering step" << setw(4) << int((segment_end() * 100) / seqlen_)
    << "% completed." << flush;
    }
    
    current_buf_index_++;

}


vector<int> Segment::extract_field_VARIANT ( const string field ) {

    vector<int> allelic_state_at_seg_end;
    assert( nsam_ == field.size() );
    if ( nsam_ != field.size() ){
        throw WrongNumberOfEntry(field);
    }

    for ( size_t i = 0; i < nsam_; i++ ) {
        int seg_content;
        switch (field[i]) {
            case '.': seg_content = -1; break;  // missing data
            case '/': seg_content = 2; break;   // unphased heterozygous genotype
            case '0': seg_content = 0; break;
            case '1': seg_content = 1; break;
            default: throw InvalidSeg("Unknown character found in .seg file; expect one of '.', '/', '0' or '1'.");
        }
        allelic_state_at_seg_end.push_back ( seg_content );
    }
    return allelic_state_at_seg_end;
}


void Segment::calculate_num_of_expected_mutations ( size_t nsam, double theta ){
    double sum_of_one_over = 0.0;
    for ( size_t i = 1; i < nsam ; i++ ){
        sum_of_one_over += 1.0 / i;
    }
    sum_of_one_over *= theta;
    this->num_of_expected_mutations_ = sum_of_one_over;
}
