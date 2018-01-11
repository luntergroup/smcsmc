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
                  long long data_start, double max_segment_length )
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
    long long next_start_pos = -1;
    
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
                long long new_seg_start = strtoll( tmp_line.c_str(), &end_ptr, 10 );
                if (*end_ptr != '\t')
                    throw InvalidSegmentStartPosition(tmp_line, to_string(new_seg_start) );
                long new_seg_len   = strtol( tmp_line.c_str() + col_starts[1], &end_ptr, 10 );
                // don't check the end_ptr -- the last length in our segment files have a decimal point
                if ( ( tmp_line[ col_starts[2] ] == 'T' || tmp_line[ col_starts[2] ]== 'F' ) &&
                     tmp_line[ col_starts[2] + 1] == '\t' ) {
                    // column 3 is T or F; ignore the value, but require the next column to also be T/F,
                    // and the column after that the chromosome
                    if ( col_starts.size() != 6 )
                        throw InvalidSeg("Require 6 (or 3) columns");
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
                        throw InvalidSeg("Require 3 (or 6) columns");
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
                        if (new_seg_start + new_seg_len > data_start_) {
                            buffer.push_back( SegDatum( new_seg_start,
                                                        new_seg_len,
                                                        state,
                                                        allele ) );
                        }
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

    // set phasing status for the sequence columns
    this->phased.clear();
    for( SegDatum sd : buffer ) {
        for( size_t idx=0; idx < sd.allele_state.size(); ++idx ) {
            if (phased.size() <= idx) phased.push_back(true);
            if (sd.allele_state[idx] == 2)
                phased[idx] = false;
        }
    }
}


void Segment::reset_data_to_first_entry() {

    current_buf_index_ = 0;
    set_end_data(false);
    segment_start_ = 0;
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
    long new_seg_start = (long)(sd.segment_start - data_start_);
    long new_seg_end = new_seg_start + sd.segment_length;
    if (new_seg_start < 0) new_seg_start = 0;
    if (new_seg_start > segment_start_)
        throw InvalidSeg("Internal error - segment computation problem (start)");
    if (new_seg_end < 0)
        throw InvalidSeg("Internal error - segment computation problem (end)");

    segment_start_ = new_seg_start;
    segment_length_ = new_seg_end - new_seg_start;
    segment_state_ = sd.segment_state;
    allelic_state_at_Segment_end = sd.allele_state;
    chrom_ = 1; // dummy value

    if ( current_buf_index_ % 1000 == 1 ) {
        cout << "\r" << " Particle filtering " << setw(4) << int((segment_end() * 100) / seqlen_)
             << "% completed." << flush;
    }

    set_lookahead();
    
    current_buf_index_++;
}


void Segment::set_lookahead() {
    /* sets variables used by the Auxiliary Particle Filter implementation */
    first_singleton_distance.clear();
    first_singleton_distance.resize( nsam_, 0.0 );
    relative_mutation_rate.clear();
    relative_mutation_rate.resize( nsam_, 0.0 );
    doubleton.clear();
    first_split_distance = -1;

    /* set distance to mutation on left */
    distance_to_mutation_ = 0;
    for (int i = current_buf_index_; i >= 0; i-- ) {
        if (!buffer[i].all_alleles_missing()) {
            break;
        }
        distance_to_mutation_ = buffer[current_buf_index_].segment_start - buffer[i].segment_start;
    }

    /* hardcoded threshold of total length of missing data; beyond this no mutation will be recorded */
    double max_missing_data = 2000000;

    vector<bool> found_doubleton( nsam_ );     // records which sequences have participated in a doubleton
    int num_singletons = 0;
    int num_unphased_singletons = 0;
    int num_doubleton_sequences = 0;
    double total_length_times_branches = 0.1;  // used to calculate weighing for mutation rate due to missing data
    double total_length_times_branches_missing = 0.1;
    double total_current_missing = 0.0;        // current streak of missingness in any lineage
    double last_singleton_distance = 0.0;      // position of last singleton used; for bailing out the doubleton loop

    double distance = 0.0;
    for ( size_t i = current_buf_index_; i < buffer.size(); i++ ) {
        // update distance to mutation
        if (!buffer[i].all_alleles_missing()) {
            distance_to_mutation_ = min( distance_to_mutation_,
                                         (double)(buffer[i].segment_start + buffer[i].segment_length 
                                                  - buffer[current_buf_index_].segment_start) );
        }
        // is mutation at  i  a singleton, doubleton, or something else?
        // also identify the samples; take the first possible in case of unphased hets
        int num_var = 0, num_missing = 0;
        int s1 = -1, s2 = -1;
	is_singleton_unphased.clear();
        for ( size_t j = 0; j < nsam_ ; j++ ) {
	    is_singleton_unphased.push_back( false );
            if (buffer[i].allele_state[j] > 0) {
                num_var++;  // count the mutation
                if (num_var==1) s1 = j;
                if (num_var==2) s2 = j;
                if (buffer[i].allele_state[j] == 2) {
		    is_singleton_unphased[j] = true;
		    is_singleton_unphased.push_back( true );
                    j++;    // skip next sample, in case of a diploid unphased het, to avoid counting het twice
		}
            }
            if (buffer[i].allele_state[j] == -1) {                                       // missing data
                num_missing++;
                if (num_missing == 1) total_current_missing += buffer[i].segment_length; // keep track of current streak of missing data

                // if stretch of missing data is too long, stop recording this lineage.
                // Mark singleton if necessary, with the negative of the distance to START of segment
                // Also mark sequence as participating (i.e. not able to participate) in doubletons
                if (total_current_missing > max_missing_data) {
                    if (first_singleton_distance[j] == 0) {
                        // no singleton yet seen in this lineage; stop looking and mark as "no mutation seen" (negative distance)
                        // however, if most of the segment was missing, ignore (mark as minus small distance)
                        const double epsilon = 1e-6;
                        first_singleton_distance[j] = -(buffer[i].segment_start - buffer[current_buf_index_].segment_start) - epsilon;
                        last_singleton_distance = -first_singleton_distance[j];
                        if (first_singleton_distance[j] < 0.5 * total_current_missing) {
                            first_singleton_distance[j] = -epsilon;
                        }
                        relative_mutation_rate[j] = total_length_times_branches_missing / total_length_times_branches;
                        num_singletons++;
                    }
                    if (!found_doubleton[j]) {
                        // no doubleton yet found that uses this lineage; stop looking
                        found_doubleton[j] = true;
                        num_doubleton_sequences++;
                    }
                }
            }
        }
        // reset missing streak if warranted
        if (num_missing == 0) 
            total_current_missing = 0.0;
        // update total length, taking account of missing data
        total_length_times_branches         += buffer[i].segment_length * nsam_;
        total_length_times_branches_missing += buffer[i].segment_length * (nsam_ - num_missing);
        // if a long streak of missing data is encountered, stop since we then are unsure of existence of singletons and doubletons
        if (total_current_missing > max_missing_data)
            continue;

        double have_doubleton = false;
        distance = buffer[i].segment_start + buffer[i].segment_length - buffer[current_buf_index_].segment_start;
        if (num_var == 1) { // a singleton
            if (first_singleton_distance[s1] == 0) {
                first_singleton_distance[s1] = distance;
                relative_mutation_rate[s1] = total_length_times_branches_missing / total_length_times_branches;
                num_singletons++;
                last_singleton_distance = first_singleton_distance[s1];
		if (is_singleton_unphased[s1]) {
		  first_singleton_distance[s1 + 1] = distance;
		  relative_mutation_rate[s1 + 1] = relative_mutation_rate[s1];
		  num_singletons++;   // actually, number of alleles for which singleton information is obtained
		  num_unphased_singletons++;
		}
            }
        } else { // not a singleton
            // now loop over doubletons, and see (i) if we already have this doubleton (if one),
            // and (ii) which existing doubletons are incompatible with this allele
            for ( Doubleton& d : doubleton ) {
                if ( ( (d.seq_idx_1 | 1) == d.seq_idx_2 &&         // cherry involves one individual
                        buffer[i].allele_state[d.seq_idx_1] == 2)  // the non-singleton mutation has a het at the individual
                     ||
                     ( (buffer[i].allele_state[d.seq_idx_1] + buffer[i].allele_state[d.seq_idx_2] == 1) &&     // condition tests
                       ((buffer[i].allele_state[d.seq_idx_1] | buffer[i].allele_state[d.seq_idx_2]) == 1))) {  // for 0,1 or 1,0
                    // either a diploid individual carrying the cherry has an unphased het which is not a singleton,
                    // or the two individuals have differing genotypes neither of which are hets -- an incompatibility
                    d.incompatible = true;
                }
                // check if mutation is a doubleton and we've already seen it ( (@) see below )
                if (num_var == 2 && d.seq_idx_1 == s1 && d.seq_idx_2 == s2) {
                    have_doubleton = true;
                    // if no incompatible evidence for existing doubleton was found yet,
                    // this is the latest evidence for presence of the doubleton we know
                    if (!d.incompatible) {
                        d.last_evidence_distance = distance;
                    }
                }
            }
        }
        // enter new doubleton, if any, if we haven't entered it already, and if it is compatible with earlier doubletons
	if (num_var == 2 && have_doubleton == false && buffer[i].allele_state[s1] > -1 && buffer[i].allele_state[s2] > -1) {
            // in the case of unphased individuals, see if we can make this doubleton compatible by using the
            // other allele of the diploid individual
            for (int d1 = 0 ; d1 <= (buffer[i].allele_state[s1] == 2) ; d1++) {
                for (int d2 = 0 ; d2 <= (buffer[i].allele_state[s2] == 2) ; d2++) {
                    if (!found_doubleton[s1+d1] && !found_doubleton[s2+d2]) {
                        // use the original s1 s2 indices to mark the doubleton, so that we can match on the original indices
                        // at (@) above
                        doubleton.push_back( Doubleton( s1, s2, distance, distance,
                                                        buffer[i].allele_state[s1] == 2,
                                                        buffer[i].allele_state[s2] == 2 ) );
                        // mark sequences participating in doubletons, and count them
                        found_doubleton[s1+d1] = true; num_doubleton_sequences++;
                        found_doubleton[s2+d2] = true; num_doubleton_sequences++;
                        // bail out
                        d1 = 1; d2 = 1;
                    }
                }
            }
        }
        // now consider splits
        if (first_split_distance == -1 && num_var > 2 && (int)nsam_ - num_var > 2) {
            first_split_distance = distance;
            allelic_state_at_first_split = buffer[i].allele_state;
            mutation_count_at_first_split = min( num_var, (int)nsam_ - num_var );
        }
        // bail out if we're done
        if ((num_singletons == (int)nsam_) && num_doubleton_sequences >= (int)nsam_ - 1) break;
        // also bail out if way beyond last singleton, to avoid very long searches (it's a quadratic algo after all!)
        if ((num_singletons == (int)nsam_) && distance > (2+num_unphased_singletons)*last_singleton_distance) break;
    }
    // If we ran off the input sequence, fill in any missing singleton entries
    if (num_singletons < (int)nsam_) {
        for ( size_t j = 0; j < nsam_ ; j++ ) {
            if (first_singleton_distance[j] == 0) {
                first_singleton_distance[j] = -distance;
                relative_mutation_rate[j] = total_length_times_branches_missing / total_length_times_branches;
            }
        }
    }
    // Debug output
    if (0) {
        cout << current_buf_index_ << " " << buffer[current_buf_index_].segment_start << " ";
        for ( size_t j = 0; j < nsam_ ; j++ ) cout << buffer[current_buf_index_].allele_state[j];
        for ( size_t j = 0; j < nsam_ ; j++ ) cout << " S:" << first_singleton_distance[j];
        for ( Doubleton d : doubleton )
            cout << " D:" << d.seq_idx_1 << "," << d.seq_idx_2 << ":"
                 << d.first_evidence_distance << "-" << d.last_evidence_distance;
        cout << endl;
        for ( size_t j = 0; j < nsam_ ; j++ ) cout << " M:" << relative_mutation_rate[j]; cout << endl;
    }
    // cout << current_buf_index_ << " " << buffer[current_buf_index_].segment_start << " dist=" << distance_to_mutation() << endl;
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
        if ( seg_content == -1 && (i % 2) == 1 ) {
            if (!(allelic_state_at_seg_end[i-1] == -1)) {
                throw InvalidSeg("Found inconsistent unphased heterozygous marks");
            }
        }
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


bool Segment::is_repeated_allele() const {
    if (current_buf_index_ == 0 || buffer[current_buf_index_ - 1].all_alleles_missing()) 
        return false;
    int idx = (int)current_buf_index_ - 2;
    while (idx > 0 && buffer[idx].all_alleles_missing()) 
        --idx;
    for (int j=0; j<(int)buffer[idx].allele_state.size(); j++) {
        if (buffer[idx].allele_state[j] != buffer[current_buf_index_ - 1].allele_state[j])
            return false;
    }
    return true;
}

