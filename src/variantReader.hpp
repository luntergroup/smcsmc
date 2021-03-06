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

#ifndef VARIANTREADER
#define VARIANTREADER


enum INPUT_FILETYPE {EMPTY, VCF, GVCF, RGVCF};
enum Seq_State { ZERO_SEG,      /* point; may be variant */
                 SEQ_INVARIANT, /* segment; invariant */
                 MISSING  };    /* segment; no data */ 
enum Variant_State {SNP, INVARIANT, OTHER_VARIANT};


/*!
 * P1 ----------- P2 ___________ P3
 * 
 * ----: INVARIANT
 * ____: MISSING
 * 
 * Let P denote Variant position, at position P, the state can be a SNP, INVARIANT or OTHER_VARIANT(indel, insertion, or deletion)
 * Between two positions P1 and P2, the sequence state (from position P1 to P2) can be missing data, or INVARIANT
 */ 
class VariantPosition{

    public:
        int chrom() const { return this ->chrom_; }
        int site() const { return this->site_; } 
        vector <int> int_vec_of_sample_alt;
        size_t nsam() const { return this->nsam_; } 

    protected:
        
        VariantPosition(){};
        Variant_State current_variant_state;
        Variant_State previous_variant_state;
        void set_nsam( size_t nsam ) { this->nsam_ = nsam; }

        //all these numbers can not be negative
        int pervious_chrom_;
        int previous_site_at_;
        
        int site_;
        int chrom_;
                
        size_t nsam_; // number of diploid sample

        string ref;
        vector <string> alt;
        vector < vector<string> > sample_alt;
        vector <string> sample_names;

        vector <string> vec_of_sample_alt;
        vector <bool> phased; // True if it is phased, which has '|'
        
};


class VariantSegment: public VariantPosition{
        
    protected:
        Seq_State previous_seg_state;        
        Seq_State current_seg_state;
        void reset_chrom_site(){ 
            this->chrom_ = 0; 
            this->site_ = 0; 
            this->pervious_chrom_ = 0; 
            this->previous_site_at_ = 0; 
            this->previous_seg_state = ZERO_SEG;
            this->previous_variant_state = INVARIANT;
            }
        int seg_end_site() const { return this->seg_end_site_; } 
        int seg_end_site_;
};


/*! \brief VCF file reader @ingroup group_data */
class VariantReader: public VariantSegment{
    friend class PfParam;
    friend class ParticleContainer;
    #ifdef UNITTEST
    friend class TestVCF;
    friend class TestGVCF;
    friend class TestParam;
    #endif

    public:
        // Constructors and Destructors
        VariantReader(string file_name, INPUT_FILETYPE FileType_in = EMPTY, int buffer_length = 100);
        ~VariantReader(){};
        
        // Methods
        void read_new_line();
        void reset_data_to_first_entry();
        void force_to_end_data() { this->end_data_ = true; }
        
        // DEBUG
        void print_vcf_line();
        void print();

        bool end_data() const { return this->end_data_; }
        INPUT_FILETYPE FileType;

    private:

        // Methods
        void empty_block();
        void read_new_block();
        void init();
        
        //  Setters and getters:
        //size_t nfield() const { return this->nfield_; }    

        //bool eof() const { return this->eof_; }
        //int even_interval() const { return this-> even_interval_; }
        void set_even_interval( int interval ) { this->even_interval_ = interval; }
        
        
        string tmp_line;
        string tmp_str;
        void skip_Header( );
        void check_feilds( );
        
        void set_empty_line_entry();
        void initialize_read_newLine();                
        void check_and_update_block();
        void check_and_update_newLine();
        void finalize_read_new_line();
        
        void extract_field_CHROM   ( );
        void extract_field_POS     ( );
        void extract_field_ID      ( );
        void extract_field_REF     ( );
        void extract_field_ALT     ( );
        void extract_field_QUAL    ( );
        void extract_field_FILTER  ( );
        void extract_field_INFO    ( );
        void extract_field_FORMAT  ( );
        void extract_field_VARIANT ( );
        
        string extract_field_ALT_str( size_t start, size_t end );


        
        bool print_sample_name();
        
        // Members

        // Line related
        size_t feild_start;
        size_t field_end;
        int field_index;

        
        // FILE related
        bool eof_;
        string file_name_;
        size_t current_line_index_; // line counter in the entire vcf file
        size_t file_length_;
        bool end_data_;
        size_t end_pos_;
        int buffer_max_number_of_lines;
        
        // Header related
        size_t header_end_pos_;
        //size_t nfield_;
        size_t header_end_line;        
        
        // Block related
        size_t current_block_line_;  // line counter in the current data block
        size_t empty_file_line_counter_;

        size_t ghost_num_mut;
        int filter_window_;  // If two snps are too close, i.e. difference between the site is less than this window, should skip to the next read.
        int missing_data_threshold_; // if two snps are too far away apart, treat as missing data
        int even_interval_; 
        
        //size_t vcf_file_length;

        vector <string> buffer_lines;
        
        bool skip_tmp_line;
    };

#endif
