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

#ifndef VARIANTREADER
#define VARIANTREADER



class Variant{
    friend class VariantReader;

    public:
        int chrom() const { return this ->chrom_; }
        int site() const { return this->site_; } 
        vector <bool> vec_of_sample_alt_bool;
        size_t nsam() const { return this->nsam_; } 
        bool invariant() const { return this->invariant_; }

    private:

        Variant(){};

        void set_nsam( size_t nsam ) { this->nsam_ = nsam; }

        size_t nsam_;

        int site_;
        int chrom_;

        string ref;
        vector <string> alt;
        vector < vector<string> > sample_alt;
        vector <string> sample_names;

        vector <string> vec_of_sample_alt;
        vector <bool> phased; // True if it is phased, which has '|'
        
        bool invariant_;
    };

enum INPUT_FILETYPE {EMPTY, VCF, GVCF, RGVCF};

/*! \brief VCF file reader @ingroup group_data */
class VariantReader: public Variant{
    friend class PfParam;
    friend class ParticleContainer;
    #ifdef UNITTEST
    friend class TestVariantReader;
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
        
        //
        // Getters:
        //
        //bool withdata() const { return this->withdata_; }
        void set_missding_data ( bool TRUEorFALSE ) { this->missing_data_ = TRUEorFALSE ; }
        bool missing_data() const { return this->missing_data_; }
        //bool empty_file() const { return this->empty_file_; }

        bool end_data() const { return this->end_data_; }

    private:

        //
        // Methods
        //   
        void empty_block();
        void read_new_block();
        void init(string infile_name, int buffer_length, INPUT_FILETYPE FileType_in);

        //  Setters and getters:
        size_t nfield() const { return this->nfield_; }    

        bool eof() const { return this->eof_; }
        int even_interval() const { return this-> even_interval_; }
        void set_even_interval( int interval ) { this->even_interval_ = interval; }
                 
        void check_feilds(string line);
        string extract_alt_(string tmp_str, size_t start, size_t end);
        bool print_sample_name();
        
        //
        // Members
        //   
        
        INPUT_FILETYPE FileType;
        
        
        bool eof_;
        string file_name_;
        size_t current_line_index_; // line counter in the entire vcf file
        //bool withdata_;
        bool missing_data_;
        //bool empty_file_;
        size_t vcf_length_;
        bool end_data_;
        size_t end_pos_;
        
        
        //HEADER
        size_t header_end_pos_;
        size_t nfield_;
        
        //VCFBODY
        size_t current_block_line_;  // line counter in the current data block
        size_t empty_file_line_counter_;
        
        
        //all these numbers can not be negative
        int pervious_chrom_;
        int previous_site_at_;
        int ghost_num_mut;
        int filter_window_;  // If two snps are too close, i.e. difference between the site is less than this window, should skip to the next read.
        int missing_data_threshold_; // if two snps are too far away apart, treat as missing data
        int even_interval_; 
        
        //size_t vcf_file_length;
        size_t header_end_line;        
        int buffer_max_number_of_lines;

        vector <string> buffer_lines;
        
        bool skip;

    };

#endif
