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
#include<boost/assign/std/vector.hpp>
#include<cassert>       // assert
#include<stdexcept>      // std::invalid_argument
#include<stdlib.h>     /* strtol, strtod */
using namespace std;

#ifndef NDEBUG
#define dout std::cout
#else
#pragma GCC diagnostic ignored "-Wunused-value"
#pragma clang diagnostic ignored "-Wdeprecated-register"
#define dout 0 && std::cout
#endif

#ifndef VCF
#define VCF

/*! @ingroup group_data */
class Vcf{
    friend class PfParam;
    
    #ifdef UNITTEST
    friend class TestVcf;
    friend class TestParam;
    #endif

    public:
    
        /*!
         * Constructors and Destructors
         */          
        Vcf(){};
        Vcf(string file_name, int buffer_length = 100);
        ~Vcf(){};        
        
        /*!
         * Methods
         */ 
        void read_new_line();
        void reset_VCF_to_data();
        void force_to_end_data() { this->end_data_ = true; }
        
        // DEBUG
        void print_vcf_line(vector<string> sample_names);
        void print();
        
        /*!
         *  Getters:
         */ 
        double site() const { return this->site_; }
        bool withdata() const { return this->withdata_; }
        bool end_data() const { return this->end_data_; }

        /*!
         * Members
         */       
        size_t vcf_file_length;
        vector <string> sample_names;
        size_t header_end_line;        
        int buffer_max_number_of_lines;

        vector <string> buffer_lines;
        
        bool skip;
        string ref;
        vector <string> alt;
        vector <string> vec_of_sample_alt;
        vector <bool> vec_of_sample_alt_bool;
        vector < vector<string> > sample_alt;
        vector <bool> phased; // True if it is phased, which has '|'
            
    private:

        /*!
         * Methods
         */ 
        void empty_block();
        void read_new_block();
        void init(string infile_name, int buffer_length);

        /*!
         *  Setters and getters:
         */ 
        size_t nsam() const { return this->nsam_; }
 
        void set_nsam( size_t nsam ) { this->nsam_ = nsam; }
        size_t nfield() const { return this->nfield_; }    
        int chrom() const { return this ->chrom_; }

        bool eof() const { return this->eof_; }
        double even_interval() const { return this-> even_interval_; }
        void set_even_interval( double interval ) { this->even_interval_ = interval; }
                 
        void check_feilds(string line);
        string extract_alt_(string tmp_str, size_t start, size_t end);
        bool print_sample_name();
        /*!
         * Members
         */ 
        bool eof_;
        string file_name_;
        size_t current_line_index_;
        bool withdata_;
        size_t vcf_length_;
        bool end_data_;
        size_t end_pos_;
        double even_interval_; 
        
        //HEADER
        size_t header_end_pos_;
        size_t nsam_;
        size_t nfield_;
        
        //VCFBODY
        size_t current_block_line_;
        size_t empty_file_line_counter_;
        
        //VCFLINE    
        double site_;
        double previous_site_at_;
        int chrom_;
        int pervious_chrom_;
    };

#endif
