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

#include "scrm/param.h"
#include "vcf.hpp"
//#include <boost/lexical_cast.hpp> 

using namespace std;

#ifndef PFARGPfParam
#define PFARGPfParam

/*!
 * \brief pf-ARG parameters 
 */ 
class PfParam{
    #ifdef UNITTEST
    //friend class TestVcf;
    friend class TestParam;
    #endif
    friend class ParticleContainer;
    public:
        
        //
        // Constructors and Destructors
        //              
        PfParam(int argc, char *argv[]);            
        ~PfParam();
        
        //
        // Methods
        //
        int  log( );
        void appending_Ne_file( bool hist = false );

        void print_help();
        void print_option();
        void print_example();

        //
        // Members
        //
        
        // ------------------------------------------------------------------
        // PfParameters 
        // ------------------------------------------------------------------                         
        size_t N;            /*!< \brief Number of particles */
        int    EM_steps;     /*!< \brief Number of EM iterations */
        double ESSthreshold; /*!< \brief Effective sample size, scaled by the number of particles = ESS * N , 0 < ESS < 1 */

        // ------------------------------------------------------------------
        // Action 
        // ------------------------------------------------------------------
        double lag;            
        bool online_bool;

        //bool hist_bool;
        bool heat_bool;
        //bool finite_bool;

        Vcf * VCFfile;
        Param *SCRMparam;
        Model *model;
        MersenneTwister *rg ;

        double default_loci_length;

        double ESS () const { return this-> ESS_;} // scaled between zero and one
        
    private:

        //
        // Methods
        //   
        
        void init();
        void insert_mutation_rate_in_scrm_input ( );
        void insert_recomb_rate_and_seqlen_in_scrm_input ( );
        void insert_sample_size_in_scrm_input ( );
        void finalize_scrm_input ( );
        void finalize ( );
        void convert_scrm_input();
        //void log_param( double inferred_recomb_rate, vector < vector<double> > migrate );
        void log_param( );


        void nextArg(){
            ++argc_i;            
            if (argc_i >= argc_) {
                throw std::invalid_argument(std::string("Not enough parameters when parsing options: ") + argv_[argc_i-1]);
                }
            }

        template<class T> T readNextInput() {
            this->nextArg();
            
            char c;
            T input;
            std::stringstream ss(argv_[argc_i]);
            ss >> input;
            if (ss.fail() || ss.get(c)) {
                throw std::invalid_argument(std::string("Failed to parse option: ") + argv_[argc_i-1]);
                }
            return input;
            }
            
        //
        // Members
        //  
        const int argc_;
        int argc_i;
        char * const* argv_;

        // ------------------------------------------------------------------
        // PfParameters 
        // ------------------------------------------------------------------                         
        bool   ESS_default_bool;
        string scrm_input;
        bool   EM_bool;
        double original_recombination_rate_;
        
        // ------------------------------------------------------------------
        // Default values 
        // ------------------------------------------------------------------        
        int    default_nsam;
        double default_mut_rate; 
        double default_recomb_rate; 
        size_t default_num_mut;

        // ------------------------------------------------------------------
        // Input 
        // ------------------------------------------------------------------
        string pattern;     /*! population segement pattern */
        double top_t;
        double ESS_;   // scaled between zero and one

        // VCF related
        string vcf_NAME; // vcf file name
        int buff_length; // number of lines vcf file read at once
        int ghost; // ghost snp, when no data is used, this is used for debugging
        int filter_window_;  // If two snps are two close, i.e. difference between the site is less than this window, should skip to the next read.
        int missing_data_threshold_;

        // ------------------------------------------------------------------
        // Output 
        // ------------------------------------------------------------------
        bool log_bool;            
        string out_NAME_prefix;            
        string HIST_NAME;
        string Ne_NAME;
        string log_NAME;
        string TMRCA_NAME;
        string WEIGHT_NAME;
        //string BL_NAME;
        string SURVIVOR_NAME;
        int heat_seq_window;
        

    };
    
//template<class T>
 //T readInput(char input[])
 //{
   //return boost::lexical_cast<T>(input);
 //}

#endif
