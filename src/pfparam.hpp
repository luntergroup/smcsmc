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
        //int  log( double inferred_recomb_rate, vector < vector<double> > migrate );
        //void appending_Ne_file(Model *model, bool hist = false);
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
        double ESSthreshold; /*!< \brief Effective sample size respect to the number of particles = ESS * N , 0 < ESS < 1 */

        // ------------------------------------------------------------------
        // Action 
        // ------------------------------------------------------------------
        double lag;            
        bool online_bool;

        //bool hist_bool;
        bool mid_bool;
        bool heat_bool;
        bool finite_bool;

        Vcf * VCFfile;
        Param *SCRMparam;
        Model *model;
        MersenneTwister *rg ;

        double ESS () const { return this-> ESS_;}
    private:

        // ------------------------------------------------------------------
        // PfParameters 
        // ------------------------------------------------------------------                         
        double ESS_; 
        bool   ESS_default_bool;
        double top_t;
        string scrm_input;
        bool   EM_bool;
        double original_recombination_rate_;
        //
        // Methods
        //   
        void nextArg(std::string option);
        void init();
        void insert_mutation_rate_in_scrm_input ( );
        void insert_recomb_rate_and_seqlen_in_scrm_input ( );
        void insert_sample_size_in_scrm_input ( );
        void finalize_scrm_input ( );
        void finalize ( );
        void convert_scrm_input();
        //void log_param( double inferred_recomb_rate, vector < vector<double> > migrate );
        void log_param( );
        
        //
        // Members
        //  
        const int argc_;
        int argc_i;
        char * const* argv_;
        
        // ------------------------------------------------------------------
        // Default values 
        // ------------------------------------------------------------------        
        int    default_nsam;
        double default_mut_rate; 
        double default_recomb_rate; 
        double default_loci_length;
        // ------------------------------------------------------------------
        // Input 
        // ------------------------------------------------------------------
        string vcf_NAME;
        int buff_length;
        string pattern;     /*! population segement pattern */

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
        double window;

    };

#endif
