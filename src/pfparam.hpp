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
//#include "usage.hpp"
#include "vcf.hpp"

using namespace std;

#ifndef PFARGPfParam
#define PFARGPfParam


class PfParam{
    public:
        
        /*!
         * Constructors and Destructors
         */              

        PfParam(int argc, char *argv[]);            
        ~PfParam();
        
        /*!
         * Methods
         */ 
        
        //int  log(Model *model, size_t random_seed, pfTime * runningtime, double inferred_recomb_rate);
        int  log(Model *model, size_t random_seed, double inferred_recomb_rate);
        void log_param(Model* model ,size_t random_seed, double inferred_recomb_rate);
        //void log_end(pfTime * running_time);
        void appending_Ne_file(Model *model, bool hist = false);

        void print_help();
        void print_option();
        void print_example();

        /*!
         * Members
         */
        
        int    default_nsam;
        double default_mut_rate; 
        double default_recomb_rate; 
        double default_loci_length;
        // ------------------------------------------------------------------
        // PfParameters 
        // ------------------------------------------------------------------                         
        size_t N; /*!< number of particles */
        double ESS; 
        double ESSthreshold; /*!< Effective sample size respect to the number of particles = ESS * N , 0 < ESS < 1 */
        bool   ESS_default_bool;
        string pattern;     /*! population segement pattern */
        double top_t;
        string scrm_input;
        bool   EM_bool;
        int    EM_steps;
        // ------------------------------------------------------------------
        // Input 
        // ------------------------------------------------------------------            
        string vcf_NAME;
        int buff_length;
        // ------------------------------------------------------------------
        // Action 
        // ------------------------------------------------------------------
        double lag;            
        bool online_bool;
        // ------------------------------------------------------------------
        // Output 
        // ------------------------------------------------------------------
        bool log_bool;            
        bool hist_bool;
        bool heat_bool;
        string out_NAME_prefix;            
        string HIST_NAME;
        string Ne_NAME;
        string log_NAME;
        string TMRCA_NAME;
        string WEIGHT_NAME;
        string BL_NAME;            
        double window;

        Vcf * VCFfile;
        //char * const*scrm_argv_;
        //char **scrm_argv_;
        //int scrm_argc_;
        
    private:
        void nextArg(std::string option);
        void init();
        void insert_mutation_rate_in_scrm_input ( );
        void insert_recomb_rate_and_seqlen_in_scrm_input ( );
        void insert_sample_size_in_scrm_input ( );
        void finalize_scrm_input ( );
        void finalize ( );
        void convert_scrm_input();
        const int argc_;
        int argc_i;
        char * const* argv_;            
    };

#endif
