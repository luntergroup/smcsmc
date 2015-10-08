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

#include "mersenne_twister.h"
#include "param.h"
#include "segdata.hpp"

using namespace std;

#ifndef PFARGPfParam
#define PFARGPfParam

/*!
 * \brief smcsmc parameters
 */
class PfParam{
    #ifdef UNITTEST
    friend class TestParam;
    #endif
    friend class ParticleContainer;
    public:

        // Constructors and Destructors
        PfParam(int argc, char *argv[]);
        ~PfParam();

        //
        // Methods
        //
        int  log( );
        void appending_Ne_file( bool hist = false );
        void append_to_count_file( size_t epoch, string label, int from_pop, int to_pop, double opportunity, double count, double weight );

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

        bool heat_bool;

        Segment * Segfile;
        Param *SCRMparam;
        Model model;
        MersenneTwister *rg ;
        vector<int> record_event_in_epoch;

        static const int RECORD_RECOMB_EVENT = 1;
        static const int RECORD_COALMIGR_EVENT = 2;

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

        int readRange(int& second) {
			this->nextArg();
			char c;
			double input;
			std::stringstream ss(argv_[argc_i]);
			ss >> input;
			second = input;
			if (ss.fail() || input < 1) throw std::invalid_argument(std::string("Failed to parse int (e.g. '1') or int range (e.g. '1-10'): ") + argv_[argc_i]);
			if (!ss.get(c)) return input;
			if (c != '-') throw std::invalid_argument(std::string("Parsing failure: expected '-' after int in range expression: ") + argv_[argc_i]);
			ss >> second;
			if (ss.fail() || second <= input || ss.get(c)) throw std::invalid_argument(std::string("Failed to parse int or int range: ") + argv_[argc_i]);
			return input;
		}

        // Members
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
        size_t default_nsam;
        double default_mut_rate;
        double default_recomb_rate;
        double default_num_mut;

        // ------------------------------------------------------------------
        // Input
        // ------------------------------------------------------------------
        string pattern;     /*! population segement pattern */
        double top_t_;
        double ESS_;   // scaled between zero and one

        // Segdata related
        string input_SegmentDataFileName;

        // ------------------------------------------------------------------
        // Output Related
        // ------------------------------------------------------------------
        bool log_bool;
        bool rescue_bool;
        string out_NAME_prefix;
        string HIST_NAME;
        string Ne_NAME;
        string Count_NAME;
        string log_NAME;
        string TMRCA_NAME;
        string WEIGHT_NAME;
        //string BL_NAME;
        string SURVIVOR_NAME;
        int heat_seq_window;

        // ------------------------------------------------------------------
        // Version Info
        // ------------------------------------------------------------------
        string smcsmcVersion;
        string scrmVersion;
        string compileTime;
};

#endif
