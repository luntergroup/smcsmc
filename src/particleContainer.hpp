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

#include"particle.hpp"
#include"pfparam.hpp"

#ifndef NDEBUG
#define resampledout (std::cout << "    RESAMPLE ")
#else
#pragma GCC diagnostic ignored "-Wunused-value"
#define resampledout 0 && (std::cout << "    RESAMPLE ")
#endif



#ifndef PARTICLECONTAINER
#define PARTICLECONTAINER

/*! \brief ForestState Container, particle filters*/
class ParticleContainer {

    friend class CountModel;

    public:
        //
        // Constructors and Destructors
        //
        ParticleContainer(Model* model,
                          MersenneTwister *rg,
                          const vector<int>& record_event_in_epoch,
                          size_t Num_of_states,
                          double initial_position,
                          bool heat_bool,
                          bool emptyFile,
                          vector <int> first_allelic_state);  // this is used to create the particle initial states
        ~ParticleContainer(); //Destuctor needs to free memory of each particles that the pointers are pointing to...

        //
        // Methods
        //
        void update_state_to_data( double mutation_rate, double loci_length, Segment * Segfile, valarray<double> & weight_cum_sum);
        void extend_ARGs( double mutation_rate, double extend_to, Segment_State segment_state );
        void set_particles_with_random_weight();
        void ESS_resampling(valarray<double> weight_cum_sum, valarray<int> &sample_count, int mutation_at, double ESSthreshold, int num_state);
        bool appendingStuffToFile(double x_end, PfParam &pfparam);
        void cumulate_recomb_opportunity_at_seq_end( double seqend );
        void normalize_probability();
        void clear();
        void print_particle_probabilities();
        void print_ln_normalization_factor();

        //
        // Debugging tools
        //
        void print();
        bool check_state_orders();
        void print_particle_newick();

    private:
        //
        // Methods
        //
        void update_weight_at_site( double mutation_rate, const vector <int> &data_at_tips);
        bool next_haplotype( vector<int>& haplotype_at_tips, const vector<int>& data_at_tips ) const;
        void store_normalization_factor();
        // Resampling
        void resample(valarray<int> & sample_count);
        void duplicate_particles ( valarray<int> & sample_count );
        void resample_for_check(valarray<int> & sample_count);
        void shifting(int number_of_particles);
        void trivial_resampling( std::valarray<int> & sample_count, size_t num_state );
        void systematic_resampling(std::valarray<double> cum_sum, std::valarray<int>& sample_count, int sample_size);
        void update_cum_sum_array_find_ESS(std::valarray<double> & weight_cum_sum);

        //
        // Setters and getters:
        //
        double ESS() const {return this->ESS_;};
        void set_ESS(double ess){this->ESS_ = ess;};
        RandomGenerator* random_generator() const { return this->random_generator_; }
        double current_printing_base() const { return this->current_printing_base_;}
        void set_current_printing_base (double base) { this->current_printing_base_ = base;}
        double ln_normalization_factor() const { return this->ln_normalization_factor_; }

        //
        // Members
        //
        vector <ForestState*> particles;
        double ESS_;
        RandomGenerator* random_generator_; // This is for particle filter only,
        double current_printing_base_;
        bool heat_bool_ ;
        double ln_normalization_factor_;
        double temp_sum_of_weights;
};

#endif
