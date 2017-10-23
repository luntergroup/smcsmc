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
        ParticleContainer(Model *model,
                          MersenneTwister *rg,
                          PfParam *pfparam,
                          size_t Num_of_states,
                          double initial_position,
                          bool emptyFile,
                          vector <int> first_allelic_state);  // this is used to create the particle initial states
        ~ParticleContainer(); //Destuctor needs to free memory of each particles that the pointers are pointing to...

        //
        // Methods
        //
        void update_state_to_data( Segment * Segfile, bool ancestral_aware, const TerminalBranchLengthQuantiles& term_branch_lengths );
        void extend_ARGs( double extend_to, const vector<int>& data_at_site );

        void extend_ARGs_importance_sampling( double extend_to, const vector<int>& data_at_site, const Segment& segment,
                                              bool ancestral_aware, const TerminalBranchLengthQuantiles& terminal_branch_lengths );
        void set_particles_with_random_weight();
        int  resample(int update_to, const PfParam &pfparam, vector<ForestState*>* particles = NULL, int to_sample = -1);
        void normalize_probability();
        void clear();
        double ln_normalization_factor() const { return this->ln_normalization_factor_; }

        //
        // Debugging tools
        //
        void print();
        bool check_state_orders();
        void print_particle_newick();
        void print_recent_recombination_histogram();

    private:
        //
        // Methods
        //
        void update_data_status_at_leaf_nodes( const vector<int>& data_at_tips );
        int calculate_initial_haplotype_configuration( const vector<int>& data_at_tips, vector<int>& haplotype_at_tips ) const;
        void update_weight_at_site( const Segment &segment,
                                    const vector<ForestState*>& particles,
                                    bool ancestral_aware,
                                    const TerminalBranchLengthQuantiles& term_brach_len);
        void update_lookahead_likelihood( const Segment& segment,
                                          const vector<ForestState*>& particles,
                                          bool ancestral_aware,
                                          const TerminalBranchLengthQuantiles& terminal_branch_lengths );
        bool next_haplotype( vector<int>& haplotype_at_tips, const vector<int>& data_at_tips ) const;
        // Resampling
        vector<ForestState*> implement_resampling( const vector<ForestState*>& particles, valarray<int> & sample_count, double total_pilot_weight) const;
        void systematic_resampling( valarray<double>& partial_sum, valarray<int>& sample_count, int N) const;

        //
        // Setters and getters:
        //
        RandomGenerator* random_generator() const { return this->random_generator_; }

        //
        // Members
        //
        vector <ForestState*> particles;
        int num_particles;
        RandomGenerator* random_generator_;
        double ln_normalization_factor_;
        Model* model;
};

#endif
