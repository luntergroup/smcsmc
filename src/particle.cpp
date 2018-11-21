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

#include <climits> // INT_MAX
#include "particle.hpp"
#include "pfparam.hpp"
#include "descendants.hpp"


inline double fastexp(double x) {
    // Based on a generalized continued fraction.  Just two divisions, two multiplications.
    // If |x|<0.71844823, the relative error is < 1e-6
    // See wikipedia Exponential_function, and L.Lorentzen and H. Waadeland, Continued Fractions, Atlantis Studies in Maths pg 268
    double xx = x*x;
    if (xx < 0.516167859) {
        return 1 + 2*x / (2 - x + xx / (6 + xx * 0.1));
    } else {
        return exp(x);
    }
}


inline double fastexp_approx(double x) {
    // Based on a generalized continued fraction.  Just one division, two multiplications.
    // The approximate branch is used only for |x|<1.44885, where the relative error is < 1e-2
    // See wikipedia Exponential_function, and L.Lorentzen and H. Waadeland, Continued Fractions, Atlantis Studies in Maths pg 268
    double xx = x*x;
    if (xx < 2.099166) {
        return 1 + 2*x / (2 - x + xx * (1.0/6));
    } else {
        return exp(x);
    }
}


inline double nchoosek(int n, int k) {
    double nCk = 1.0;
    for (int i=1; i<=k; i++) {
        nCk *= (n-i+1)/double(i);
    }
    return nCk;
}


inline double exp_digamma( double x ) {
  if (x > 10) return x - 0.5 + (x+0.5)/(24*x*x);
  double f = 0.0;
  while (x < 6) {
    f = f + 1.0/x; // psi(x) = psi(x+1) - 1/x
    x = x + 1.0;
  }
  double psi = log(x) - 1/(2*x) - 1/(12*x*x);
  return exp(psi - f);
}


/*! \brief Initialize a new ForestState.  A copy of model/random_generator is made if own_m_and_rg==True; either way the objects are not stolen
 * @ingroup group_pf_init
 * */
ForestState::ForestState( Model* model,
                          RandomGenerator* random_generator,
                          PfParam* pfparam,
                          bool own_model_and_random_generator)
            : Forest( model, random_generator ),
              lookahead_weight_(1.0),
              pfparam( *pfparam )
              // , recent_recombination_count(0) 
            {

    /*! Initialize base of a new ForestState, then do nothing, other members will be initialized at an upper level */
    setPosteriorWeight( 1.0 );
    setPilotWeight( 1.0 );
    setMultiplicity( 1 );
    total_delayed_adjustment_ = 1.0;
    setSiteWhereWeightWasUpdated(0.0);
    first_event_height_ = -1.0;
    first_coal_height_ = -1.0;
    last_coal_height_ = -1.0;
    owning_model_and_random_generator = own_model_and_random_generator;
    if (owning_model_and_random_generator) {
        // as we're owning it, we should make copies
        size_t new_seed = (size_t)random_generator_->sampleInt( INT_MAX );
        random_generator_ = new MersenneTwister( new_seed , random_generator_->ff() );
        model_ = new Model( *model_ );
    }
    save_recomb_state();
}


/*! \brief Create a newly copied ForestState
    @ingroup group_pf_resample
*/
ForestState::ForestState( const ForestState & copied_state )
            :Forest( copied_state ),
             lookahead_weight_( copied_state.lookahead_weight_), 
             pfparam( copied_state.pfparam )
             // , recent_recombination_count( copied_state.recent_recombination_count) 
           {

    setPosteriorWeight( copied_state.posteriorWeight() );
    setPilotWeight( copied_state.pilotWeight() );
    total_delayed_adjustment_ = copied_state.total_delayed_adjustment_;
    delayed_adjustments = copied_state.delayed_adjustments;
    setSiteWhereWeightWasUpdated( copied_state.site_where_weight_was_updated() );
    copyEventContainers ( copied_state );
    current_rec_ = copied_state.current_rec_;
    owning_model_and_random_generator = copied_state.owning_model_and_random_generator;
    if (owning_model_and_random_generator) {
        // as we (and copied_state) own model and rg, we should make copies
        size_t new_seed = (size_t)random_generator_->sampleInt( INT_MAX );
        random_generator_ = new MersenneTwister( new_seed , random_generator_->ff() );
        model_ = new Model( *model_ );
    }
    _current_seq_idx = copied_state._current_seq_idx;
    setMultiplicity( 1 );
}


void ForestState::copyEventContainers(const ForestState & copied_state ) {
    // Copy trees
    for (size_t i=0; i < copied_state.eventTrees.size() ; i++) {
        EvolutionaryEvent* new_event = copied_state.eventTrees[i];
        eventTrees.push_back( new_event );
        if (new_event) {
            new_event->increase_refcount();
        }
    }
}


void ForestState::init_EventContainers( Model * model ) {
    for (size_t i=0; i < model->change_times_.size() + 1; i++) {
        eventTrees.push_back( NULL );
    }
}


/*! \brief Destructor of ForestState
 * Recursively remove all the previous states, if the pointer counter is zero
 */
ForestState::~ForestState() {
    clear_eventContainer();
    if (owning_model_and_random_generator) {
        delete random_generator_;
        delete model_;
        random_generator_ = NULL;
        model_ = NULL;
    }
    this->contemporaries_.clear();
    this->rec_bases_.clear();
}


/*! Clear coalescent and recombination events recorded between two states.*/
void ForestState::clear_eventContainer() {

    for (int i = eventTrees.size()-1 ; i>=0 ; --i) {
        if (eventTrees[i] && eventTrees[i]->decrease_refcount_is_zero()) {
            // We use placement new, so need to call destructor explicitly.
            // However the destructor should recursively delete its parents,
            // and therefore must know the epoch -- but we can't pass parameters.
            // So, call a helper deleter, that can take a parameter and both
            // destructs and deallocates the memory.
            eventTrees[i]->deletethis( i );  // this recursively deletes its parents
        }
    }
}


/*! Records recombination, coalescent and migration events as they occur while simulating a tree backwards in time.
 *  Also applies Variational Bayes likelihood factor
 */
void ForestState::record_all_event(TimeIntervalIterator const &ti) {

    double start_base, end_base;
    double start_height, end_height;
    size_t start_height_epoch, end_height_epoch;
    int record_mode = pfparam.record_event_in_epoch[ model().current_time_idx_ ];

    // find extent of time interval where an event may have occurred
    start_height = (*ti).start_height();
    start_height_epoch = ti.forest().model().current_time_idx_;
    if (this->tmp_event_.isNoEvent()) {
        end_height = start_height + (*ti).length();
    } else {
        end_height = this->tmp_event_.time();
        assert (end_height > start_height);
    }
    // interval either runs to event, or to change point; in both cases the end_height has the
    // same epoch as start_height (but call to getTimeIdx will return epoch+1 if end_height ran to end of interval)
    end_height_epoch = start_height_epoch;

    // loop over the two nodes
    for (int i=0; i<2; i++) {

        if (states_[i] == 2) {
            // node i is tracing an existing non-local branch; opportunities for recombination
            if (!(record_mode & (PfParam::RECORD_RECOMB_EVENT | PfParam::RECORD_TREE_EVENT))) continue;
            if (start_height_epoch > max_epoch_to_record_) continue;
            start_base = get_rec_base(active_node(i)->last_update());
            end_base = this->current_base();
            if (end_base == start_base) continue;
            double recomb_pos = -1;
            if (tmp_event_.isRecombination() && tmp_event_.active_node_nr() == i) {
                recomb_pos = (start_base + end_base)/2;   // we should sample from [start,end], but the data isn't used
            }
            // recording recombination events?  Then add event / opportunity in tree data structure
            if (record_mode & PfParam::RECORD_RECOMB_EVENT) {
                void* event_mem = Arena::allocate( start_height_epoch );
                EvolutionaryEvent* recomb_event = new(event_mem) EvolutionaryEvent( start_height, start_height_epoch,
                                                                                    end_height, end_height_epoch,
                                                                                    start_base, end_base, 1 );
                if (recomb_pos > -1) {
                    recomb_event->set_recomb_event_pos( recomb_pos );
                    recomb_event->set_descendants( get_descendants( active_nodes_[i] ) ); // set signature for node's descendants
                }
                recomb_event->add_leaf_to_tree( &eventTrees[ model().current_time_idx_] );
            }
            // recording coalescent-tree-modifying events?
            if ((record_mode & PfParam::RECORD_TREE_EVENT) && recomb_pos > -1) {
                int tree_epoch = eventTrees.size() - 1;            // pseudo-epoch in which tree-modifying events are stored
                void* event_mem = Arena::allocate( tree_epoch );   // allocate memory
                EvolutionaryEvent* tree_event = new(event_mem) EvolutionaryEvent( start_height, start_height_epoch,
                                                                                  end_height, end_height_epoch,
                                                                                  start_base, end_base, 1 );
                tree_event->set_recomb_event_pos( recomb_pos );
                tree_event->set_descendants( get_descendants( active_nodes_[i] ) ); // set signature for node's descendants
                tree_event->add_leaf_to_tree( &eventTrees[ tree_epoch ] );
            }
            assert(recomb_event->print_event());
        } else if (states_[i] == 1) {
            // node i is tracing out a new branch; opportunities for coalescences and migration
            // consider normal (not pairwise) coalescences that occurred on this node
            bool coal_event = (tmp_event_.isCoalescence() && tmp_event_.active_node_nr() == i);
            int weight = 0;
            // account for potential pairwise coalescence opportunity and event
            if (i==0 && states_[1]==1 && active_node(0)->population() == active_node(1)->population()) {
                weight = 1;
                coal_event |= tmp_event_.isPwCoalescence();
            }
            bool migr_event = (tmp_event_.isMigration() && tmp_event_.active_node_nr() == i);
            if (coal_event || migr_event) {
                if (first_event_height_ < 0.0)              first_event_height_ = end_height;
                if (coal_event && first_coal_height_ < 0.0) first_coal_height_ = end_height;
                last_coal_height_ = end_height;
                if (model().variational_bayes_correction_) {
                    double event_count =
                        coal_event ?
                        model().coalescent_event_count( active_node(i)->population() ) :
                        model().migration_event_count( active_node(i)->population(), tmp_event_.mig_pop() );
                    adjustWeights( exp_digamma( event_count ) / event_count );
                }
            }
            if (!(record_mode & PfParam::RECORD_COALMIGR_EVENT)) continue;
            if (start_height_epoch > max_epoch_to_record_) continue;
            start_base = this->current_base();
            weight += ti.contemporaries_->size( active_node(i)->population() );
            // Record coalescence and migration opportunity (note: if weight=0, there is still opportunity for migration)
            void* event_mem = Arena::allocate( start_height_epoch );
            EvolutionaryEvent* migrcoal_event = new(event_mem) EvolutionaryEvent(start_height, start_height_epoch,
                                                                                 end_height, end_height_epoch, start_base,
                                                                                 active_node(i)->population(), weight);
            if (coal_event) {
                migrcoal_event->set_coal_event();
            }
            if (migr_event) {
                migrcoal_event->set_migr_event( tmp_event_.mig_pop() );
            }
            // add event in tree data structure
            migrcoal_event->add_leaf_to_tree( &eventTrees[ model().current_time_idx_] );
            // add coalescent-tree-modifying event in data structure.  Note - must be recording coal/migr events!
            if ((record_mode & PfParam::RECORD_TREE_EVENT) && !migrcoal_event->is_no_event()) {
                int tree_epoch = eventTrees.size() - 1;      // pseudo-epoch in which tree-modifying events are stored
                event_mem = Arena::allocate( tree_epoch );   // allocate memory
                EvolutionaryEvent* tree_event = new(event_mem) EvolutionaryEvent( PfParam::RECORD_TREE_EVENT, *migrcoal_event );
		tree_event->set_descendants( get_descendants( active_node(i) ) );   // set signature for node's descendants
                tree_event->add_leaf_to_tree( &eventTrees[ tree_epoch ] );
            }
            assert(migrcoal_event->print_event());
        }
    }
}


void ForestState::record_recomb_extension (){

    // iterate over time intervals (but do NOT prune branches at this stage)
    for (TimeIntervalIterator ti(this, this->nodes_.at(0)); ti.good(); ++ti) {
        // Create a recombination event for this slice (which may be smaller than an epoch -- but in our case it usually won't be)
        int contemporaries = ti.contemporaries_->numberOfLocalContemporaries();
        if (contemporaries > 0 && 
            (pfparam.record_event_in_epoch[ model().current_time_idx_ ] & PfParam::RECORD_RECOMB_EVENT) &&
            model().current_time_idx_ <= max_epoch_to_record_ ) {
            double start_height = (*ti).start_height();
            double end_height = (*ti).end_height();
            size_t start_height_epoch = ti.forest()->model().current_time_idx_;
            //assert( start_height_epoch == ti.forest()->model().getTimeIdx( start_height ) );
            size_t end_height_epoch = start_height_epoch;

            // try and find recombination event that could be extended
            EvolutionaryEvent* event = eventTrees[ model().current_time_idx_];
            EvolutionaryEvent* previous = NULL;
            while (event && event->ref_counter() == 1 && event->get_end_base() >= current_base()) {
                if (event->is_recomb() &&
                    event->is_no_event() &&    // cannot extend if recomb event already set, to avoid setting twice
                    event->get_end_base() == current_base() &&
                    event->get_start_height() == start_height &&
                    event->get_end_height() == end_height &&
                    event->get_weight() == contemporaries) {

                    // extendable event
                    event->set_end_base( next_base() );
                    // make event first in chain, so assumption in resample_recombination_position is met (and current fn is slightly faster)
                    if (previous != NULL) {
                        // only do anything if not already first in chain. Note, all events have ref count 1 so simple linked list
                        previous->parent() = event->parent();
                        event->parent() = eventTrees[ model().current_time_idx_ ];
                        eventTrees[ model().current_time_idx_ ] = event;
                    }
                    goto dont_make_new;   // break out of while loop; then immediately continue for loop
                } else {
                    previous = event;
                    event = event->parent();
                }
            }
            
            // couldn't find extendable recombination record, so make new
            void* event_mem = Arena::allocate( start_height_epoch );
            // no event for now
            EvolutionaryEvent* recomb_event = new(event_mem) EvolutionaryEvent( start_height, start_height_epoch,
                                                                                end_height, end_height_epoch, current_base(),
                                                                                next_base(), contemporaries );
            recomb_event->add_leaf_to_tree( &eventTrees[ model().current_time_idx_] );
        }
    dont_make_new: ;
    }
}


void ForestState::record_recomb_event( double event_height )
{
    this->writable_model()->resetTime( event_height );
    size_t epoch_i = this->model().current_time_idx_;
    int record_mode = pfparam.record_event_in_epoch[ epoch_i ];
    if (epoch_i > max_epoch_to_record_) return;    /* skip events in long missing data segments */
    if (record_mode & PfParam::RECORD_RECOMB_EVENT) {
        // find the EvolutionaryEvent to add this event to.
        EvolutionaryEvent* event = eventTrees[ epoch_i ];
        while ( !event->is_recomb() || !event->recomb_event_overlaps_opportunity_t( event_height ) ) {
            event = event->parent();
            assert (event != NULL);
        }
        if (event->start_base() == this->current_base()) { /* to deal with events on edge of long missing data segments */
            event->set_recomb_event_time( event_height );
            event->set_descendants( get_descendants( rec_point.base_node() ) );  // calculate signature for descendants subtended by node
        }
    }
    // add coalescent-tree-modifying event in data structure.  Note - must be recording recombination events!
    if ((record_mode & PfParam::RECORD_TREE_EVENT)) {
        int tree_epoch = eventTrees.size() - 1;            // pseudo-epoch in which tree-modifying events are stored
        void* event_mem = Arena::allocate( tree_epoch );   // allocate memory
        //EvolutionaryEvent* tree_event = new(event_mem) EvolutionaryEvent( PfParam::RECORD_TREE_EVENT, *event );
        EvolutionaryEvent* tree_event = new(event_mem) EvolutionaryEvent( event_height, epoch_i,
                                                                          event_height, epoch_i,
                                                                          this->current_base(), this->current_base(), 1 );
        tree_event->set_recomb_event_time( event_height );
        tree_event->set_descendants( get_descendants( rec_point.base_node() ) );  // calculate signature for descendants subtended by node
        tree_event->add_leaf_to_tree( &eventTrees[ tree_epoch ] );
    }
}


void ForestState::resample_recombination_position(void) {

    // first, obtain a fresh sequence position for the next recombination, overwriting the existing sample in next_base_
    // (the implementation, in forest.h, then calls sampleNextBase which is overloaded and implemented in this file)
    this->resampleNextBase();

    // then, create private event records, in effect re-doing the work of record_recomb_event
    for (int epoch = 0; epoch < model().change_times_.size(); epoch++) {
        if ((pfparam.record_event_in_epoch[ epoch ] & PfParam::RECORD_RECOMB_EVENT) &&
            epoch <= max_epoch_to_record_) {
            EvolutionaryEvent* old_event = eventTrees[ epoch ];      // pointer to old event to be considered for copying
            EvolutionaryEvent** new_chain = &eventTrees[ epoch ];    // ptr to ptr to current event chain
            // break out the loop if no recombination opportunity has been recorded at this epoch,
            // (implying that no opportunity has been recorded at more ancient epochs either)
            // Note: this assumes that the most recent recombination opportunity is the last opportunity recorded, and first in the chain!
            //       see record_recomb_extension, where care is taken that this assumption is met
            if ( !old_event || !old_event->is_recomb() ||
                 !old_event->recomb_event_overlaps_opportunity_x( current_base() ) ) {
                break;
            }
            // loop over the recombination opportunity records for this epochs
            // there can be several if the tree has nodes in this epoch.
            do {
                // if multiple particle records refer to this event,
                // make a copy of current event before modifying the end_base member
                if (old_event->ref_counter() > 1) {
                    void* event_mem = Arena::allocate( epoch );
                    EvolutionaryEvent* new_event = new(event_mem) EvolutionaryEvent( *old_event );
                    // splice into event tree, and update pointers (true == always treat new_event as a tree)
                    new_event->add_leaf_to_tree( new_chain, true );
                    new_chain = &(new_event->parent_);              // friend function access
                    // update record
                    new_event->set_end_base( this->next_base() );
                } else {
                    old_event->set_end_base( this->next_base() );
                    new_chain = &(old_event->parent_);
                }
                old_event = old_event->parent();

            } while ( old_event && old_event->is_recomb() &&
                      old_event->recomb_event_overlaps_opportunity_x( current_base() ) );
        }
    }
}


void ForestState::includeLookaheadLikelihood( const Segment& segment, const TerminalBranchLengthQuantiles& terminal_branch_lengths ) {

    if (pfparam.auxiliary_particle_filter == 0)
        return;

    // extract rho and mu and nsam
    // note: this assumes the current \rho is constant -- will not work well with hotspots and focused sampling!
    int cur_rec_idx = model().get_position_index();
    double cur_recomb_segment_start = model().change_position(cur_rec_idx);
    if (current_base() <= cur_recomb_segment_start) --cur_rec_idx;
    double recomb_rate = model().recombination_rate( cur_rec_idx );
    if (pfparam.recomb_bias.using_posterior_biased_sampling()) {
        recomb_rate = pfparam.recomb_bias.get_true_rate();
    }
    const double mut_rate = model().mutation_rate();
    const int nsam = model().sample_size();
    double rel_rho[] = {1.0, 0.5};
    double rel_rho_p[] = {0.5, 0.5};
    
    // compute the likelihood.  First singletons
    double likelihood = 1.0;
    double rho_tbl = 2 * recomb_rate * (nsam - 1) / nsam;
    Node* node = nodes()->first();
    int num_leaves = 0;
    vector<const Node*> leaves( nsam );
    vector<double> mut_prob( nsam );
    while (num_leaves < nsam) {
        if (node->in_sample()) {
            num_leaves++;
            int i = node->label() - 1;
            leaves[i] = node;
	}
	node = node->next();
    }
    for (int i=0; i<num_leaves; i++) {
      double p = 0;
      double si = segment.first_singleton_distance[i];
      double li = leaves[i]->getLocalParent()->height();
      if (segment.is_singleton_unphased[i]) {
	// we don't know which branch contained the singleton, so add length of alternative branch
	li += leaves[i+1]->getLocalParent()->height();
      }
      double rel_mut_rate = segment.relative_mutation_rate[i];
      double li_mu = li * mut_rate * rel_mut_rate;
      mut_prob[i] = li_mu;
      if (segment.is_singleton_unphased[i]) {
	// NOTE: adding mutation rates is not strictly correct for usage in the cherry model below, but
	// as the terms is there to provide a soft landing in an edge case, this will not affect overall behaviour.
	mut_prob[i+1] = li_mu;
      }
      for (int r = 0; r < sizeof(rel_rho)/sizeof(*rel_rho); r++) {
	// integrate over rate of change: expected, and half expected, to
	// model autocorrelation of branch lengths across recombination events
	double li_rho = li * rho_tbl * rel_rho[r];
	double fe = fastexp_approx(-(li_rho + li_mu) * abs(si));
	for (int q = 0; q < terminal_branch_lengths.quantiles.size(); ++q) {
	  // integrate over the distribution of terminal branch lengths, for this branch i
	  double qbot = (q == 0 ? 0.0 : terminal_branch_lengths.quantiles[q-1]);
	  double qtop = (q == terminal_branch_lengths.quantiles.size()-1 ? 1.0 : terminal_branch_lengths.quantiles[q]);
	  double l_prime = terminal_branch_lengths.lengths[ i ][ q ];
	  double lprime_mu = l_prime * mut_rate * rel_mut_rate;
	  double div = (li_rho + li_mu - lprime_mu);
	  if (abs(div) < (li_rho + li_mu + lprime_mu) * 1e-5) {
	    lprime_mu = lprime_mu * 1.0001;
	  }
	  if (si > 0) {
	    // mutation present
	    p += rel_rho_p[r] * (qtop-qbot) * (( li_rho * lprime_mu * fastexp_approx(-lprime_mu * si) +
						 (li_mu - lprime_mu) * (li_rho + li_mu) * fe )
					       / (li_rho + li_mu - lprime_mu) );
	  } else {
	    // missing data - no mutation
	    p += rel_rho_p[r] * (qtop-qbot) * (( li_rho * fastexp_approx(-lprime_mu * (-si)) +
						 (li_mu - lprime_mu) * fe )
					       / (li_rho + li_mu - lprime_mu) );
	  }
	}
      }
      //cout << " Leaf " << i << " d=" << si << " ht=" << li << " p=" << p;
      likelihood *= p;
      // do not double-count unphased singletons
      if (segment.is_singleton_unphased[i]) {
	i++;
      }
    }

    // next doubletons
    if (pfparam.auxiliary_particle_filter >= 2) {
        double l_mean = 0.0;
        for (int i=0; i<nsam; i++) l_mean += terminal_branch_lengths.lengths[ i ].back() / nsam;
        double rho_c = 4 * recomb_rate * (nsam - 2) / nsam;
        double rhoprime_c = recomb_rate * (nsam - 1);
        double p_equilibrium = 2.0/(3*(nsam-1));
        int ph1, ph2;
        for (const Segment::Doubleton& d : segment.doubleton) {
            // do we have the cherry?  Check greedily by considering all phasings
            for (ph1 = 0; ph1 <= d.unphased_1; ph1++) {
                for (ph2 = 0; ph2 <= d.unphased_2; ph2++) {
                    if (leaves[d.seq_idx_1 + ph1]->parent() == leaves[d.seq_idx_2 + ph2]->parent()) {
                        // we do
                        double l = leaves[d.seq_idx_1 + ph1]->getLocalParent()->height();
                        double p = 0;
                        for (int r = 0; r < sizeof(rel_rho)/sizeof(*rel_rho); r++) {
                            // integrate over rate of change: expected, and half expected, to
                            // model autocorrelation of branch lengths across recombination events
                            double exp_rho = fastexp_approx(-rho_c * rel_rho[r] * l * d.last_evidence_distance);
                            // note - do not include a factor for a mutation; we don't want to model the length of the cherry root branch
                            p += rel_rho_p[r] * exp_rho + p_equilibrium * (1.0 - exp_rho);
                        }
                        likelihood *= p;
                        //cout << " Ch " << d.seq_idx_1+ph1 << d.seq_idx_2+ph2 << " p=" << p;
                        ph1 = ph2 = 99;
                    }
                }
            }
            if (ph1 < 99) {
                // we don't have a cherry.
                // include a term for a double mutation (i.e. a single one, since one mutation isn't included normally), to cover the case that first_evidence_distance == 0
                // the mutation rate is taken to be average of the mutation rates in the two terminal branches
                double mutprob = (mut_prob[d.seq_idx_1] + mut_prob[d.seq_idx_2]) * 0.5;
                double p = 0;
                for (int r = 0; r < sizeof(rel_rho)/sizeof(*rel_rho); r++) {
                    // integrate over rate of change: expected, and half expected, to
                    // model autocorrelation of branch lengths across recombination events
                    p += rel_rho_p[r] * (mutprob + (1.0 - mutprob) * p_equilibrium * (1.0 - fastexp_approx(-rhoprime_c * rel_rho[r] * l_mean * d.first_evidence_distance)));
                }
                likelihood *= p;
                //cout << " NoCh " << d.seq_idx_1+ph1 << d.seq_idx_2+ph2 << " fed=" << d.first_evidence_distance << " p=" << p;
            }
        }
    }
        
    // finally, splits
    if (segment.first_split_distance > -1 && pfparam.auxiliary_particle_filter >= 3) {
        // probability of no change to the current topology.
        // Assume 1/2 of recombinations cause change in split
        double rate_of_change = getLocalTreeLength() * recomb_rate / 2;
        double p_nochange = fastexp_approx( -rate_of_change * segment.first_split_distance );

        // probability of split data given current tree.
        bool ancestral_aware = false;
        double p_splitdata = calculate_likelihood( ancestral_aware, segment.allelic_state_at_first_split );

        // probability of correct split.  This would be n-choose-k at equilibrium.  But the current tree, if it is not
        // supporting the split, is likely to be -almost- supporting the split; and in addition we need to have a heavy-tailed
        // distribution.  So instead, assume that current tree has one lineage at the wrong side of the split.  We need a
        // recombination in the correct branch (1/n) and need a coalescence into the split (k/n), together probability k/n^2,
        // and recombiations and coalescences occuring quickly -- factor 1/4. (justify??)
        // -- consider n-choose-k ?
        int k = segment.mutation_count_at_first_split;
        double p_correct_split = k / double(4.0 * nsam * nsam);
        // for testing: let apf=2 do split distance, apf=4 uses n-choose-k factor
        if (pfparam.auxiliary_particle_filter == 4)
            p_correct_split = 1.0 / nchoosek( nsam, k );

        // length of the split branch. Under constant-pop-size coalescent model (CCM) total length is l(n) = 2(1 + ... + 1/(n-1)).
        // The split branch has k descendants; since hanging configurations are uniformly distributed (Xavier Didelot) other
        // branches at that time have the same expected number of descendants, so the number of branches is in expectation
        // n/k.  The length of the next branch to coalesce into this upper subtree is l( n/k + 1 ) - l( n/k ) = 2k / n, which
        // is a measure of the length of the subtree supporting the split.  After all branches have coalesced into the tree,
        // the branch will have undergone further coalescences, so we divide this length by 2 (justification?)
        // Now, our model is not CCM, so scale by the ratio of total branch length under CCM, 2(gamma + ln(n)), to the
        // empirical total branch length ETBL, to get for the  approximate length of the split branch,
        //  k * ETBL / ( 2 * n * ( gamma + ln(n) ) ).
        double etbl = terminal_branch_lengths.mean_total_branch_length;
        double split_branch_length = k * etbl / ( 2 * nsam * ( 0.577 * log(nsam) ) );
        double p = p_nochange * p_splitdata + (1.0 - p_nochange) * p_correct_split * mut_rate * split_branch_length;
        likelihood *= p;
        //double mut_branch_length = p_splitdata / (this->model().mutation_rate()); // only for display
        //cout << " Split L=" << p_splitdata << " true brlen=" << mut_branch_length << " dist=" << segment.first_split_distance
        //     << " p_noch=" << p_nochange << " p_avg=" << p_correct_split * mut_rate * split_branch_length << " p=" << p << endl;
    }
        
    //cout << "Lookahead: " << likelihood << endl;
    
    // actually include the lookahead likelihood into the pilot weight
    lookahead_weight_ *= likelihood;
    pilot_weight_ *= likelihood;
}


/*!
 * Calculate the marginal likelihood of a node recursively.
 *
 * * */

inline void ForestState::cal_partial_likelihood_infinite(Node * node, const vector<int>& haplotype, double marginal_likelihood[2]) {

    // deal with the case that this node is a leaf node
    if ( node->first_child() == NULL ) {
        assert ( node->in_sample() );                                // we only traverse the local tree, therefore leaf node must be in sample
        assert ( node->label() <= haplotype.size() );
        assert ( node->second_child() == NULL );                     // leaf (if a node has no first child, it won't have a second child)
        int mutation_state = haplotype[ node->label() - 1 ];
        marginal_likelihood[0] = mutation_state == 1 ? 0.0 : 1.0;    // also encode state==-1 (missing data) as 1.0
        marginal_likelihood[1] = mutation_state == 0 ? 0.0 : 1.0;    // also encode state==-1 (missing data) as 1.0
      return;
    }

    // Genealogy branch lengths are in number of generations, the mutation rate is unit of per site per generation
    double mutation_rate = this->model().mutation_rate();
    Node *left = trackLocalNode(node->first_child()); 
    Node *right = trackLocalNode(node->second_child());
    double t_left = node->height() - left->height();
    double t_right = node->height() - right->height();
    assert (t_left >= 0 && t_right >= 0 && mutation_rate > 0);
    double p_left =fastexp(-t_left * mutation_rate);   // probability that no mutation occurred along branch (infinite sites model)
    double p_right = fastexp(-t_right * mutation_rate); // probability that no mutation occurred along branch (infinite sites model)
    double marg_like_left[2];
    double marg_like_right[2];
    cal_partial_likelihood_infinite(left, haplotype, marg_like_left);
    cal_partial_likelihood_infinite(right, haplotype, marg_like_right);

    marginal_likelihood[0] = ( marg_like_left[0]*p_left + marg_like_left[1]*(1-p_left) )
        * ( marg_like_right[0]*p_right + marg_like_right[1]*(1-p_right) );
    marginal_likelihood[1] = ( marg_like_left[1]*p_left + marg_like_left[0]*(1-p_left) )
        * ( marg_like_right[1]*p_right + marg_like_right[0]*(1-p_right) );

    return;
}



/*!
 * \brief Calculate the marginal likelihood of the genealogy at data site i,
 *  If there is no data given at the site i, return likelihood as 1.
 * @ingroup group_pf_resample
 */
double ForestState::calculate_likelihood( bool ancestral_aware, const vector<int>& haplotype ) {
    double marginalLikelihood[2];
    cal_partial_likelihood_infinite(this->local_root(), haplotype, marginalLikelihood);
    double prior[2];
    if ( ancestral_aware ) {
        prior[0] = 1;
        prior[1] = 0;
    } else {
        prior[0] = 0.5;
        prior[1] = 0.5;
    }
    double likelihood = marginalLikelihood[0]*prior[0] + marginalLikelihood[1]*prior[1];
    return likelihood;
}


Node* ForestState::trackLocalNode(Node *node) const {
    assert( node->local() );
    if (node->countChildren() == 0) return node;
    if (node->countChildren() == 1) return trackLocalNode(node->first_child());
    
    assert( node->countChildren() == 2 );
    assert( node->first_child()->local() || node->second_child()->local() );
    
    if ( node->first_child()->local() ) {
        if (node->second_child()->local()) return node;
        else return trackLocalNode(node->first_child());
    }
    else return trackLocalNode(node->second_child());
}


double ForestState::trackLocalTreeBranchLength( const vector<int>& data_at_site ) {
    total_local_branch_length_ = 0.0;
    trackSubtreeBranchLength( local_root(), data_at_site );
    return total_local_branch_length_;
}


double ForestState::trackSubtreeBranchLength( Node* currentNode, const vector<int>& data_at_site ) {
    // returns total length of branches that carry data below currentNode, which is local;
    // or -1 if currentNode has no data-carrying descendants.
    // Side effect: sets total_local_branch_length_ if currentNode is a potential data root (both children carry data)

    if (currentNode->in_sample() ) {                             // current node is a leaf node
        if (data_at_site[currentNode->label()-1] >= 0) return 0; // leaf node carries data
        else return -1;                                          // leaf node carries no data
    }
    Node* child_left = trackLocalNode(currentNode->first_child());
    Node* child_right = trackLocalNode(currentNode->second_child());    
    double stbl_left = trackSubtreeBranchLength( child_left, data_at_site );
    double stbl_right = trackSubtreeBranchLength( child_right, data_at_site );
    if (stbl_left >= 0.0)  stbl_left  += currentNode->height() - child_left->height();
    if (stbl_right >= 0.0) stbl_right += currentNode->height() - child_right->height();
    // potential data root?
    if (stbl_left >= 0.0 && stbl_right >= 0.0) {
        total_local_branch_length_ = stbl_left + stbl_right;
        return total_local_branch_length_;
    }
    if (stbl_left >= 0.0)
        return stbl_left;
    else
        return stbl_right;
}


double ForestState::find_delay( double coal_height ) {
    size_t indx = 0;
    while (indx+1 < model().change_times().size() &&
           model().change_times()[indx+1] <= coal_height )
        indx++;
    double delay = model().application_delays[indx];
    return delay;
}


void ForestState::extend_ARG ( double extend_to, int leaf_status, const vector<int>& data_at_site,
                               vector<ForestState*>* pParticleContainer) {

    double updated_to = this->site_where_weight_was_updated();
    double track_local_tree_branch_length;
    switch (leaf_status) {
    case -1:
        track_local_tree_branch_length = 0;
        break;  // all data missing; no observable mutations
    case 1:
        track_local_tree_branch_length = getLocalTreeLength();
        break;  // no data missing; all observable
    default:
        track_local_tree_branch_length = trackLocalTreeBranchLength( data_at_site );
        break;  // mixed case -- need to compute
    }

    assert (updated_to >= this->current_base());
    while ( updated_to < extend_to ) {

        // First, update the likelihood up to either extend_to or the end of this state
        double new_updated_to = min( extend_to, this->next_base() );

        // Calculate the total tree length of the subtree over leaf nodes that carry data.
        // For leaf nodes that carry no data, there is no evidence for
        // presence or absence of mutations on corresponding branches.
        // These branches do not contribute to localTreeBranchLength
        adjustWeights( fastexp( -model().mutation_rate()
                                * track_local_tree_branch_length
                                * (new_updated_to - updated_to) ) );

        assert (abs(track_local_tree_branch_length - trackLocalTreeBranchLength( data_at_site )) < 1e-4);
        
        // Update importance sampling correction for the weight of the particle.  This requires
        // some explanation.
        //
        // Suppose we're looking at a segment  [0,T)  without recombinations, followed by a
        // recombination at a point  y  in a tree Y.  Let  mu(Y)  be the tree's total branch
        // length, and  rho(y) = rho  the (constant) recombination rate per generation per
        // nucleotide at the recombination point  y.  The probability density of this event is
        // then
        //
        //      exp( -\int_{t=0}^T \int_{y in Y}  rho(y) dt dy ) * rho(y) dt dy
        //   =  exp( -rho T mu(Y) ) * rho dt dy
        //
        // Now suppose we sample from a different process where  rho'(y)  is not constant across
        // the tree  Y  but depends on  y.  The imporance weight for bringing a sample from that
        // density to the desired probability density above, is the ratio of the two densities:
        //
        //   exp( -int_{t=0}^T \int_{y in Y}  (rho(y) - rho'(y)) dt dy ) * (rho(y) / rho'(y))
        //
        // We implemented  rho'(y)  as the normal recombination process with rate  rho  but
        // occurring on a "weighted" tree, so that the weighted tree length is  mu'(Y) =
        //  \int_{y in Y} w(y) dy,  where the weight  w(y)  depends only on the height in the
        // tree.  Another view is that the new recombination rate  rho'(y) = w(y) rho.  Using
        // this notation, the importance weight becomes
        //
        //   IW = exp( -T (mu(Y) - mu'(Y)) / w(y)
        //
        // The first factor is implemented by extend_ARG, and applied immediately to the particle.
        // This factor tends not to depend strongly on the location of  y.  However, we bias
        // towards recent recombinations, and for those  w(y)  tends to be large, so  IW  tends
        // to be small.  In order not to immediately "undo" (through sampling) the effect of
        // biasing towards recent recombinations, the factors  1/w(y)  are applied with a delay.
        // These factors are calculated by the function  SampleNextGenealogy,  which amonst others
        // samples a point  y  on the current tree.
        //
        // Note that the importance weight is for a single particle, even if multiplicity() > 1.
        if(model().biased_sampling) {
          adjustWeights( importance_weight_over_segment( updated_to, new_updated_to ) );
        }

        // Rescue the invariant
        updated_to = new_updated_to;
        setSiteWhereWeightWasUpdated( updated_to );

        // Next, if we haven't reached extend_to now, sample new genealogy, and a new recomb point
        if ( updated_to < extend_to ) {

	    if (updated_to == model().getCurrentSequencePosition()) {

	        this->sampleNextGenealogy( true );    // only advances recombination event pointer in this case
	        this->sampleNextBase( true );         // obtain nest recomb rate change event, and push it onto vector
                record_recomb_extension();

	    } else {

		// a recombination has occurred
                int mult = multiplicity();
		first_event_height_ = -1.0;
		first_coal_height_ = -1.0;

                // We shall implement the change.  Spawn a new particle if necessary
                if (mult > 1) {
		      
                    // A recombination has occurred on a particle with multiplicity > 1.
                    // Spawn a new particle with multiplicity one less, resample its recombination
                    //  position, and deal with it later in the loop in particleContainer.  At the
                    //  same time, make the current particle have multiplicity one, and deal with it
                    //  the normal way.
                    assert( pParticleContainer != NULL );
                    pParticleContainer->push_back( spawn() );
                    pParticleContainer->back()->set_current_base( updated_to );
                    pParticleContainer->back()->setMultiplicity( mult - 1 );
                    // save recombination index, owned by singleton Model, into *this Particle, and restore from back() particle;
                    // then sample new recombination locus for the new particle; and recover recombination index
                    save_recomb_state();
                    pParticleContainer->back()->restore_recomb_state();
                    pParticleContainer->back()->resample_recombination_position();
                    pParticleContainer->back()->save_recomb_state();
                    restore_recomb_state();
                }
		    
                // set new multiplicity
                setMultiplicity( 1 );

                // sample new tree
                recombination_bias_importance_weight_ = 1.0;
                double importance_weight = this->sampleNextGenealogy( true );

                if (leaf_status == 0) track_local_tree_branch_length = trackLocalTreeBranchLength( data_at_site );
                if (leaf_status == 1) track_local_tree_branch_length = getLocalTreeLength();
		
                // Implement importance weight, with appropriate delay
                if (first_event_height_ == -1.0) throw std::runtime_error("No coalescence found where one was expected");
		    
                // If the coalescence occurred above top bias height, i.e. the sampling failed
                // to produce an early coalescence as intended, then apply the importance weight
                // factor due to the recombination bias immediately, so that we don't pollute
                // the particles.  However, always apply the importance weight due to guiding
                // with the appropriate delay
                double delay_height = first_event_height_;
                if (pfparam.delay_type == PfParam::RESAMPLE_DELAY_COAL)   delay_height = first_coal_height_;
                if (pfparam.delay_type == PfParam::RESAMPLE_DELAY_RECOMB) delay_height = rec_point.height();
                int idx = 0;
                while (idx+1 < model().bias_heights().size() && model().bias_heights()[idx+1] < delay_height) ++idx;
                if (model().bias_strengths()[idx] == 1.0) {
                    // the height of the focal event (recombination, coalescence, or migration, as per delay_type)
                    // is not biased, so apply the importance weight for the recombination bias immediately, and factor
                    // it out of the overall importance weight (which may include a factor for recombination guiding)
                    adjustWeights( recombination_bias_importance_weight_ );
                    importance_weight /= recombination_bias_importance_weight_;
                }
                double delay = find_delay( delay_height );
                // enter the importance weight, and apply it semi-continuously:
                // in 3 equal factors, at geometric intervals (double each time)
                // (See implementation in particle.hpp: void applyDelayedAdjustment(), and
                //  class DelayedFactor)
                adjustWeightsWithDelay( importance_weight, delay, 3 );
		    
                // sample a new recombination position (this calls the virtual overloaded
                // sampleNextBase())   Note: it returns an importance "rate", which is ignored;
                // see comments in sampleNextBase below.
                // The new recombination position does take account of multiplicity(), so that
                // particles with high multiplicity() have higher effective recombination rates.
                this->sampleNextBase( true );
                this->record_recomb_extension();

                // If this is an extension after an actual recombination, record the recomb event
                record_recomb_event( rec_point.height() );

	    }
        }
        // record current position as recombination, so that resampling of recomb. position starts from here
        set_current_base( updated_to );
    }

    // apply weights that we passed during the extension above
    while ( !delayed_adjustments.empty()
            && delayed_adjustments.top().application_position < extend_to ) {
        
        this->applyDelayedAdjustment();
        
    }

}


std::string ForestState::newick(Node *node) {
  std::stringstream tree;
  if (node->in_sample()) tree << node->label();
  else {
    Node *left = node->getLocalChild1();
    Node *right = node->getLocalChild2();

    tree << "(" << this->newick(left) << ":" <<
           (node->height() - left->height()) * this->model().scaling_factor() <<
           "," << this->newick(right) << ":" <<
           (node->height() - right->height()) * this->model().scaling_factor() << ")";
  }

  return tree.str();
}


//
// biased sampling
//

double ForestState::getWeightedLocalTreeLength( const bool recomb_bias,
                                                const bool recomb_guide ) {
    double total_length = 0.0;
    double dummy = 0.0;
    double rate1 = 0.0, rate2 = 0.0;
    bool branch1 = local_root()->first_child() && local_root()->first_child()->samples_below() > 0;
    bool branch2 = local_root()->second_child() && local_root()->second_child()->samples_below() > 0;

    if ( branch1 )
        rate1 = sampleOrMeasureWeightedTree( local_root()->first_child(), total_length, dummy, recomb_bias, recomb_guide );
    if ( branch2 )
        rate2 = sampleOrMeasureWeightedTree( local_root()->second_child(), total_length, dummy, recomb_bias, recomb_guide );
    // most of the recombinations in the branches from the root node cannot be assigned with confidence to either.
    // therefore make sure that the rate is equal in both.  This avoids the issue that for a cluster of deep
    // coalescences, each recombination gets biased with an importance weight of approx. 0.5, potentially leading
    // to a very large importance weight, which when postponed for a while would lead to overflows.
    double rate = max(rate1, rate2);
    if ( branch1 ) accumulateBranchLengths( rate, recomb_bias, local_root()->first_child(), total_length, dummy );
    if ( branch2 ) accumulateBranchLengths( rate, recomb_bias, local_root()->second_child(), total_length, dummy );
    return total_length;
}


// Function to calculate the total length of a tree (if cumul_length == 0), or sample
// a point on it (if length_left < 0).  The tree can be weighted for recombinations
// at particular heights (recomb_bias == true) and/or weighted to reflect the
// recombination guide (recomb_guide == true), or neither.
//
// \return the recombination rate above the node
//
double ForestState::sampleOrMeasureWeightedTree( const Node* node,
                                                 double& cumul_length,
                                                 double& local_weight,
                                                 const bool recomb_bias,
                                                 const bool recomb_guide ) {
    double rate1 = 0, rate2 = 0;
    int branches_below = 0;
    //
    // recursively enter the child nodes below, measure or sample, and obtain guide rates on descending branches
    //
    if ( node->first_child() && node->first_child()->samples_below() > 0 ) {
        rate1 = sampleOrMeasureWeightedTree( node->first_child(), cumul_length, local_weight, recomb_bias, recomb_guide );
        accumulateBranchLengths( rate1, recomb_bias, node->first_child(), cumul_length, local_weight );
        branches_below += 1;
    }
    if ( node->second_child() && node->second_child()->samples_below() > 0 ) {
        rate2 = sampleOrMeasureWeightedTree( node->second_child(), cumul_length, local_weight, recomb_bias, recomb_guide );
        accumulateBranchLengths( rate2, recomb_bias, node->second_child(), cumul_length, local_weight );
        branches_below += 1;
    }
    //
    // calculate guide rate above node and return
    //
    if (branches_below == 0) {                // leaf node
        if (recomb_guide) {
            const RecombBiasSegment& rbs = pfparam.recomb_bias.get_recomb_bias_segment(model().get_position_index());
            if (current_base() < rbs.get_locus() || current_base() > rbs.get_end()) {
                cout << "Problem: current base: " << current_base() << endl;
                cout << " rbs segment: [" << rbs.get_locus() << "," << rbs.get_end() << ")" << endl;
                throw std::runtime_error("Current base not in expected rbs segment!  Should never happen!");
            }
            return rbs.get_leaf_rate( node->label()-1 );
        } else {
            return 1.0;  // relative rate - uniform across tree
        }
    } else if (branches_below == 1) {         // single child branch -- transmit rate unchanged
        return rate1 + rate2;
    } else {
        // two child branches.  Calculate (unbiased, unnormalized) rate on branch doing down from node,
        // as the arithmetic average of the rates on the two incoming branches.  Previously we
        // tried the hyperbolic average, but this causes recombinations to be focused onto
        // recent branches when the tree doesn't fit the data well.  These carry large importance
        // weights, which are delayed for a long time, leading to overflows.
        return (rate1 + rate2) * 0.5;
    }
}


void inline ForestState::accumulateBranchLengths( double rate_above, 
                                                  bool recomb_bias, 
                                                  const Node* node, 
                                                  double& cumul_length, 
                                                  double& local_weight ) {
    size_t time_idx = 0;
    double lower_end = 0.0;            // this is correct irrespective of whether recomb_bias is used or not
    double upper_end = 1e100;          // default for !recomb_bias
    double local_weight_ = rate_above; // default for !recomb_bias
    do {
        if (recomb_bias) {
            local_weight_ = rate_above * model().bias_strengths()[time_idx];
            upper_end = model().bias_heights()[time_idx+1];
        }
        lower_end = max( lower_end, node->height() );
        upper_end = min( upper_end, node->parent_height() );
        double weighted_branch_length = local_weight_ * max(0.0, upper_end - lower_end);
        if (cumul_length < 0 && cumul_length + weighted_branch_length >= 0) {
            // a sample fell on the branch.  Note, 0 < initial cumul_length <= total branch length
            local_weight = local_weight_;                                              // store local weight
            double sampled_height = lower_end + (-cumul_length) / local_weight_;       // compute height
            rec_point = TreePoint(const_cast<Node*>(node), sampled_height, false);     // store node and height
            if (model().bias_heights().size() >= 2 && sampled_height < model().bias_heights()[1]) {
                //recent_recombination_count++;
            }
        }
        cumul_length += weighted_branch_length;
        lower_end = upper_end;
        ++time_idx;
    } while ( upper_end < node->parent_height() );
}


//
//
// The actual functions doing the sampling, and calculating of importance weights, are below:
//
//


double ForestState::samplePoint( bool record_and_bias ) {

    double dummy_rate = 0.0;   _unused(dummy_rate);
    bool use_recomb_guide = record_and_bias && pfparam.recomb_bias.using_posterior_biased_sampling();
    bool use_recomb_bias  = record_and_bias && model().biased_sampling;
    double full_weighted_local_tree_length = getWeightedLocalTreeLength( use_recomb_bias, use_recomb_guide );
    double sample_length = (random_generator()->sample() - 1.0) * full_weighted_local_tree_length;  // guaranteed < 0
    double recomb_weighted_local_tree_length;   // tree length with only recombination weighting
    if (use_recomb_bias) {
        if (use_recomb_guide) {
            recomb_weighted_local_tree_length = getWeightedLocalTreeLength( use_recomb_bias, false );
        } else {
            recomb_weighted_local_tree_length = full_weighted_local_tree_length;
        }
    } else {
        recomb_weighted_local_tree_length = getLocalTreeLength();
    }
    // sample; updates full_local_weight and put sample in rec_point
    double full_local_weight = 0.0;
    double rate1 = 0.0, rate2 = 0.0;
    bool branch1 = local_root()->first_child() && local_root()->first_child()->samples_below() > 0;
    bool branch2 = local_root()->second_child() && local_root()->second_child()->samples_below() > 0;

    if ( branch1 )
        rate1 = sampleOrMeasureWeightedTree( local_root()->first_child(), sample_length, full_local_weight, use_recomb_bias, use_recomb_guide );
    if ( branch2 )
        rate2 = sampleOrMeasureWeightedTree( local_root()->second_child(), sample_length, full_local_weight, use_recomb_bias, use_recomb_guide );
    // most of the recombinations in the branches from the root node cannot be assigned with confidence to either.
    // therefore make sure that the rate is equal in both.  This avoids the issue that for a cluster of deep
    // coalescences, each recombination gets biased with an importance weight of approx. 0.5, potentially leading
    // to a very large importance weight, which when postponed for a while would lead to overflows.
    double rate = max(rate1, rate2);

    if ( branch1 ) accumulateBranchLengths( rate, use_recomb_bias, local_root()->first_child(), sample_length, full_local_weight );
    if ( branch2 ) accumulateBranchLengths( rate, use_recomb_bias, local_root()->second_child(), sample_length, full_local_weight );

    if ( sample_length < 0 ) {
        throw std::runtime_error("Failure to sample point -- should never happen!");
    }

    // compute importance weight.  Note -- this returns the ratio of densities sampling the
    // tree point  y  conditional on a recombination sequence position  x, as the current
    // implementation samples recombination positions according to the guide recombination rate,
    // but disregards topological guiding and recombination biasing (and uses importance weights
    // to correctfor the guide recombination rate only.)  It is the responsibility of this
    // code to compute the importance weight to account for biasing on the tree.
    double sampled_density = full_local_weight / full_weighted_local_tree_length;
    double target_density  = 1.0               / getLocalTreeLength();
    double importance_weight = target_density / sampled_density;

    // compute importance weight due to recombination biasing only (as opposed to full biasing, which
    // also includes recombination guiding), to allow the delay for this importance weight to be
    // applied independently to that for guiding.  First calculate the recombination bias weight.
    double recomb_bias_wt = 1.0;
    if (use_recomb_bias) {
        int idx=0;
        while (model().bias_heights()[idx+1] < rec_point.height())
            ++idx;
        recomb_bias_wt = model().bias_strengths()[idx];
    }
    double recombination_density = recomb_bias_wt / recomb_weighted_local_tree_length;
    recombination_bias_importance_weight_ = target_density / recombination_density;
    
    // done -- importance weight (and delay) are implemented in extend_ARG,
    // via sampleNextGenealogy (scrm/forest.cc)
    return importance_weight;
}




/**
 * Function for adjusting the importance weight predata across a segment
 *
 * This function must reflect the actual sampling algorithm in sampleNextBase
 * In addition, it accounts for the difference between the posterior-weighted recombination rate
 *   and the true recombination rate
 */
double ForestState::importance_weight_over_segment( double previously_updated_to,
                                                    double update_to) {

    int cur_rec_idx = model().get_position_index();
    double cur_recomb_segment_start = model().change_position(cur_rec_idx);
    // the update segment may, or may not, overlap the "current" recombination rate segment,
    // since that is moved along as soon as the simulated segment hits a change point.
    if (update_to <= cur_recomb_segment_start) {
        // The "current" recombination rate segment does not overlap; move the index one back
        assert (cur_rec_idx > 0);
        --cur_rec_idx;
    }
    double sequence_distance = update_to - previously_updated_to;
    double posterior_weighted_recombination_rate = model().recombination_rate( cur_rec_idx );
    double true_recomb_rate = posterior_weighted_recombination_rate;

    if (pfparam.recomb_bias.using_posterior_biased_sampling()) {
        true_recomb_rate = pfparam.recomb_bias.get_true_rate();
        const RecombBiasSegment* rbs = &pfparam.recomb_bias.get_recomb_bias_segment( cur_rec_idx );
        assert(update_to > rbs->get_locus());
        assert( rbs->get_rate() == posterior_weighted_recombination_rate );
        if (rbs->get_rate() != posterior_weighted_recombination_rate ) {
            cout << "Problem: "
                 << rbs->get_rate()<< " != " << posterior_weighted_recombination_rate
                 << " at idx " << cur_rec_idx << endl;
            throw std::runtime_error("Guide and model recombination rate differ.  This should never happen!");
        }
    }

    double target_rate =  (sequence_distance
                           * true_recomb_rate
                           * getLocalTreeLength());
    double sampled_rate = (sequence_distance
                           * posterior_weighted_recombination_rate
                           * getLocalTreeLength() );

    // return target_density / sampled_density, that is, exp( -target_rate ) / exp( -sampled_rate )
    // note that, similar to the posterior weight, the importance weight is independent of the particle's
    // multiplicity; the importance weight is for a single underlying particle, while the multiplicity
    // simulates a virtual ensemble of such particles.
    double importance_weight = fastexp( sampled_rate - target_rate );

    return importance_weight;
}


/**
 * Function for sampling the sequence position of the next recombination event
 * under a biased sampling procedure.
 *
 * \return the importance 'rate' needed to compute importance weight  exp(-rate L)  correcting
 *         for biased sampling.
 *
 * Note: the function samplePoint()
 *
 */

double ForestState::sampleNextBase( bool record_and_bias ) {

    // three sources of biasing:
    // "focus": varying rates across epochs, used to increase recombinations very low (and very high?) in tree
    // "tree guide": varying rates across the tree, to guide recombinations towards branches with posterior support
    // "position guide": varying rates across the sequence, to guide recombinations towards branches with posterior support
    
    double distance_until_rate_change  = model().getNextSequencePosition() - current_base();
    // we use conditional weighting (in samplePoint) to implement focus and tree guiding, so that while locally within the
    // tree rates change, the overall rate remains that without focus and tree guiding.  So the "weighted" local tree length
    // is the same as the normal local tree length
    double wtd_local_tree_len          = getLocalTreeLength();
    // position guiding is implemented by making the recombination rate position-dependent within model()
    double recomb_rate                 = model().recombination_rate();    // posterior weighted

    if (!record_and_bias) {
        // for the purpose of calibrating the lag times, do not use any biasing.  So, use the true recombination rate
        // rather than the posterior weighted one.  (Importance_weight_over_segment above is not used when calibrating,
        //  so it doesn't need to know the record_and_bias setting.)
        if (pfparam.recomb_bias.using_posterior_biased_sampling()) {
            recomb_rate                = pfparam.recomb_bias.get_true_rate();
            distance_until_rate_change = 1e100;
        }
    }
    // total recombination rate per nt, with all three biasing sources factored in
    double weighted_pernuc_recomb_rate = recomb_rate * wtd_local_tree_len;
    // total recombination rate per nt, with additionally the particle's multiplicity factored in
    double weighted_pernuc_particle_recomb_rate = weighted_pernuc_recomb_rate * multiplicity();
    // sample a next event position
    double length = random_generator()->sampleExpoLimit(weighted_pernuc_particle_recomb_rate,
                                                        distance_until_rate_change);

    if (length == -1) {          // No recombination until the model changes

        set_next_base(model().getNextSequencePosition());
        if (next_base() < model().loci_length()) {
            writable_model()->increaseSequencePosition();
        }

    } else {                     // A recombination occurred in the sequence segment

        double current_position = current_base();
        double next_position = current_base() + length;
        if (current_position == next_position) {
            cerr << "Warning: Next recombination position is identical due to underflow at "
                 << "pos=" << current_position << "; length=" << length << endl;
            cerr << "Check that recombination rates are sensible! (at this position, rate = " << recomb_rate << " per nt;"
                 << " tree length " << wtd_local_tree_len << " generations; multiplicity = " << multiplicity() << ")" << endl;
            next_position = nextafter( current_base(), current_base()*2 );
        }
        set_next_base(next_position);

    }
    assert(next_base() > current_base());
    assert(next_base() <= model().loci_length());

    // the importance weight for this event is calculated by importance_weight_over_segment
    // the multiplicity is dealt with within extend_ARG
    return 1.0;
}



/**
 * Function to modify the tree after we encountered a recombination on the
 * sequence. Also samples a place for this recombination on the tree, marks the
 * branch above as non-local (and updates invariants) if needed, cuts the
 * subtree below away and starts a coalescence from it's root.
 * @ingroup group_scrm_next
 * @ingroup group_pf_update
 */
double ForestState::sampleNextGenealogyWithoutImplementing() {

    // advance current_rec_ counter, so that current_base() is up to date, which will ensure
    // that TimeIntervalIterator does all the necessary pruning, which ensures that this
    // routine sees the exact same tree as sampleNextGenealogy
    current_rec_++;
  
    // Prune the tree (test if this is necessary; it may not be for coalescences since
    //  the tree may be pruned before branches are sampled from)
    // Also set the primary root, which may get un-set due to pruning.
    for (TimeIntervalIterator ti(this, this->nodes_.at(0)); ti.good(); ++ti) {
      if ((*ti).time_interval_iterator()->end_node() != NULL &&
	  (*ti).time_interval_iterator()->end_node()->is_root())
	set_primary_root( (*ti).time_interval_iterator()->end_node() );
    }

    if (current_base() == model().getCurrentSequencePosition()) {
        current_rec_--;   // undo the advance above
        return -1;        // recombination rate change; signal that no recombination has occurred
    }

    contemporaries_.clear(true);

    // Sample the recombination point into TreePoint rec_point member
    double importance_weight = samplePoint( true );

    // Virtual (local) root node that will be coalescing into the existing tree
    Node start_node( rec_point.base_node()->height() );     // is_root() and local() by default
    start_node.set_population( rec_point.base_node()->population() );
    start_node.set_first_child( rec_point.base_node() );    // store actual node is first_child()

    // We can have one or active local nodes: If the coalescing node passes the
    // local root, it also starts a coalescence.
    Node second_node( local_root()->height() );
    second_node.set_population(  local_root()->population() );
    second_node.set_parent(      local_root()->is_root() ? NULL : local_root()->parent() );      // sets is_root()
    second_node.make_nonlocal(   local_root()->last_update() ); // sets local()
    second_node.set_first_child( local_root() );                // store actual node in first_child()

    // Initialize Temporary Variables
    // tmp_event_.time() is the time of the current event, and may not be
    // part of the tree (the "new" events will not be, as the tree won't be
    // modified).  The time intervals processed in the main loop below will
    // start at tmp_event_time().
    tmp_event_ = Event();
    coalescence_finished_ = false;

    if (start_node.height() > second_node.height()) {
        tmp_event_.set_time( second_node.height() );
        set_active_node(0, &second_node);
        set_active_node(1, &start_node);
    } else {
        tmp_event_.set_time( rec_point.height() );   // recombination height, not height of start_node
        set_active_node(0, &start_node);
        set_active_node(1, &second_node);
    }

    // Start iterator from rec_point.base_node(), possibly several nodes
    // below recombination
    for (TimeIntervalIterator ti(this, active_node(0)->first_child()); ti.good(); ++ti) {

        // Skip intervals before first event;
        // remain processing this TimeInterval until no new events occur within it.
        while ( tmp_event_.time() < (*ti).end_height() ) {

            // Update States & Rates (see their declaration for explanation);
            states_[0] = getNodeState(active_node(0), tmp_event_.time());
            states_[1] = getNodeState(active_node(1), tmp_event_.time());

            // Fixed time events (e.g pop splits/merges & single migration events first
            if (model().hasFixedTimeEvent( tmp_event_.time() ))
                dontImplementFixedTimeEvent(ti);

            // Calculate the rates of events in this time interval
            TimeInterval interval( &ti, tmp_event_.time(), (*ti).end_height() );
            calcRates(interval);

            // Sample the time at which the next event happens (if any)
            // If no event, tmp_event_time_ and tmp_event_.time() are set to -1
            sampleEvent(interval, tmp_event_time_, tmp_event_);

            // Implement the event
            if ( tmp_event_.isNoEvent() ) {
                this->dontImplementNoEvent(*ti, coalescence_finished_);
            }

            else if ( tmp_event_.isPwCoalescence() ) {
                last_coal_height_ = tmp_event_time_;
                if (first_event_height_ < 0.0) first_event_height_ = tmp_event_time_;
                if (first_coal_height_ < 0.0)  first_coal_height_  = tmp_event_time_;
		tmp_event_time_ = primary_root()->height();          // Disable buffer for next genealogy.
                this->coalescence_finished_ = true;                  // we're done
            }

            else if ( tmp_event_.isRecombination() ) {
                this->dontImplementRecombination(tmp_event_, ti);
            }

            else if ( tmp_event_.isMigration() ) {
                if (first_event_height_ < 0.0) first_event_height_ = tmp_event_time_;
                tmp_event_.node()->set_population(tmp_event_.mig_pop());  // implement
            }

            else if ( tmp_event_.isCoalescence() ) {
                last_coal_height_ = tmp_event_time_;
                if (first_event_height_ < 0.0) first_event_height_ = tmp_event_time_;
                if (first_coal_height_ < 0.0)  first_coal_height_  = tmp_event_time_;
                this->dontImplementCoalescence(tmp_event_, ti);
            }

            if (coalescence_finished()) {
		tmp_event_time_ = primary_root()->height(); // Disable buffer for next genealogy.
		current_rec_--;   // undo the advance above
                return importance_weight;
            }
        }
    }
    throw std::logic_error("No final coalescence event was sampled!");
}


void ForestState::dontImplementFixedTimeEvent(TimeIntervalIterator &ti) {

    double sample;
    bool migrated;
    size_t chain_cnt, pop_number = model().population_number();

    for (size_t i = 0; i < 2; ++i) {
        if (states_[i] != 1) continue;
        chain_cnt = 0;
        while (true) {
            migrated = false;
            sample = random_generator()->sample();
            for (size_t j = 0; j < pop_number; ++j) {
                sample -= model().single_mig_pop(active_node(i)->population(), j);
                if (sample < 0) {
                    tmp_event_ = Event((*ti).start_height());
                    tmp_event_.setToMigration(active_node(i), i, j);
                    tmp_event_.node()->set_population(tmp_event_.mig_pop());    // implement
                    migrated = true;
                    break;
                }
            }

            // Stop if no migration occurred
            if (!migrated) break;
            
            // Resolve a maximum of 10k chained events for each node
            if (chain_cnt == 10000) throw std::logic_error("Cycle detected when moving individuals between populations");
            ++chain_cnt;
        }
    }
}


void ForestState::dontImplementNoEvent(const TimeInterval &ti, bool &coalescence_finished) {

    // set start-of-interval to end height
    tmp_event_.set_time( ti.end_height() );

    if (ti.end_height() == DBL_MAX) throw std::logic_error("Lines did not coalesce.");
    if (states_[0] == 2) {
        set_active_node(0, virtualPossiblyMoveUpwards(active_node(0), ti));
        if (active_node(0)->local()) {
	    coalescence_finished = true;
	    tmp_event_time_ = primary_root()->height(); // Disable buffer for next genealogy.
	    return;
        }
    }
    
    // There are no local node above the local root, which is the lowest node
    // that active_node(1) can be.
    if (states_[1] == 2) set_active_node(1, virtualPossiblyMoveUpwards(active_node(1), ti));
    
    if (active_node(0)->first_child() == active_node(1)->first_child()) {
        coalescence_finished = true;
	tmp_event_time_ = primary_root()->height(); // Disable buffer for next genealogy.
    }
}


void ForestState::dontImplementRecombination(const Event &event, TimeIntervalIterator &ti) {

    Node* recombining_node = event.node();
    recombining_node->make_local();                // it is now coalescing
    recombining_node->set_parent(NULL);            // make root
    recombining_node->set_height( event.time() );
}


void ForestState::dontImplementCoalescence(const Event &event, TimeIntervalIterator &tii) {

  // Coalescence: sample target point and implement the coalescence
  Node* coal_node = event.node();
  Node* target = contemporaries_.sample(coal_node->population());

  coal_node->set_parent( target->parent() );
  coal_node->make_local();
  coal_node->make_nonlocal( target->last_update() );
  coal_node->set_first_child( target );
  // not necessary to update population, or set active node

  if ( getOtherNodesState() == 2 ) {
      // If the coalescing node coalesced into the branch directly above
      // a recombining node, we are done.  The test here is a bit different
      // from the test in forest.cc, because we don't implement the event,
      // so we test for equality of the base nodes in the actual tree
      if ( getOtherNode()->first_child() == getEventNode()->first_child() ) {
          coalescence_finished_ = true;
	  tmp_event_time_ = primary_root()->height(); // Disable buffer for next genealogy.
          return;
      }
  }

  if ( target->local() ) {
    // Only active_node(0) can coalescence into local nodes. active_node(1) is
    // at least the local root and hence above all local nodes.
    // If active_node(0) coalescences into the local tree, there are no more
    // active nodes and we are done.
    tmp_event_time_ = primary_root()->height(); // Disable buffer for next genealogy.
    coalescence_finished_ = true;
  }
}


/**
 * Helper function for doing a coalescence.
 * Moves the 'active' flag (i.e. the node stored in root_1 or root_2 in sampleCoalescence)
 * from a node to it's parent if the branch above the node
 * ends this the current time interval.
 *
 * This function is used the pass the active flag upwards in the tree if the
 * node is active, but neither coalescing nor a recombination happens on the
 * branch above, e.g. after a local-branch became active because it was hit by a
 * coalescence or a non-local branch was active and no recombination occurred.
 *
 * Also updates the active node if it moves up.
 *
 * \param node An active node
 * \param time_interval The time interval the coalescence is currently in.
 *
 * \return  Either the parent of 'node' if we need to move upwards or 'node'
 *          itself
 */
Node* ForestState::virtualPossiblyMoveUpwards(Node* node, const TimeInterval &time_interval) {
  if ( node->parent_height() == time_interval.end_height() ) {
      Node* newnode = node->parent();
      node->set_height( newnode->height() );
      node->set_population( newnode->population() );
      node->make_local();
      node->make_nonlocal( newnode->last_update() );
      node->set_first_child( newnode );
      node->set_parent( newnode->parent() );
  }
  return node;
}
