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

#include <climits> // INT_MAX

#include "particle.hpp"
#include "pfparam.hpp"


void ForestState::making_copies( int number_of_copies ){
    for ( size_t i = 0; i < this->ForestState_copies.size(); i++){
        this->ForestState_copies[i] = NULL;
    }
    this->ForestState_copies.clear(); // clear previous pointers to the new copies
    if ( number_of_copies == 0 ) return;
    this->ForestState_copies.push_back( this );
    if ( number_of_copies == 1 ) return;
    
    for (int ii = 2; ii <= number_of_copies; ii++) { 
        ForestState* new_copy_state = new ForestState( *this );
        // Give the new particle its own random generator (for multithreading)
        size_t new_seed = (size_t) this->random_generator_->sampleInt( INT_MAX );
        new_copy_state->random_generator_ = new MersenneTwister( new_seed , this->random_generator_->ff() ); 
        new_copy_state->model_ = new Model( *this->model_);
        this->ForestState_copies.push_back(new_copy_state);
    }
    assert( this->ForestState_copies.size() == (size_t)number_of_copies ); 
}



/*! \brief Initialize a new ForestState 
 * @ingroup group_pf_init
 * */
ForestState::ForestState( Model* model, RandomGenerator* random_generator, const vector<int>& record_event_in_epoch)
            :Forest( model, random_generator ),
             record_event_in_epoch(record_event_in_epoch) {
    /*! Initialize base of a new ForestState, then do nothing, other members will be initialized at an upper level */
    //this->init_EventContainers( model );    
	//this->buildInitialTree();    
    //cout << "this->TmrcaHistory.size() = "<<this->TmrcaHistory.size()<<endl;
    // Initialize the initial particles at base 0, with equal weights
    this->setParticleWeight( 1.0 );
    this->setSiteWhereWeightWasUpdated(0.0);
    new_forest_counter++; // DEBUG
}


/*! \brief Create a newly copied ForestState 
    @ingroup group_pf_resample 
*/
ForestState::ForestState( const ForestState & copied_state )
            :Forest( copied_state ),
             record_event_in_epoch( copied_state.record_event_in_epoch ) {
    this->setParticleWeight( copied_state.weight() );
    this->setSiteWhereWeightWasUpdated( copied_state.site_where_weight_was_updated() );    
    this->setAncestor ( copied_state.ancestor() );
    this->copyEventContainers ( copied_state );
    this->segment_count_ = copied_state.segment_count_;
    
    for (size_t i = 0 ; i < copied_state.TmrcaHistory.size(); i++ ){
        this->TmrcaHistory.push_back(copied_state.TmrcaHistory[i]);
    }        
	dout << "current particle's weight is " << this->weight()<<endl;
    new_forest_counter++;  // DEBUG
}
    
    

void ForestState::copyEventContainers(const ForestState & copied_state ){

	// Copy trees
	for (size_t i=0; i < copied_state.eventTrees.size() ; i++) {
		EvolutionaryEvent* new_event = copied_state.eventTrees[i];
		eventTrees.push_back( new_event );
		if (new_event) 
			new_event->increase_refcount();
	}
}



void ForestState::init_EventContainers( Model * model ) {

    for (size_t i=0; i < model->change_times_.size(); i++) {
		eventTrees.push_back( NULL );
	}
}



/*! \brief Destructor of ForestState
 * Recursively remove all the previous states, if the pointer counter is zero
 */
ForestState::~ForestState(){

	this->clear_eventContainer();
    delete this->random_generator_;     //MULTITRHREADING
    delete this->model_;
    delete_forest_counter++;
	dout << "A Foreststate is deleted" << endl;
}


/*! Clear coalescent and recombination events recorded between two states.*/
void ForestState::clear_eventContainer(){ 

    for (int i = eventTrees.size()-1 ; i>=0 ; --i) {
		if (eventTrees[i] && eventTrees[i]->decrease_refcount_is_zero()) {
			// We use placement new, so need to call destructor explicitly.
			// However the destructor should recursively delete its parents,
			// and therefore must know the epoch -- but we can't pass parameters.
			// So, call a helper deleter, that can take a parameter and both
			// destructs and deallocates the memory.
			/*
			 * delete eventTrees[i];  // this recursively deletes its parents
			 */
			eventTrees[i]->deletethis( i );  // this recursively deletes its parents
		}
	}
	//eventTrees.clear();
}



void ForestState::record_all_event(TimeInterval const &ti, double &recomb_opp_x_within_scrm){

    dout << endl;
    ForestStatedout << "Start Recording " << endl;
    double recomb_opp_x_within_smcsmc = 0; // DEBUG
	double start_base, end_base;
	double start_height, end_height;
	size_t start_height_epoch, end_height_epoch;

    // find extent of time interval where an event may have occurred
    start_height = ti.start_height();
    start_height_epoch = ti.forest().model().current_time_idx_;
    //assert( start_height_epoch == ti.forest().model().getTimeIdx( start_height ) );
    if (this->tmp_event_.isNoEvent()) {
		end_height = start_height + ti.length();
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
            if (!(record_event_in_epoch[ writable_model()->current_time_idx_ ] & PfParam::RECORD_RECOMB_EVENT)) continue;
            start_base = active_node(i)->last_update();
            end_base = this->current_base();
			if (end_base == start_base) continue;
			//EvolutionaryEvent* recomb_event = new EvolutionaryEvent( start_height, start_height_epoch, end_height, end_height_epoch, start_base, end_base, 1 );
			// create event (on stack), so that we can append events to existing ones.  (Only done for recombinations; should check whether it happens at any frequency)
			void* event_mem = Arena::allocate( start_height_epoch );
			EvolutionaryEvent* recomb_event = new(event_mem) EvolutionaryEvent( start_height, start_height_epoch, end_height, end_height_epoch, start_base, end_base, 1 );
			double recomb_pos = -1;
			if (tmp_event_.isRecombination() && tmp_event_.active_node_nr() == i) {
				recomb_pos = (start_base + end_base)/2;   // we should sample from [start,end], but the data isn't used
				recomb_event->set_recomb_event_pos( recomb_pos );
			}
			// add event in tree data structure
			recomb_event->add_leaf_to_tree( &eventTrees[ writable_model()->current_time_idx_] );
            recomb_opp_x_within_smcsmc += end_base - start_base;
		} else if (states_[i] == 1) {
            // node i is tracing out a new branch; opportunities for coalescences and migration
            if (!(record_event_in_epoch[ writable_model()->current_time_idx_ ] & PfParam::RECORD_COALMIGR_EVENT)) continue;
            start_base = this->current_base();
            int weight = ti.numberOfContemporaries( active_node(i)->population() );
			// consider normal (not pairwise) coalescences that occurred on this node
	        bool coal_event = (tmp_event_.isCoalescence() && tmp_event_.active_node_nr() == i);
	        bool migr_event = (tmp_event_.isMigration() && tmp_event_.active_node_nr() == i);
	        // account for potential pairwise coalescence opportunity and event
	        if (i==0 && states_[1]==1 && active_node(0)->population() == active_node(1)->population()) {
				weight += 1;
				coal_event |= tmp_event_.isPwCoalescence();
			}
			// Record coalescence and migration opportunity (note: if weight=0, there is still opportunity for migration)
			void* event_mem = Arena::allocate( start_height_epoch );
			EvolutionaryEvent* migrcoal_event = new(event_mem) EvolutionaryEvent(start_height, start_height_epoch, end_height, end_height_epoch, start_base, active_node(i)->population(), weight);
			// Record any events
			if (coal_event) migrcoal_event->set_coal_event();
			if (migr_event) migrcoal_event->set_migr_event( tmp_event_.mig_pop() );
			// add event in tree data structure
			migrcoal_event->add_leaf_to_tree( &eventTrees[ writable_model()->current_time_idx_] );
		}
	}
	
	assert ( recomb_opp_x_within_smcsmc == recomb_opp_x_within_scrm );
    dout << endl;        

    return;
}



void ForestState::record_Recombevent_b4_extension (){ 
    #ifdef _RecombRecordingOff // if the recombination recording off macro is defined, then return without recording the event
        return;
    #endif
    
    dout << endl;
    ForestStatedout << " Build new genealogy, compute the recombination opportunity at each level " << endl;

    // iterate over time intervals (but do NOT prune branches at this stage)
    for (TimeIntervalIterator ti(this, this->nodes_.at(0), false); ti.good(); ++ti) {
        ForestStatedout << " * * Time interval: " << (*ti).start_height() << " - " << (*ti).end_height() << " " ;
        dout << ", with " << ti.numberOfLocalContemporaries() << " local Contemporaries, " << ti.numberOfLocalContemporaries() << " * " << "( " << (*ti).end_height() << " - " << (*ti).start_height() << ")" << std::endl;
        // Create a recombination event for this slice (which may be smaller than an epoch -- but in our case it usually won't be)
        int contemporaries = ti.numberOfLocalContemporaries();
        if (contemporaries > 0 && (record_event_in_epoch[ writable_model()->current_time_idx_ ] & PfParam::RECORD_RECOMB_EVENT)) {
            double start_height = (*ti).start_height();
            double end_height = (*ti).end_height();
            size_t start_height_epoch = ti.forest().model().current_time_idx_;
            //assert( start_height_epoch == ti.forest().model().getTimeIdx( start_height ) );
            size_t end_height_epoch = start_height_epoch;
            void* event_mem = Arena::allocate( start_height_epoch );
            EvolutionaryEvent* recomb_event = new(event_mem) EvolutionaryEvent( start_height, start_height_epoch, end_height, end_height_epoch, this->current_base(), this->next_base_, contemporaries );  // no event for now
            recomb_event->add_leaf_to_tree( &eventTrees[ writable_model()->current_time_idx_] );
        }
    }
}


void ForestState::record_Recombevent_atNewGenealogy ( double event_height ){
    #ifdef _RecombRecordingOff // if the recombination recording off macro is defined, then return without recording the event
        return;
    #endif
    recombination_counter++; // DEBUG
    this->writable_model()->resetTime( event_height );
    size_t epoch_i = this->writable_model()->current_time_idx_;
    if (!(record_event_in_epoch[ epoch_i ] & PfParam::RECORD_RECOMB_EVENT)) 
        return;
    dout << "current_time_idx_ =  " << epoch_i << " = [" << this->writable_model()->getCurrentTime() << "," << this->writable_model()->getNextTime() << endl;
    // find the EvolutionaryEvent to add this event to.
	EvolutionaryEvent* event = eventTrees[ epoch_i ];
	while ( !event->is_recomb() || !event->recomb_event_overlaps_opportunity( event_height ) ) {
		event = event->parent();
		assert (event != NULL);
	}
	assert (event->start_base() == this->current_base());
	event->set_recomb_event_time( event_height );
}




void ForestState::include_haplotypes_at_tips(vector <int> &haplotypes_at_tips){
	for (size_t j=0; j < haplotypes_at_tips.size();j++){		
		this->nodes()->at(j)->set_mutation_state(haplotypes_at_tips[j]);
	    }
    }


/*! 
 * Calculate the marginal likelihood of a node recursively, let X denote the state of Node *, and Y denote the state of the 
 *  first child, and Z denote the state of the second child. Let t1 denote the time from X to Y, and t2 is the time from X to Z.
 *  Let u(t) be the probability function of the no mutation occurs within time t. 
 * 
  \verbatim
      X
   t1/ \ t2
    Y   Z
  \endverbatim
 * 
 * 	Suppose that X, Y and Z only take values 0 or 1.
 *  When X=0, consider the likelihood of the tree in four cases of (y,z) pairs (0,0), (0,1), (1,0) and (1,1)
 *  Similarliy X=1, we have four cases too. 
 * X[0]=y[0] * ut1 * z[0] * ut2 + y[0] * ut1 * z[1] * (1-ut2) + y[1] * (1-ut1) * z[0] * ut2 + y[1] * (1-ut1) * z[1] * (1-ut2);
 * X[1]=y[1] * ut1 * z[1] * ut2 + y[1] * ut1 * z[0] * (1-ut2) + y[0] * (1-ut1) * z[1] * ut2 + y[0] * (1-ut1) * z[0] * (1-ut2);
 * 
 * After simplification, 
 * X[0]=(y[0]*ut1 + y[1]*(1-ut1)) * (z[0]*ut2 + z[1]*(1-ut2)) 
 * X[1]=(y[1]*ut1 + y[0]*(1-ut1)) * (z[1]*ut2 + z[0]*(1-ut2)) 
 * 
 * Note: even though there are only two states, x[0]+x[1] is not 1!!! becasue x[0] is a marginal probability, but still conditional on all the possiblities of the variants 
 *  
 *  If X, Y and Z take values A, T, C and G. The computation would be more complex.  Four equations, with 16 terms in each.
 * 
 * If we consider JC69 model, then 
 * X[0]= y[0] * ut1 * z[0] * ut2  + y[0] * ut1 * sum(z[1,2,3]) * (1-ut2) +  sum(y[1,2,3]) * (1-ut1) * z[0] * ut2 + sum(y[1,2,3]) * (1-ut1) * sum(z[1,2,3]) * (1-ut2) 
 * ...
 * 
 * @ingroup group_resample
 * * */


inline valarray<double> ForestState::cal_partial_likelihood_infinite(Node * node) {

	valarray<double> part_lik(2);    // partial likelihood of data under this node conditional on current node being in state 0 or 1

	// deal with the case that this node is a leaf node
	if ( node->first_child() == NULL ) {
	  assert ( node->second_child() == NULL );          // if a node has no first child, it won't have a second child
	  assert ( node->in_sample() );                     // we only traverse the local tree, therefore the leaf node must be in our sample
	  part_lik[0] = node->mutation_state() == 1 ? 0.0 : 1.0;	// also encode state==-1 (missing data) as 1.0
	  part_lik[1] = node->mutation_state() == 0 ? 0.0 : 1.0;   // also encode state==-1 (missing data) as 1.0
	  return part_lik;
    }

    // Genealogy branch lengths are in number of generations, the mutation rate is unit of per site per generation, often in the magnitude of 10 to the power of negative 8.
	double mutation_rate = this->model().mutation_rate();
	Node *left = trackLocalNode(node->first_child());   // jz_stable; for jz use Node *left = node->getLocalChild1();
	Node *right = trackLocalNode(node->second_child()); // jz_stable; for jz use Node *right = node->getLocalChild2();
	double t_left = node->height() - left->height();
	double t_right = node->height() - right->height();
	assert (t_left >= 0 && t_right >= 0 && mutation_rate > 0);
	double p_left = exp(-t_left * mutation_rate);   // probability that no mutation occurred along left branch.  Assume infinite sites model
	double p_right = exp(-t_right * mutation_rate); // probability that no mutation occurred along left branch.  Assume infinite sites model
	valarray<double> part_lik_left = cal_partial_likelihood_infinite(left);     // partial likelihoods of data under left and right children
	valarray<double> part_lik_right = cal_partial_likelihood_infinite(right);

	part_lik[0] = ( part_lik_left[0]*p_left + part_lik_left[1]*(1-p_left) ) * ( part_lik_right[0]*p_right + part_lik_right[1]*(1-p_right) );
	part_lik[1] = ( part_lik_left[1]*p_left + part_lik_left[0]*(1-p_left) ) * ( part_lik_right[1]*p_right + part_lik_right[0]*(1-p_right) );

	return part_lik;

}

	
/*! 
 * \brief Calculate the marginal likelihood of the genealogy at data site i, 
 *  If there is no data given at the site i, return likelihood as 1. 
 * @ingroup group_pf_resample
 */
double ForestState::calculate_likelihood( ) {
    //dout << "calculate_likelihood function, root is " << this->local_root() << endl;
    valarray<double> marginalLikelihood = cal_partial_likelihood_infinite(this->local_root());
    dout << "marginal likelihood is " << marginalLikelihood[0] <<  "," << marginalLikelihood[1] << endl;
    double prior[2] = {0.5,0.5};
    double likelihood = marginalLikelihood[0]*prior[0] + marginalLikelihood[1]*prior[1];
    dout << "ForestState::calculate_likelihood( ) likelihood is " << likelihood << endl;
    return likelihood;
}


double ForestState::trackLocalTreeBranchLength() {
    BranchLengthData bld = trackSubtreeBranchLength( this->local_root() );
    if (bld.subtreeBranchLength == -1) {
		// none of the leaves carry data -- total length is 0
		return 0;
	}
	// return branch length of the subtree subtending leaves carrying data
	return bld.subtreeBranchLength;
}


BranchLengthData ForestState::trackSubtreeBranchLength ( Node * currentNode ) {

	if (currentNode->in_sample() ) {
		// current node is a leaf node
		if (currentNode->mutation_state() >= 0) {
			// leaf node carries data
			return BranchLengthData( 0, 0 );
		} else {
			// leaf node carries no data
			return BranchLengthData( 0, -1 );
		}
	}

	Node* left_local_child = trackLocalNode(currentNode->first_child());
	Node* right_local_child = trackLocalNode(currentNode->second_child());
	
    BranchLengthData bld_left  = this->trackSubtreeBranchLength( left_local_child );
    BranchLengthData bld_right = this->trackSubtreeBranchLength( right_local_child );

	// calculate branch length of partial tree, including the branch from this node to the child node
    double leftBL = bld_left.partialBranchLength + (currentNode->height() - left_local_child->height());
    double rightBL = bld_right.partialBranchLength + (currentNode->height() - right_local_child->height());

	// return correct partial tree branch length, and subtree branch length.  The calculation depends on
	// whether left and right subtrees carry data or not.
    if (bld_left.subtreeBranchLength >= 0 && bld_right.subtreeBranchLength >= 0)
		// both left and right subtrees carry data, so current node is a possible root node
        return BranchLengthData( leftBL+rightBL, leftBL+rightBL );

    if (bld_left.subtreeBranchLength >= 0)
		// left subtree carries data, but right one doesn't -- keep left root as possible root node
		return BranchLengthData( leftBL, bld_left.subtreeBranchLength );

    if (bld_right.subtreeBranchLength >= 0)
		// same for right subtree
		return BranchLengthData( rightBL, bld_right.subtreeBranchLength );
		
	// neither subtree contains data, so just return the length data of either
	return bld_left;
}


double ForestState::extend_ARG ( double mutation_rate, double extend_to, Segment_State segment_state, bool updateWeight, bool recordEvents ) {

    double updated_to = this->site_where_weight_was_updated();
    dout << "Particle current base is at " << this->current_base() << " weight is updated to " << updated_to <<endl;
    assert (updated_to >= this->current_base());
    double likelihood = 1.0;
    
    while ( updated_to < extend_to ) {
        
        dout << "  Now at " <<this->current_base()<< " updated_to " << updated_to << " and extending to " << extend_to << endl;            
        /*!
         * First, update the likelihood up to either extend_to or the end of this state
         */
        double update_to = min( extend_to, this->next_base() );
        // calculate the total tree length of the subtree over leaf nodes that carry data.  
        double localTreeBranchLength = this->trackLocalTreeBranchLength();
        // for leaf nodes that carry no data, there is no evidence for presence or absence of mutations on corresponding branches.
		// This is accounted for in the calculation of localTreeBranchLength.  In addition, a segment can be explicitly marked as
		// having 'missing data' (segment_state != SEGMENT_INVARIANT), in which case ALL branches are considered uninformative.
        assert ( (segment_state == SEGMENT_INVARIANT) || (localTreeBranchLength == 0 ));
        double likelihood_of_segment = exp( -mutation_rate * localTreeBranchLength * (update_to - updated_to) );
        
        dout << " Likelihood of no mutations in segment of length " << (update_to - updated_to) << " is " << likelihood_of_segment ;
        dout << ( ( segment_state == SEGMENT_INVARIANT ) ? ", as invariant.": ", as missing data" ) << endl;
        likelihood *= likelihood_of_segment;
        updated_to = update_to;  // rescues the invariant
        /*!
         * Next, if we haven't reached extend_to now, add a new state and iterate
         */
        if ( updated_to < extend_to ) {
            this->sampleNextGenealogy( recordEvents );
        }        
    }
    assert (updated_to == extend_to);        
    if (updateWeight) {
        this->setParticleWeight( this->weight() * likelihood );
    }
    this->setSiteWhereWeightWasUpdated( extend_to );
    return likelihood;
}


std::string ForestState::newick(Node *node) {
  if(node->in_sample()){
    std::ostringstream label_strm;
    label_strm<<node->label();
    return label_strm.str();
  }
  else{
    Node *left = this->trackLocalNode(node->first_child());
    double t1 = node->height() - left->height();
    std::ostringstream t1_strm;
    t1_strm << t1 / (4 * this->model().default_pop_size);

    Node *right = this->trackLocalNode(node->second_child());
    double t2 = node->height() - right->height();
    std::ostringstream t2_strm;
    t2_strm << t2 / (4 * this->model().default_pop_size);

    return "("+this->newick(left)+":"+t1_strm.str()+","+ this->newick(right)+":"+t2_strm.str() +")";
  }
}
