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

#include"particle.hpp"
#include <climits> // INT_MAX

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
        //new_copy_state->random_generator_ = new MersenneTwister( new_seed ); 
        new_copy_state->model_ = new Model( *this->model_);
        //dout << "new random seed is " << this->particles[old_state_index]->random_generator_->seed() + this->particles.size() << endl;
        this->ForestState_copies.push_back(new_copy_state);
        }
    assert( this->ForestState_copies.size() == (size_t)number_of_copies ); 
    }



/*! \brief Initialize a new ForestState 
 * @ingroup group_pf_init
 * */
ForestState::ForestState( Model* model, RandomGenerator* random_generator )
            :Forest( model, random_generator ) {
    /*! Initialize base of a new ForestState, then do nothing, other members will be initialized at an upper level */
    //this->init_EventContainers( model );    
	//this->buildInitialTree();    
    //cout << "this->TmrcaHistory.size() = "<<this->TmrcaHistory.size()<<endl;
    // Initialize the initial particles at base 0, with equal weights
  	//this->setParticleWeight( 1.0 );
	//this->setSiteWhereWeightWasUpdated(0.0);
    //new_forest_counter++; // DEBUG
    }


/*! \brief Create a newly copied ForestState 
    @ingroup group_pf_resample 
*/
ForestState::ForestState( const ForestState & copied_state )
            :Forest( copied_state ) {
    this->setParticleWeight( copied_state.weight() );
	this->setSiteWhereWeightWasUpdated( copied_state.site_where_weight_was_updated() );    
    this->setAncestor ( copied_state.ancestor() );
    this->copyEventContainers ( copied_state );
    this->segment_count_ = copied_state.segment_count_;
    
    for (size_t i = 0 ; i < copied_state.TmrcaHistory.size(); i++ ){
        this->TmrcaHistory.push_back(copied_state.TmrcaHistory[i]);
        }        
	dout << "current particle's weight is " << this->weight()<<endl;
    //new_forest_counter++;  // DEBUG
    }
    
    

void ForestState::copyEventContainers(const ForestState & copied_state ){
    // Copy Coalescent events
    for (size_t i = 0 ; i < copied_state.CoaleventContainer.size() ; i++ ){ 
        deque < Coalevent* > CoaleventContainer_i;   
        for (size_t ii = 0 ; ii < copied_state.CoaleventContainer[i].size(); ii++){
            Coalevent* new_coalevent = copied_state.CoaleventContainer[i][ii];
            new_coalevent->pointer_counter_++;
            CoaleventContainer_i.push_back ( new_coalevent ) ;
            }
        CoaleventContainer.push_back( CoaleventContainer_i );        
        }
    // Copy Recombination events
    for (size_t i = 0 ; i < copied_state.RecombeventContainer.size() ; i++ ){ 
        deque < Recombevent* > RecombeventContainer_i;  
        for (size_t ii = 0 ; ii < copied_state.RecombeventContainer[i].size(); ii++){
            Recombevent* new_recombevent = copied_state.RecombeventContainer[i][ii];
            new_recombevent->pointer_counter_++;
            RecombeventContainer_i.push_back (new_recombevent) ;
            }
        RecombeventContainer.push_back( RecombeventContainer_i );        
        }
    // Copy Migration events
    for (size_t i = 0 ; i < copied_state.MigreventContainer.size() ; i++ ){ 
        deque < Migrevent* > MigreventContainer_i;   /*!< \brief Coalescent events recorder */
        for (size_t ii = 0 ; ii < copied_state.MigreventContainer[i].size(); ii++){
            Migrevent* new_migrevent = copied_state.MigreventContainer[i][ii];
            new_migrevent->pointer_counter_++;
            MigreventContainer_i.push_back (new_migrevent) ;
            }
        MigreventContainer.push_back( MigreventContainer_i );        
        }    
    }



void ForestState::init_EventContainers( Model * model ){
    for (size_t i = 0 ; i < model->change_times_.size() ; i++ ){ 
        deque < Coalevent*  > CoaleventContainer_i;   /*!< \brief Coalescent events recorder */
        CoaleventContainer.push_back( CoaleventContainer_i );        
        deque < Recombevent* > RecombeventContainer_i; /*!< \brief Recombination events recorder */
        RecombeventContainer.push_back( RecombeventContainer_i );
        deque < Migrevent* > MigreventContainer_i;   /*!< \brief Migration events recorder */
        MigreventContainer.push_back( MigreventContainer_i );                
        }
    }



/*! \brief Destructor of ForestState
 * Recursively remove all the previous states, if the pointer counter is zero
 */
ForestState::~ForestState(){
	this->clear_CoaleventContainer();   
	this->clear_RecombeventContainer(); 
    this->clear_MigreventContainer();    
    delete this->random_generator_;     //MULTITRHREADING
    delete this->model_;
    //delete_forest_counter++;
	dout << "A Foreststate is deleted" << endl;
    }


void ForestState::record_all_event(TimeInterval const &ti){
    double coal_opportunity = 0.0;
    double recomb_opportunity = 0.0;
    double migr_opportunity = 0.0;

    //opportunity_y is the branch length of the conterporaries within the interval
    // if there is no events, then take the full length of this time interval
    //             otherwise, take the distance between the event and the bottom of this interval
    double opportunity_y = this->tmp_event_.isNoEvent() ? ti.length() : (this->tmp_event_.time() - ti.start_height());

    for (int i = 0; i < 2; i++) {
        //if ( states_[i] == 0 ) continue; //NEW Only Nodes in state 1 or 2 can do something
        if (states_[i] == 2) {
            // node i is tracing a non-local branch; opportunities for recombination
            assert(active_node(i)->last_update() >= model().getCurrentSequencePosition());
			// If the interval [ active_node(i)->last_update(), this->current_base() ] straddles one (or more) recombination rate changes,
			// we should do something more complicated:
			// 1. sample a recombination point somewhere along the interval (adjust rates for the changing recombination rates)
			// 2. generate appropriate "no event" and "event" recombination events, with appropriate opportunities.
            recomb_opportunity += ( this->current_base() - active_node(i)->last_update() ) * opportunity_y;
            //if (this->current_base() > active_node(i)->last_update())
                //this->record_Recombevent( 0,
                                      //( this->current_base() - active_node(i)->last_update() ) * opportunity_y,
                                      //( tmp_event_.isRecombination() && ( tmp_event_.node() == active_node(i) ) ) ? EVENT : NOEVENT, 
                                      //active_node(i)->last_update()); // DEBUG 0 for now, this is tricky, as the last update is base different
            dout << "Tracing non-local branch along " << opportunity_y << " generations, from " << active_node(i)->last_update() << " to " << this->current_base() << endl;
            }
        else if (states_[i] == 1) {
            // node i is tracing out a new branch; opportunities for coalescences and migration
            coal_opportunity += ti.numberOfContemporaries( active_node(i)->population() ) * opportunity_y; // jz_stable
            //coal_opportunity += contemporaries_.size( active_node(i)->population() ) * opportunity_y; // jz
            migr_opportunity += opportunity_y;
            }
        }

    // only coalescences into contemporaries were considered; pairwise coalescence between active nodes could also occur
    if ((states_[0] == 1) && (states_[1] == 1) && (active_node(0)->population() == active_node(1)->population() ) ) {
        coal_opportunity += opportunity_y;
        }
        
    if (coal_opportunity > 0) {
        //this->record_event(active_node(0)->population(), size_t(-1), ti.start_height(), ti.start_height() + opportunity_y, coal_opportunity, (tmp_event_.isCoalescence() || tmp_event_.isPwCoalescence()) ? COAL_EVENT : COAL_NOEVENT );
        //// check: do we need to consider different cases? is it possible to coalesce in active_node(1)->population() ???? Joe: I dont think so ...
        if ((states_[0] == 1) && (states_[1] == 1) && (active_node(0)->population() == active_node(1)->population() ) ) {
            this->record_Coalevent(active_node(0)->population(), 
                                    //ti.start_height(), 
                                    //ti.start_height() + opportunity_y, 
                                    coal_opportunity, 
                                    (tmp_event_.isCoalescence() || tmp_event_.isPwCoalescence()) ? EVENT : NOEVENT );
            } 
        else if (states_[0] == 1){
            this->record_Coalevent(active_node(0)->population(), 
                                    //ti.start_height(), 
                                    //ti.start_height() + opportunity_y, 
                                    coal_opportunity, 
                                    (tmp_event_.isCoalescence() || tmp_event_.isPwCoalescence()) ? EVENT : NOEVENT ); 
            } 
        else if (states_[1] == 1){
            this->record_Coalevent(active_node(1)->population(), 
                                    //ti.start_height(), 
                                    //ti.start_height() + opportunity_y, 
                                    coal_opportunity, 
                                    (tmp_event_.isCoalescence() || tmp_event_.isPwCoalescence()) ? EVENT : NOEVENT );
            }
        }

    if (migr_opportunity > 0 && this->model().population_number()>1) {
        if ( tmp_event_.isMigration() ){
            this->record_Migrevent(tmp_event_.node()->population(), 
                                    //ti.start_height(), 
                                    //ti.start_height() + opportunity_y, 
                                    migr_opportunity, EVENT , tmp_event_.mig_pop());
            } 
        else {
            this->record_Migrevent(active_node(0)->population(),                                    
                                    //ti.start_height(), 
                                    //ti.start_height() + opportunity_y, 
                                    migr_opportunity, NOEVENT, size_t(-1) );
            }
        }
    
    
    // Record recombination event count and opportunity separately
    
    if (recomb_opportunity > 0) {
        // do not store the population, as the recomb_opportunity is also calculated 
        // over the full tree rather than per-population
        this->record_Recombevent( 0, recomb_opportunity, tmp_event_.isRecombination() ? EVENT : NOEVENT, tmp_event_.node()->last_update() );
        //if (tmp_event_.isRecombination()) {
            //this->record_Recombevent( 0,
                                     ////ti.start_height(), 
                                     ////ti.start_height() + opportunity_y, 
                                     //recomb_opportunity,
                                     //EVENT, tmp_event_.node()->last_update()); 
        //} else {
            //this->record_Recombevent( 0,
                                      ////ti.start_height(), 
                                      ////ti.start_height() + opportunity_y, 
                                      //recomb_opportunity,
                                      //NOEVENT, tmp_event_.node()->last_update()); 
        //}
    }
    return;
}


/*! \brief Record Coalescent events
* @ingroup group_count_coal
*/
void ForestState::record_Coalevent(
                  size_t pop_i,
                  //double start_time, 
                  //double end_time, 
                  double opportunity, 
                  eventCode event_code) {
    //cout<<"current_time_index = " << this->writable_model()->current_time_idx_<<endl;
    Coalevent * new_event = new Coalevent( pop_i,
                          //start_time,
                          //end_time, 
                          opportunity,
                          event_code);	
	new_event->set_epoch_index ( this->writable_model()->current_time_idx_ ) ;
    new_event->set_end_base ( this->current_base() );
    assert(new_event->print_event());
    this->CoaleventContainer[this->writable_model()->current_time_idx_].push_back(new_event);
    
    }


/*! \brief Record Recombination events
* @ingroup group_count_coal
*/
void ForestState::record_Recombevent(size_t pop_i,
                          //double start_time, 
                          //double end_time, 
                          double opportunity, 
                          eventCode event_code,
                          double base){
    Recombevent* new_event = new Recombevent( pop_i,
                          //start_time,
                          //end_time, 
                          opportunity,
                          event_code,
                          base);
	new_event->set_epoch_index ( this->writable_model()->current_time_idx_ ) ;
    new_event->set_end_base ( this->current_base() );
    //if (event_code==EVENT){recombination_counter++;} // DEBUG
    assert(new_event->print_event());
    this->RecombeventContainer[this->writable_model()->current_time_idx_].push_back(new_event);
    }
    
void ForestState::compute_opportunity_y_s ( ){
    this->opportunity_y_s.clear();
    this->writable_model()->resetTime( 0 );
    
    for (TimeIntervalIterator ti(this, this->nodes_.at(0)); ti.good(); ++ti) {
        dout << "* * Time interval: " << (*ti).start_height() << " - "
             << (*ti).end_height() << " " ;
        // Need to consider multiple populations
        dout << ti.numberOfLocalContemporaries() << " local Contemporaries"   << std::endl;
        this->writable_model()->resetTime( (*ti).start_height() );
        if ( (this->writable_model()->current_time_idx_ ) >= this->opportunity_y_s.size() ){
            this->opportunity_y_s.push_back( ti.numberOfLocalContemporaries() * ( (*ti).end_height() - (*ti).start_height() ) );
        } else{
            this->opportunity_y_s[this->writable_model()->current_time_idx_] += ti.numberOfLocalContemporaries() * ( (*ti).end_height() - (*ti).start_height() ) ;
        }        
    }
    
    // Checking 
    double total_bl = 0;
    for (size_t i =0 ; i < opportunity_y_s.size(); i++){
        total_bl +=  this->opportunity_y_s[i];
        }
    
    dout << "total_bl " << total_bl<< ", this->local_tree_length() =" << this->local_tree_length()<<endl;
    assert ( (total_bl - this->local_tree_length() ) < 1e-8);
    assert ( opportunity_y_s.size() == this->model().change_times_.size() );
}


void ForestState::record_Recombevent_atNewGenealogy ( double recomb_opportunity_x, bool record_event ){
    double start_base = this->current_base() - recomb_opportunity_x ;    
    for ( size_t epoch_i = 0 ; epoch_i < this->opportunity_y_s.size() ; epoch_i ++ ){
        Recombevent* new_event = new Recombevent( 0, this->opportunity_y_s[epoch_i] * recomb_opportunity_x, (epoch_i == this->writable_model()->current_time_idx_ && record_event ) ? EVENT : NOEVENT, start_base);
        new_event->set_epoch_index ( epoch_i ) ;
        new_event->set_end_base ( this->current_base() );
        assert(new_event->print_event());
        this->RecombeventContainer[this->writable_model()->current_time_idx_].push_back(new_event);
    } 
}
    
/*! \brief Record Migration events
* @ingroup group_count_coal
*/
void ForestState::record_Migrevent(size_t pop_i,
                          //double start_time, 
                          //double end_time, 
                          double opportunity, 
                          eventCode event_code, size_t mig_pop) {
    Migrevent* new_event = new Migrevent( pop_i,
                          //start_time,
                          //end_time, 
                          opportunity,
                          event_code, mig_pop );
    new_event->set_epoch_index ( this->writable_model()->current_time_idx_ ) ;
    new_event->set_end_base ( this->current_base() );
    assert(new_event->print_event());
	this->MigreventContainer[this->writable_model()->current_time_idx_].push_back(new_event);
    }    


/*! Clear coalescent and recombination events recorded between two states.*/
void ForestState::clear_CoaleventContainer(){ 
    for (size_t time_i = 0; time_i < this->CoaleventContainer.size(); time_i++ ){
        for (size_t i=0; i < this->CoaleventContainer[time_i].size(); i++){
            CoaleventContainer[time_i][i]->pointer_counter_ --;
            if (CoaleventContainer[time_i][i]->pointer_counter_ == 0){
                delete CoaleventContainer[time_i][i];
                }
            }
        }
    }


/*! Clear recombination events recorded between two states.*/
void ForestState::clear_RecombeventContainer(){ 
    for (size_t time_i = 0; time_i < this->RecombeventContainer.size(); time_i++){
        for (size_t i=0; i < this->RecombeventContainer[time_i].size(); i++){
            RecombeventContainer[time_i][i]->pointer_counter_ --;
            if (RecombeventContainer[time_i][i]->pointer_counter_ == 0){
                delete RecombeventContainer[time_i][i];
                }
            }
        }
    }
    
    
/*! Clear migration events recorded between two states.*/
void ForestState::clear_MigreventContainer(){ 
    for (size_t time_i = 0; time_i < this->MigreventContainer.size(); time_i++){
        for (size_t i=0; i < this->MigreventContainer[time_i].size(); i++){
            MigreventContainer[time_i][i]->pointer_counter_ --;
            if ( MigreventContainer[time_i][i]->pointer_counter_ == 0 ){
                delete MigreventContainer[time_i][i];
                }
            }
        }
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


inline valarray<double> ForestState::cal_marginal_likelihood_infinite(Node * node) {

	valarray<double> x(2);
	//dout << "subtree at " << node << " first child is " << node->first_child() <<" second child is " <<  node->second_child()<<endl;

	// deal with the case that this node is a leaf node
	if ( node->first_child() == NULL ) {
	  assert ( node->second_child() == NULL );     // if a node has no first child, it won't have a second child
	  assert ( node->label() > 0 );                // we only traverse the local tree, therefore the leaf node must be contemporary
	  x[0] = node->mutation_state() == 1 ? 0.0 : 1.0;	// also encode state==-1 (missing data) as 1.0
	  x[1] = node->mutation_state() == 0 ? 0.0 : 1.0;   // also encode state==-1 (missing data) as 1.0
	  //dout << "Marginal probability at " << node->label() << " is " << x[0] << "," << x[1] << endl;
	  return x;
    }

        // Genealogy branch lengths are in number of generations, the mutation rate is unit of per site per generation, often in the magnitude of 10 to the power of negative 8.
	double mutation_rate = this->model().mutation_rate();
	Node *left = trackLocalNode(node->first_child());   // jz_stable
	Node *right = trackLocalNode(node->second_child()); // jz_stable
    //Node *left = node->getLocalChild1();  // jz
    //Node *right = node->getLocalChild2(); // jz
	double t1 = node->height() - left->height();
	double t2 = node->height() - right->height();
    double ut1 = 1 - exp(-t1 * mutation_rate); // assume infinite site
	double ut2 = 1 - exp(-t2 * mutation_rate); // assume infinite site
	assert (ut1>=0 && ut1<=1);
        assert (ut2>=0 && ut2<=1);
	valarray<double> y = cal_marginal_likelihood_infinite(left);
	valarray<double> z = cal_marginal_likelihood_infinite(right);		
	x[0] = (y[0] * ut1 + y[1] * (1-ut1) ) * (z[0] * ut2 + z[1] * (1-ut2) );
	x[1] = (y[0] * (1-ut1) + y[1] * ut1 ) * (z[0] * (1-ut2) + z[1] * ut2 );
	//dout << "Marginal probability at " << node->label() << " is " << x[0] << "," << x[1] << endl;
	return x;
}

	
/*! 
 * \brief Calculate the likelihood of the genealogy at data site i, 
 *  If there is no data given at the site i, return likelihood as 1. Since all particles at this site are equally probable 
 * @ingroup group_pf_resample
 */
double ForestState::calculate_likelihood( ) {
    //dout << "calculate_likelihood function, root is " << this->local_root() << endl;
    valarray<double> marginal_likelihood = cal_marginal_likelihood_infinite(this->local_root());
    //dout << "marginal likelihood is " << marginal_likelihood[0] <<  "," << marginal_likelihood[1] << endl;
    double prior[2] = {0.5,0.5};
    double likelihood = marginal_likelihood[0]*prior[0] + marginal_likelihood[1]*prior[1];
    //dout << "likelihood is " << likelihood << endl;
    return likelihood;
}

