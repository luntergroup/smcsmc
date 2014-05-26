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

/*! \brief Create a newly copied ForestState 
    @ingroup group_pf_resample 
*/



ForestState::ForestState( const ForestState & copied_state )
            :Forest( copied_state ) {
    this->setParticleWeight( copied_state.weight() );
	this->setSiteWhereWeightWasUpdated( copied_state.site_where_weight_was_updated() );
    this->setAncestor ( copied_state.ancestor() );
    this->copyEventContainers ( copied_state );
    
    for (size_t i = 0 ; i < copied_state.TmrcaHistory.size(); i++ ){
        this->TmrcaHistory.push_back(copied_state.TmrcaHistory[i]);
        }
        
    //cout<<this->TmrcaHistory.size()<<endl;
	dout << "current particle's weight is " << this->weight()<<endl;
	//copied_state->pointer_counter++;
    new_forest_counter++;
    }
    
    

void ForestState::copyEventContainers(const ForestState & copied_state ){
    for (size_t i = 0 ; i < copied_state.CoaleventContainer.size() ; i++ ){ 
        deque < Starevent* > CoaleventContainer_i;   /*!< \brief Coalescent events recorder */
        for (size_t ii = 0 ; ii < copied_state.CoaleventContainer[i].size(); ii++){
            Starevent* new_coalevent = copied_state.CoaleventContainer[i][ii];
            new_coalevent->pointer_counter_++;
            CoaleventContainer_i.push_back ( new_coalevent ) ;
            }
        CoaleventContainer.push_back( CoaleventContainer_i );        
        }
    for (size_t i = 0 ; i < copied_state.RecombeventContainer.size() ; i++ ){ 
        deque < Starevent* > RecombeventContainer_i;   /*!< \brief Coalescent events recorder */
        for (size_t ii = 0 ; ii < copied_state.RecombeventContainer[i].size(); ii++){
            Starevent* new_recombevent = copied_state.RecombeventContainer[i][ii];
            new_recombevent->pointer_counter_++;
            RecombeventContainer_i.push_back (new_recombevent) ;
            }
        RecombeventContainer.push_back( RecombeventContainer_i );        
        }

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

/*! \brief Initialize a new ForestState 
 * @ingroup group_pf_init
 * */
ForestState::ForestState( Model* model, RandomGenerator* random_generator )
            :Forest( model,random_generator ) {/*! Initialize base of a new ForestState */
	//this->init(); // initialize members of ForestState
  	this->setParticleWeight( 1.0 );
	this->setSiteWhereWeightWasUpdated(0.0);

    this->init_EventContainers( model );
    
	this->buildInitialTree();    
    new_forest_counter++;
    }


void ForestState::init_EventContainers( Model * model ){
    for (size_t i = 0 ; i < model->change_times_.size() ; i++ ){ 
        deque < Starevent*  > CoaleventContainer_i;   /*!< \brief Coalescent events recorder */
        CoaleventContainer.push_back( CoaleventContainer_i );        
        deque < Starevent* > RecombeventContainer_i; /*!< \brief Recombination events recorder */
        RecombeventContainer.push_back( RecombeventContainer_i );
        deque < Migrevent* > MigreventContainer_i;   /*!< \brief Migration events recorder */
        MigreventContainer.push_back( MigreventContainer_i );                
        }
    }



/*! \brief Destructor of ForestState
 * Recursively remove all the previous states, if the pointer counter is zero
 */
ForestState::~ForestState(){
	dout << "State between " << this->current_base() << " and " 
                             << this->next_base() 
                             << " is about to be removed" << endl;	
	//assert( this->pointer_counter == 0 );
	this->clear_CoaleventContainer();
	this->clear_RecombeventContainer();
    this->clear_MigreventContainer();
	
    //Remove any of the previous states, if the counter is equal to zero.
	//ForestState* prior_state = this->previous_state;
	//if ( prior_state != NULL ) {
		//prior_state->pointer_counter--;
		//if ( prior_state->pointer_counter == (int)0 ){ 
            //delete prior_state; 
            //}
        //}
        
    delete_forest_counter++;
    //cout<<"Forest state destructor is called " << delete_forest_counter<<endl;    
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

    dout << " Total rate is "<< this->rates_[0] <<endl;
    for (int i=0; i<2; i++) {
        dout << "Active node["<<i<<"] has state: "<<states_[i]<<endl;
        //dout << " calcRecombinationRate(active_node(0)) = "<< calcRecombinationRate(active_node(0)) <<" calcRecombinationRate(active_node(1)) = "<<calcRecombinationRate(active_node(1)) <<endl;
        if (states_[i] == 2) {
            dout << "Active node["<<i<<"] last updated at " << active_node(i)->last_update() <<", current_base = "<< this->current_base()<<" differ "<< ( this->current_base() - active_node(i)->last_update() ) <<endl;
            dout << "recomb_opportunity increase " << ( this->current_base() - active_node(i)->last_update() ) <<"*"<< opportunity_y <<"="<<( this->current_base() - active_node(i)->last_update() ) * opportunity_y<<endl;
            // node i is tracing a non-local branch; opportunities for recombination
            recomb_opportunity += ( this->current_base() - active_node(i)->last_update() ) * opportunity_y;
            //recomb_opportunity += ( this->next_base() - active_node(i)->last_update() ) * opportunity_y;
            }
        if (states_[i] == 1) {
            // node i is tracing out a new branch; opportunities for coalescences and migration
            coal_opportunity += ti.numberOfContemporaries( active_node(i)->population() ) * opportunity_y;
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
                                    ti.start_height(), 
                                    ti.start_height() + opportunity_y, 
                                    coal_opportunity, 
                                    (tmp_event_.isCoalescence() || tmp_event_.isPwCoalescence()) ? EVENT : NOEVENT );
            } 
        else if (states_[0] == 1){
            this->record_Coalevent(active_node(0)->population(), 
                                    ti.start_height(), 
                                    ti.start_height() + opportunity_y, 
                                    coal_opportunity, 
                                    (tmp_event_.isCoalescence() || tmp_event_.isPwCoalescence()) ? EVENT : NOEVENT ); 
            } 
        else if (states_[1] == 1){
            this->record_Coalevent(active_node(1)->population(), 
                                    ti.start_height(), 
                                    ti.start_height() + opportunity_y, 
                                    coal_opportunity, 
                                    (tmp_event_.isCoalescence() || tmp_event_.isPwCoalescence()) ? EVENT : NOEVENT );
            }
        }

    if (migr_opportunity > 0 && this->model().population_number()>1) {
        if ( tmp_event_.isMigration() ){
            this->record_Migrevent(tmp_event_.node()->population(), 
                                    tmp_event_.mig_pop(), 
                                    ti.start_height(), 
                                    ti.start_height() + opportunity_y, 
                                    migr_opportunity, EVENT );    
            } 
        else {
            this->record_Migrevent(active_node(0)->population(),    
                                    size_t(-1),           
                                    ti.start_height(), 
                                    ti.start_height() + opportunity_y, 
                                    migr_opportunity, NOEVENT );
            }
        }
    
    if (recomb_opportunity > 0) {
        // do not store the population, as the recomb_opportunity is also calculated 
        // over the full tree rather than per-population
        if (tmp_event_.isRecombination()) {
            this->record_Recombevent( 0, // active_node( tmp_event_.getActiveNode() )->population(),
                                     ti.start_height(), 
                                     ti.start_height() + opportunity_y, 
                                     recomb_opportunity,
                                     EVENT );
        } else {
            this->record_Recombevent( 0,
                                      ti.start_height(), 
                                      ti.start_height() + opportunity_y, 
                                      recomb_opportunity,
                                      NOEVENT );
        }
    }
             
        //if (states_[0] == 2){
            //this->record_Recombevent(active_node(0)->population(), 
                                        ////ti.start_height(), 
                                        ////ti.start_height() + opportunity_y, 
                                        //recomb_opportunity, 
                                        //tmp_event_.isRecombination() ? EVENT : NOEVENT );
            //} 
        //else if (states_[1] == 2){
            //this->record_Recombevent(active_node(1)->population(), 
                                        ////ti.start_height(), 
                                        ////ti.start_height() + opportunity_y, 
                                        //recomb_opportunity, 
                                        //tmp_event_.isRecombination() ? EVENT : NOEVENT );
            //}
        //}
    
    //assert(active_node(0)->population() == tmp_event_.node()->population());
  
    //// CHECKING OPPORTUNITIES
    //double act0coal_rate =  1 / ( 2.0 * this->model().population_size(active_node(0)->population()) );
    ////  CHECK ACTIVE NODE 0 AND ACTIVE NODE 1 SHOUDL BE IN THE SAME POPULATION
    //assert(active_node(0)->population() == active_node(1)->population() );
    
    //double recomb_rate = this->model().recombination_rate();
    
    //double actmig_rate[2];
    //actmig_rate[0] = 0.0;
    //actmig_rate[1] = 0.0;
    //for (int i=0; i<2; i++) {
        //if (states_[i] == 1) {
            //// node i is tracing out a new branch; opportunities for coalescences and migration
            //actmig_rate[i] = model().total_migration_rate(active_node(i)->population());
            //}
        //}
    
    //double opportunity_rate = (coal_opportunity / opportunity_y * act0coal_rate + recomb_opportunity / opportunity_y * recomb_rate + actmig_rate[0] + actmig_rate[1]);
    //if (abs(rates_[0] - opportunity_rate) > 0.00000000000001){
    ////if (abs(rates_[0] - opportunity_rate) != 0){
        //cout << "difference is "<<rates_[0]-opportunity_rate<<endl;
        //cout<<"rates_[0]= "<<rates_[0]<<endl;
        //cout << "opportunity_rate = " << opportunity_rate<<endl;
        //cout<<"                         recomb_rate = " << recomb_opportunity * recomb_rate/ opportunity_y<<endl;
        //cout<<"                         coal_rate = " << coal_opportunity * act0coal_rate / opportunity_y<<endl;
        //cout<<"                         mig_rate0 = " << model().total_migration_rate(active_node(0)->population())<<endl;
        //cout<<"                         mig_rate1 = " << model().total_migration_rate(active_node(1)->population())<<endl;
        ////cout<<" calcCoalescenceRate(active_node(0)->population(), ti);  = "<<calcCoalescenceRate(active_node(0)->population(), ti) <<endl;
        ////cout<<" calcCoalescenceRate(active_node(1)->population(), ti);  = "<<calcCoalescenceRate(active_node(1)->population(), ti) <<endl;
        ////cout<<" calcPwCoalescenceRate(active_node(1)->population()) = "<<calcPwCoalescenceRate(active_node(1)->population())<<endl; 
        //}
        
    //if (model().total_migration_rate(active_node(0)->population()) != model().total_migration_rate(active_node(0)->population())){
        //cout<<"model().total_migration_rate(active_node(0)->population()) != model().total_migration_rate(active_node(0)->population())"<<endl;
        //exit(1);
        //}
    return;
    }


/*! \brief Record Coalescent events
* @ingroup group_count_coal
*/
void ForestState::record_Coalevent(
                  size_t pop_i,
                  double start_time, 
                  double end_time, 
                  double opportunity, 
                  eventCode event_code) {
    //cout<<"current_time_index = " << this->writable_model()->current_time_idx_<<endl;
    Starevent * new_event = new Starevent( pop_i,
                          start_time,
                          end_time, 
                          opportunity,
                          event_code);	
	new_event->set_change_time_i ( this->writable_model()->current_time_idx_ ) ;
    new_event->set_base ( this->current_base() );
    assert(new_event->print_event("Coalsecent"));
    this->CoaleventContainer[this->writable_model()->current_time_idx_].push_back(new_event);
    
    }


/*! \brief Record Recombination events
* @ingroup group_count_coal
*/
void ForestState::record_Recombevent(size_t pop_i,
                          double start_time, 
                          double end_time, 
                          double opportunity, 
                          eventCode event_code){
    Starevent* new_event = new Starevent( pop_i,
                          start_time,
                          end_time, 
                          opportunity,
                          event_code);	
	new_event->set_change_time_i ( this->writable_model()->current_time_idx_ ) ;
    new_event->set_base ( this->current_base() );
    //if (event_code==EVENT){cout << "recomb+1" <<endl;}
    //cout<<"opportunity: "<<opportunity<<endl;
    assert(new_event->print_event("Recombination"));
    this->RecombeventContainer[this->writable_model()->current_time_idx_].push_back(new_event);
    }
    
    
/*! \brief Record Migration events
* @ingroup group_count_coal
*/
void ForestState::record_Migrevent(size_t pop_i,
                          size_t mig_pop,
                          double start_time, 
                          double end_time, 
                          double opportunity, 
                          eventCode event_code) {
    Migrevent* new_event = new Migrevent( pop_i,
                          mig_pop,
                          start_time,
                          end_time, 
                          opportunity,
                          event_code);	                                          
    new_event->set_change_time_i ( this->writable_model()->current_time_idx_ ) ;
    new_event->set_base ( this->current_base() );
	this->MigreventContainer[this->writable_model()->current_time_idx_].push_back(new_event);
    }    


/*! Clear coalescent and recombination events recorded between two states.*/
void ForestState::clear_CoaleventContainer(){ 
    for (size_t time_i = 0; time_i < this->CoaleventContainer.size(); time_i++ ){
        if ( this->CoaleventContainer[time_i].size() >0 ){
            cout << " This should not happen!!!"<<endl;
            cout<<"this->CoaleventContainer[time_i].size()"<<this->CoaleventContainer[time_i].size()<<endl;
        }
        for (size_t i=0; i < this->CoaleventContainer[time_i].size(); i++){
            CoaleventContainer[time_i][i]->pointer_counter_ --;
            if (CoaleventContainer[time_i][i]->pointer_counter_ == 0){
                cout<<"this should have been removed before"<<endl;
                delete CoaleventContainer[time_i][i];
                }
            }
        //this->CoaleventContainer[time_i].clear(); 
        }
	//this->CoaleventContainer.clear();
    }


/*! Clear recombination events recorded between two states.*/
void ForestState::clear_RecombeventContainer(){ 
    for (size_t time_i = 0; time_i < this->RecombeventContainer.size(); time_i++){
        if ( this->RecombeventContainer[time_i].size() >0 ){
            cout << " This should not happen!!!"<<endl;
            cout<<"this->RecombeventContainer[time_i].size()"<<this->RecombeventContainer[time_i].size()<<endl;
        }
        for (size_t i=0; i < this->RecombeventContainer[time_i].size(); i++){
            RecombeventContainer[time_i][i]->pointer_counter_ --;
            if (RecombeventContainer[time_i][i]->pointer_counter_ == 0){
                cout<<"this should have been removed before"<<endl;
                delete RecombeventContainer[time_i][i];
                }
            }
        //this->RecombeventContainer.clear();
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
        //this->MigreventContainer.clear();
        }
    }
    

void ForestState::include_haplotypes_at_tips(vector <bool> haplotypes_at_tips){
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


inline valarray<double> ForestState::cal_marginal_likelihood_infinite(Node * node){// Genealogy branch lengths are in number of generations, the mutation rate is unit of per site per generation, often in the magnitute of 10 to the power of negative 8.
	double mutation_rate = this->model().mutation_rate();
	valarray<double> marginal_likelihood(2);
	dout << "subtree at " << node << " first child is " << node->first_child() <<" second child is " <<  node->second_child()<<endl;
	if ( node->first_child() == NULL && ((node->label())>0) ){
        marginal_likelihood[1] = node->mutation_state() ? 1.0 : 0.0;
        marginal_likelihood[0] = node->mutation_state() ? 0.0 : 1.0;	
		dout << "Marginal probability at " << node->label() << " is " << marginal_likelihood[0]<<"," << marginal_likelihood[1]<<endl;
		return marginal_likelihood;
        }
	else{ // this is an interior node, but need to check if it is real, i.e. any of its children is a local
		Node *left = trackLocalNode(node->first_child());
		double t1=node->height()- left->height();
        double ut1 = 1 - exp(-t1 * mutation_rate); // assume infinite site
		assert(ut1>=0 && ut1<=1);
		valarray<double> y = cal_marginal_likelihood_infinite(left);
		Node *right = trackLocalNode(node->second_child());
		double t2=node->height()- right->height();
		double ut2 = 1 - exp(-t2 * mutation_rate); // assume infinite site
        assert(ut2>=0 && ut2<=1);
		valarray<double> z = cal_marginal_likelihood_infinite(right);		
		marginal_likelihood[0] = (y[0]*ut1 + y[1]*(1-ut1)) * (z[0]*ut2 + z[1]*(1-ut2)) ;
		marginal_likelihood[1] = (y[1]*ut1 + y[0]*(1-ut1)) * (z[1]*ut2 + z[0]*(1-ut2)) ;
		dout << "Marginal probability at " << node->label() << " is " << marginal_likelihood[0]<<"," << marginal_likelihood[1]<<endl;
		
		//dout << "node is " << node<<", t1=" << t1<<", t2=" << t2<<endl;
		//dout << "prob is " << ", ut1=" << ut1<<", ut2=" << ut2<<endl;
		//dout << "prob is " << ", y[0]=" << y[0]<<", y[1]=" << y[1]<<endl;
		//dout << "prob is " << ", z[0]=" << z[0]<<", z[1]=" << z[1]<<endl;
		//dout << "marginal_likelihood[0] =(" << y[0]*ut1 <<"+" <<  y[1]*(1-ut1) <<")*(" <<  z[0]*ut2 <<"+" <<  z[1]*(1-ut2)<<")" << endl ;
        //dout << ", marginal_likelihood  = " << marginal_likelihood[0]<<", " <<  marginal_likelihood[1]<<endl;
		return marginal_likelihood;
        }				
    }


	
/*! 
 * \brief Calculate the likelihood of the genealogy at data site i, 
 *  If there is no data given at the site i, return likelihood as 1. Since all particles at this site are equally probable 
 * @ingroup group_pf_resample
 */
double ForestState::calculate_likelihood( bool withdata ) {
	if (withdata){
		//double mutation_rate = this->model().mutation_rate();
		dout << "calculate_likelihood function, root is " <<  this->local_root()<<endl;
		valarray<double> marginal_likelihood = cal_marginal_likelihood_infinite(this->local_root());
		dout << "marginal likelihood is " << marginal_likelihood[0]<< "," << marginal_likelihood[1] <<endl;
		double prior[2] = {0.5,0.5};
		double likelihood = marginal_likelihood[0]*prior[0] + marginal_likelihood[1]*prior[1];
		dout << "likelihood is " << likelihood<<endl;
		return likelihood;
        }
	else{ return 1;}
    }

