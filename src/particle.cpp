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
//#include"usage.hpp"

/*! \brief Create a newly copied ForestState 
    @ingroup group_pf_resample 
*/
ForestState::ForestState( ForestState * copied_state )
            :Forest( copied_state ) {
	this->init( copied_state->weight(), 
                copied_state->site_where_weight_was_updated(), 
                copied_state); // initialize members of ForestState
	dout << "current particle's weight is " << this->weight()<<endl;
	copied_state->pointer_counter++;
    }


/*! \brief Initialize a new ForestState 
 * @ingroup group_pf_init
 * */
ForestState::ForestState( Model* model, RandomGenerator* random_generator )
            :Forest( model,random_generator ) {/*! Initialize base of a new ForestState */
	this->init(); // initialize members of ForestState
	this->buildInitialTree();
    assert( this->print_Coalevent() );
    assert( this->print_Migrevent() );
    }


/*! \brief Destructor of ForestState
 * Recursively remove all the previous states, if the pointer counter is zero
 */
ForestState::~ForestState(){
	dout << "State between " << this->current_base() << " and " 
                             << this->next_base() 
                             << " is about to be removed" << endl;	
	assert( this->pointer_counter == 0 );
	this->clear_CoaleventContainer();
	this->clear_RecombeventContainer();
    this->clear_MigreventContainer();
	//Remove any of the previous states, if the counter is equal to zero.
	ForestState* prior_state = this->previous_state;
	if ( prior_state != NULL ) {
		prior_state->pointer_counter--;
		if ( prior_state->pointer_counter == 0 ){ delete prior_state; }
        } 	
    delete_forest_counter++;
    //cout<<"Forest state destructor is called " << delete_forest_counter<<endl;    
	dout << "A Foreststate is deleted" << endl;
    }


/*! \brief Initialize members of an ForestState 
*/
void ForestState::init(double weight, double site , ForestState* previous_state){
	this->setParticleWeight(weight);
	this->setSiteWhereWeightWasUpdated(site);
	this->previous_state=previous_state;
	this->pointer_counter=0;
	this->CoaleventContainer.clear();
    new_forest_counter++;
    }

//ForestState::ForestState(){
	//init();
//}



void ForestState::record_all_event(TimeInterval const &ti){
    double coal_opportunity = 0.0;
    double recomb_opportunity = 0.0;
    double migr_opportunity = 0.0;

    //opportunity_y is the branch length of the conterporaries within the interval
    // if there is no events, then take the full length of this time interval
    //             otherwise, take the distance between the event and the bottom of this interval
    double opportunity_y = this->tmp_event_.isNoEvent() ? ti.length() : (this->tmp_event_.time() - ti.start_height());

    for (int i=0; i<2; i++) {
        if (states_[i] == 2) {
            // node i is tracing a non-local branch; opportunities for recombination
            recomb_opportunity += ( this->current_base() - active_node(i)->last_update() ) * opportunity_y;
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
            this->record_Coalevent(active_node(0)->population(), ti.start_height(), ti.start_height() + opportunity_y, coal_opportunity, (tmp_event_.isCoalescence() || tmp_event_.isPwCoalescence()) ? EVENT : NOEVENT );
            } 
        else if (states_[0] == 1){
            this->record_Coalevent(active_node(0)->population(), ti.start_height(), ti.start_height() + opportunity_y, coal_opportunity, (tmp_event_.isCoalescence() || tmp_event_.isPwCoalescence()) ? EVENT : NOEVENT ); 
            } 
        else if (states_[1] == 1){
            this->record_Coalevent(active_node(1)->population(), ti.start_height(), ti.start_height() + opportunity_y, coal_opportunity, (tmp_event_.isCoalescence() || tmp_event_.isPwCoalescence()) ? EVENT : NOEVENT );
            }
        }

    if (migr_opportunity > 0 && this->model().population_number()>1) {
        if ( tmp_event_.isMigration() ){
            this->record_Migrevent(tmp_event_.node()->population(), tmp_event_.mig_pop(), ti.start_height(), ti.start_height() + opportunity_y, migr_opportunity, EVENT );    
            } 
        else {
            this->record_Migrevent(active_node(0)->population(),    size_t(-1),           ti.start_height(), ti.start_height() + opportunity_y, migr_opportunity, NOEVENT );
            }
        }
    
    if (recomb_opportunity > 0) {
        if (states_[0] == 2){
            this->record_Recombevent(active_node(0)->population(), ti.start_height(), ti.start_height() + opportunity_y, recomb_opportunity, tmp_event_.isRecombination() ? EVENT : NOEVENT );
            } 
        else if (states_[1] == 2){
            this->record_Recombevent(active_node(1)->population(), ti.start_height(), ti.start_height() + opportunity_y, recomb_opportunity, tmp_event_.isRecombination() ? EVENT : NOEVENT );
            }
        }
    
    //assert(active_node(0)->population() == tmp_event_.node()->population());
  
    // CHECKING OPPORTUNITIES
    double act0coal_rate =  1 / ( 2.0 * this->model().population_size(active_node(0)->population()) );
    ////  CHECK ACTIVE NODE 0 AND ACTIVE NODE 1 SHOUDL BE IN THE SAME POPULATION
    //assert(active_node(0)->population() == active_node(1)->population() );
    
    double recomb_rate = this->model().recombination_rate();
    
    double actmig_rate[2];
    actmig_rate[0] = 0.0;
    actmig_rate[1] = 0.0;
    for (int i=0; i<2; i++) {
        if (states_[i] == 1) {
            // node i is tracing out a new branch; opportunities for coalescences and migration
            actmig_rate[i] = model().total_migration_rate(active_node(i)->population());
            }
        }
    
    double opportunity_rate = (coal_opportunity / opportunity_y * act0coal_rate + recomb_opportunity / opportunity_y * recomb_rate + actmig_rate[0] + actmig_rate[1]);
    if (abs(rates_[0] - opportunity_rate) > 0.00000000000001){        
    //if (abs(rates_[0] - opportunity_rate) != 0){        
        cout << "difference is "<<rates_[0]-opportunity_rate<<endl;
        cout<<"rates_[0]= "<<rates_[0]<<endl;
        cout << "opportunity_rate = " << opportunity_rate<<endl;
        cout<<"                         recomb_rate = " << recomb_opportunity * recomb_rate/ opportunity_y<<endl;
        cout<<"                         coal_rate = " << coal_opportunity * act0coal_rate / opportunity_y<<endl;
        cout<<"                         mig_rate0 = " << model().total_migration_rate(active_node(0)->population())<<endl;
        cout<<"                         mig_rate1 = " << model().total_migration_rate(active_node(1)->population())<<endl;
        //cout<<" calcCoalescenceRate(active_node(0)->population(), ti);  = "<<calcCoalescenceRate(active_node(0)->population(), ti) <<endl;
        //cout<<" calcCoalescenceRate(active_node(1)->population(), ti);  = "<<calcCoalescenceRate(active_node(1)->population(), ti) <<endl;
        //cout<<" calcPwCoalescenceRate(active_node(1)->population()) = "<<calcPwCoalescenceRate(active_node(1)->population())<<endl; 
        }
        
    if (model().total_migration_rate(active_node(0)->population()) != model().total_migration_rate(active_node(0)->population())){
        cout<<"model().total_migration_rate(active_node(0)->population()) != model().total_migration_rate(active_node(0)->population())"<<endl;
        exit(1);
        }
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
    Coalevent* new_event = new Coalevent( pop_i,
                                          start_time,
                                          end_time, 
                                          opportunity,
			                              event_code);	
	this->CoaleventContainer.push_back(new_event);
    }


/*! \brief Record Recombination events
* @ingroup group_count_coal
*/
void ForestState::record_Recombevent(size_t pop_i,
                          double start_time, 
                          double end_time, 
                          double opportunity, 
                          eventCode event_code){
    Recombevent* new_event = new Recombevent( pop_i,
                                          start_time,
                                          end_time, 
                                          opportunity,
			                              event_code);	
	this->RecombeventContainer.push_back(new_event);
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
	this->MigreventContainer.push_back(new_event);
    }    


/*! Clear coalescent and recombination events recorded between two states.*/
void ForestState::clear_CoaleventContainer(){ 
	for (size_t i=0; i < this->CoaleventContainer.size(); i++){
		delete CoaleventContainer[i];
    	}
	this->CoaleventContainer.clear();
    }


/*! Clear recombination events recorded between two states.*/
void ForestState::clear_RecombeventContainer(){ 
	for (size_t i=0; i < this->RecombeventContainer.size(); i++){
		delete RecombeventContainer[i];
    	}
	this->RecombeventContainer.clear();
    }
    
    
/*! Clear migration events recorded between two states.*/
void ForestState::clear_MigreventContainer(){ 
	for (size_t i=0; i < this->MigreventContainer.size(); i++){
		delete MigreventContainer[i];
    	}
	this->MigreventContainer.clear();
    }
    

void ForestState::include_haplotypes_at_tips(vector <bool> haplotypes_at_tips){
	for (size_t j=0; j < haplotypes_at_tips.size();j++){		
		this->nodes()->at(j)->set_mutation_state(haplotypes_at_tips[j]);
	    }
    }


/*! Calculate the marginal likelihood of a node recursively, let X denote the state of Node *, and Y denote the state of the 
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
valarray<double> ForestState::cal_marginal_likelihood(Node * node){// Genealogy branch lengths are in number of generations, the mutation rate is unit of per site per generation, often in the magnitute of 10 to the power of negative 8.
	double mutation_rate = this->model().mutation_rate();
	valarray<double> marginal_likelihood(2);
	dout << "subtree at " << node << " first child is " << node->first_child() <<" second child is " <<  node->second_child()<<endl;
	if ( node->first_child() == NULL && ((node->label())>0) ){
        marginal_likelihood[1] = node->mutation_state() ? 1.0 : 0.0;
        marginal_likelihood[0] = node->mutation_state() ? 0.0 : 1.0;	
		//if (node->mutation_state()){
			//marginal_likelihood[1]=1.0;
			//marginal_likelihood[0]=0.0;	
		//}
		//else{
			//marginal_likelihood[1]=0.0;
			//marginal_likelihood[0]=1.0;		
		//}
		dout << "Marginal probability at " << node->label() << " is " << marginal_likelihood[0]<<"," << marginal_likelihood[1]<<endl;
		return marginal_likelihood;
        }
	else{ // this is an interior node, but need to check if it is real, i.e. any of its children is a local
		Node *left = trackLocalNode(node->first_child());
		double t1=node->height()- left->height();
		double ut1 = 0.5 + 0.5*exp(-t1*2*mutation_rate); // let ut1 be the probability that either end of the branch to the first child carries the same state
        //double ut1 = 1 - exp(-t1 * mutation_rate); // assume infinite site
		assert(ut1>=0 && ut1<=1);
		valarray<double> y = cal_marginal_likelihood(left);
		Node *right = trackLocalNode(node->second_child());
		double t2=node->height()- right->height();
		double ut2 = 0.5 + 0.5*exp(-t2*2*mutation_rate); // let ut2 be the probability that either end of the branch to the second child carries the same state
		//double ut2 = 1 - exp(-t2 * mutation_rate); // assume infinite site
        assert(ut2>=0 && ut2<=1);
		valarray<double> z = cal_marginal_likelihood(right);		
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
double ForestState::calculate_likelihood(bool withdata) {
	if (withdata){
		//double mutation_rate = this->model().mutation_rate();
		dout << "calculate_likelihood function, root is " <<  this->local_root()<<endl;
		valarray<double> marginal_likelihood=cal_marginal_likelihood(this->local_root());
		dout << "marginal likelihood is " << marginal_likelihood[0]<< "," << marginal_likelihood[1] <<endl;
		double prior[2] = {0.5,0.5};
		double likelihood = marginal_likelihood[0]*prior[0] + marginal_likelihood[1]*prior[1];
		dout << "likelihood is " << likelihood<<endl;
		return likelihood;
        }
	else{ return 1;}
    }


//ParticleContainer::ParticleContainer(){};


/*! \brief Particle filtering Initialization 
 * Create particle initial states in the simulation
 * 
 * \ingroup group_pf_init
 */ 
ParticleContainer::ParticleContainer(
                    Model* model, 
                    RandomGenerator* rg, 
                    size_t Num_of_states, 
                    vector <bool> data_for_init_states, 
                    bool withdata,
                    double initial_position ){
    this->set_random_generator(rg);
    this->set_ESS(0);
    this->set_current_printing_base(0);    
	dout << " --------------------   Particle Initial States   --------------------" << std::endl;	
	for (size_t i=0; i<Num_of_states ;i++){
		ForestState *  new_state = new ForestState(model,rg);  // create a new state, using scrm; scrm always starts at 0.
        new_state->setSiteWhereWeightWasUpdated( initial_position );
		this->push(new_state, 1.0/Num_of_states );
	    }	    
    }


void ParticleContainer::ESS_resampling(valarray<double> weight_cum_sum, valarray<int> &sample_count, int mutation_at, double ESSthreshold, int num_state){    
    //runningtime->stopwatch_start();    
    dout << "ESS is " <<  this->ESS() <<", number of particle is " <<  num_state <<endl;
    if (this->ESS() < ESSthreshold){ // resample if the effective sample size is small
        dout << " ### PROGRESS: ESS resampling" << endl;
        this->systemetic_resampling( weight_cum_sum, sample_count, num_state);
        //this->trivial_resampling(pfARG_para.N,sample_count);             
        this->resample(sample_count);  
        }    
        //runningtime->stopwatch_end(1);
    }


/*! \brief Particle filtering Resampling step
 * 
 *  Update particles according to the sample counts
 * 
 * \ingroup group_resample
 */ 
void ParticleContainer::resample(valarray<int> & sample_count){	
	dout << " ****************************** Start making list of new states ****************************** " << std::endl;
	dout << " will make total of " << sample_count.sum()<<" particle states" << endl;
	for (size_t old_state_index = 0; old_state_index < sample_count.size(); old_state_index++){
        if ( (sample_count[old_state_index] > 0)){
            for (int ii=1; ii <= sample_count[old_state_index]; ii++){
			  	dout << "       Make a copy of the " << old_state_index<<"th particle" << endl;
                // create new copy of the resampled particle        
				ForestState * new_copy_state= new ForestState(this->particles[old_state_index]);
				this->push(new_copy_state);
                }
            } 
        else {
			delete this->particles[old_state_index];			
			this->particles[old_state_index]=NULL;
            }
        }
    
    this->shifting(sample_count.sum());
    
    dout<<"sample_count.sum() = "<<sample_count.sum()<<endl;
    assert(sample_count.sum() == int(this->particles.size()));
    dout<<this->particles.size()<<endl;
	dout << " ****************************** End of making list of new particles ****************************** " << std::endl;
	assert(this->check_state_orders());
}


void ParticleContainer::shifting(int number_of_particles){	
    /*! \todo SHIFTING */
    // SHIFTING AND RESIZING THE PARTICLE ARRAY
    // SHIFTING, this is not pretty, there must be some other ways to do this ...
    for (int i = 0; i < number_of_particles; i++){
        this->particles[i] = this->particles[i+number_of_particles];
        }
    this -> particles.resize(number_of_particles);
    this -> particles.shrink_to_fit();
    // SHIFTING AND RESIZING THE PARTICLE ARRAY FINISHED    
    }


/*! \brief Update particle weight according to the haplotype data
 *	@ingroup group_pf_update
 */
void ParticleContainer::update_state_weights_at_A_single_site(
    double mutation_at,
    double mutation_rate, 
    bool withdata,
    vector <bool> haplotypes_at_tips){
			
	// now update the weights of all particles, by calculating the likelihood of the data over the previous segment	
	for (size_t particle_i=0; particle_i < this->particles.size(); particle_i++){
		this->particles[particle_i]->include_haplotypes_at_tips(haplotypes_at_tips);

		double likelihood_of_haplotypes_at_tips = this->particles[particle_i]->calculate_likelihood(withdata);
        dout << "updated weight =" << this->particles[particle_i]->weight()  << "*" <<  likelihood_of_haplotypes_at_tips <<endl;

        this->particles[particle_i]->setParticleWeight( this->particles[particle_i]->weight() * likelihood_of_haplotypes_at_tips);
		dout << "particle " <<  particle_i<<" done" << endl;
        }
    
    this->normalize_probability(); // It seems to converge slower if it is not normalized ...
	dout << endl;
    }




ParticleContainer::~ParticleContainer(){
    // The following message may show up very often if it is passed by reference ...
    //dout << "ParticleContainer destructor is called" << endl;
}


void ParticleContainer::clear(){
	dout << "ParticleContainer clear() is called" << endl;
	for (size_t i = 0; i < this->particles.size(); i++){
		if (this->particles[i]!=NULL){
			delete this->particles[i];
			this->particles[i]=NULL;
            }
        }
	particles.clear();
	dout << "Particles are deleted" << endl;
    }
 

void ParticleContainer::push(ForestState* state, double weight){
	state->setParticleWeight(weight);
	this->particles.push_back(state);              
    }




/*! 
 * @ingroup group_pf_resample
 * @ingroup group_pf_update
 * \brief Calculate the effective sample size, and update the cumulative weight of the particles
 */
void ParticleContainer::update_cum_sum_array_find_ESS(std::valarray<double> & weight_cum_sum){
	double wi_sum=0;
	double wi_sq_sum=0;
	double Num_of_states = this->particles.size();
	weight_cum_sum=0; //Reinitialize the cum sum array 
	
	for (size_t i=0; i<Num_of_states ;i++){
		//update the cum sum array
		double w_i=this->particles[i]->weight();
		weight_cum_sum[i+1]=weight_cum_sum[i]+w_i;
		wi_sum = wi_sum + w_i;
		//wi_sq_sum = wi_sq_sum + pow(this->particles[i]->weight,2);
		wi_sq_sum = wi_sq_sum + w_i * w_i;
	}

	//check for the cum weight
	dout << "### particle weights ";
    for (size_t i=0;i<Num_of_states;i++){
        dout << this->particles[i]->weight()<<"  ";
    }
    dout << std::endl<<std::endl;
	
    dout << "### updated cum sum of particle weight ";
    for (size_t i=0;i<weight_cum_sum.size();i++){
        dout << weight_cum_sum[i]<<"  ";
    }dout << std::endl;

	this->set_ESS(wi_sum * wi_sum / wi_sq_sum);	
}



/*! 
 * @ingroup group_pf_update
 * \brief Update the current state to the next state, at the given site, update all particles to it's latest genealogy state.  Also include the likelihood for no mutations.
 */
void ParticleContainer::extend_ARGs(double mutation_at, double mutation_rate, bool withdata){

 	for (size_t particle_i=0;particle_i < this->particles.size(); particle_i++){
        dout << "We are updating particle " << particle_i << endl;
        /*! 
         * For each particle, extend the current path until the the site such that the next genealogy change is beyond the mutation
         * Invariant: the likelihood is correct up to 'updated_to'
         */
        double updated_to = this->particles[particle_i]->site_where_weight_was_updated();
        while (updated_to < mutation_at) {
            
            dout << "  Now at " << updated_to << " and extending to " << mutation_at << endl;            
            /*!
             * First, update the likelihood up to either mutation_at or the end of this state
             */
            double update_to = min( mutation_at, this->particles[particle_i]->next_base() );
            double length_of_local_tree = this->particles[particle_i]->local_tree_length(); // in generations
            double likelihood_of_segment;
            if (withdata) {
                likelihood_of_segment = exp( -mutation_rate * length_of_local_tree * (update_to - updated_to) ); // assume infinite site model
                //likelihood_of_segment = pow(0.5 - 0.5*exp(-length_of_local_tree*2*mutation_rate ), (update_to - updated_to)); // Finite site model
                dout << " Likelihood of no mutations in segment of length" << (update_to - updated_to) << " is " << likelihood_of_segment << endl;
                } 
            else {
                likelihood_of_segment = 1;
                dout << " no data" << endl;
                }
            this->particles[particle_i]->setParticleWeight( this->particles[particle_i]->weight() * likelihood_of_segment);
            updated_to = update_to;  // rescues the invariant
            
            /*!
             * Next, if we haven't reached mutation_at now, add a new state and iterate
             */
            if (updated_to < mutation_at) {
                ForestState * median_state = new ForestState(this->particles[particle_i]);
                this->particles[particle_i] = median_state;
                this->particles[particle_i]->sampleNextGenealogy();
                }

            assert(this->particles[particle_i]->print_Recombevent());
            assert(this->particles[particle_i]->print_Coalevent());
            assert(this->particles[particle_i]->print_Migrevent());
            }
        
        assert (updated_to == mutation_at);        
        this->particles[particle_i]->setSiteWhereWeightWasUpdated( mutation_at );
        }
    
    /*! normalize the probability upon until the mutation */
    //this->normalize_probability(); // This normalization doesn't seem to do much ...
    }


void ParticleContainer::normalize_probability(){
    double total_probability = 0;
    for ( size_t particle_i = 0;particle_i < this->particles.size(); particle_i++ ){
        total_probability += this->particles[particle_i]->weight();
        }
    for ( size_t particle_i = 0; particle_i < this->particles.size(); particle_i++ ){
        this->particles[particle_i]->setParticleWeight( this->particles[particle_i]->weight() / total_probability);
        }
    }


void ParticleContainer::update_state_to_data(Vcf * VCFfile, Model * model, valarray<double> & weight_cum_sum){
    dout <<  " ******************** Update the weight of the particles  ********** " <<endl;
    dout << " ### PROGRESS: update weight at " << VCFfile->site()<<endl;
    double mutation_at = VCFfile->site();
    bool withdata = VCFfile->withdata();
    double mutation_rate = model->mutation_rate();
    
    /*!
     * \verbatim
                 
                                  previous VCFfile->site()
                              p0  mut0  p1  p2  p3         p4  mut1
                              .   .     .   .   .          .   . 
                              .   .     .   .   .          .   VCFfile->site()
                              .   .         .   .          .   .
                                  .     1       .              .
                                  .     x---o   .          4   .
                              0   .     |   |   .          x-------o
                              x---------o   |   .          |   .
                              |   .         |              |   .
                              |   .         x---o          |   .
         x--------o           |   .         2   |          |   .
                  |           |   .             x----------o   .
                  |           |   .             3              .
                  x-----------o   .                            .
                                  .                            .
                                  .                            .
      
      \endverbatim
     * At the start of calling this function, the tail ForestState is state 0.
     * Let p_i denotes the position of the start of State i.
     * mut0 was the previous mutation position
     * mut1 is the current mutation position
     * 
     * At the start of this function, likelihood of this particle is updated upto mut0
     * When extending the ARG, we first update the likelihood from mut0 to p1 for 
     * not observing any mutation on this segment of sequence. The new state occurs at p1,
     * until p4. The probability of not observing mutation is updated until mut1.
     * 
     * Call function update_state_weights_at_A_single_site(), update the particle weight for
     * observe mutation at site mut1. The probability of all particles are updated until mut1.
     * 
     * Calculate the cumulated likelihood of all particles and the effective sample size (ESS).
     * 
     * 
     */
    
    //Extend ARGs and update weight for not seeing mutations along the equences
    this->extend_ARGs(mutation_at, mutation_rate, withdata);
    
    //Update weight for seeing mutation at the position 
    this->update_state_weights_at_A_single_site(mutation_at, mutation_rate, withdata, VCFfile->vec_of_sample_alt_bool); 
    
    //Update the cumulated probabilities, as well as computing the effective sample size                    
    this->update_cum_sum_array_find_ESS(weight_cum_sum); 
    //runningtime->stopwatch_end(2);
    }


/*! 
 * Remove previous states that are no longer used
 */
void ParticleContainer::clean_old_states(double xstart){
    dout << "remove states before " << xstart <<endl;    
    for (size_t i = 0; i < this->particles.size(); i++){        
        ForestState* current_state = this->particles[i];
        dout << " End particle [" << i << "] lasted from " << current_state->current_base() << " to base " << current_state->next_base() << endl;        
        ForestState* prior_state = current_state->previous_state;        
        
        /*!
         * Check if the current state is pointing a previous state
         *      if yes, check if the previous state should be removed
         *          if yes, remove all the previous state
         *          if no, move on the previous state, and check again
         */ 
         
        while (prior_state!= NULL){
            if ( prior_state->next_base() < xstart ){
                prior_state->pointer_counter--;
                current_state->previous_state = NULL;
                delete prior_state;                
                break;
                }
            current_state = prior_state;
            prior_state = prior_state->previous_state;
            }      
        }
    }


/*! 
 * @ingroup group_naive 
 * \brief Use simple random sampling to resample
 */
void ParticleContainer::trivial_resampling(size_t N, std::valarray<int> & sample_count){
    for (size_t i=0; i<N ;i++){
        size_t index = random_generator()->sampleInt(N); 
        sample_count[index]=sample_count[index]+1;
        }
    }
    

/*! 
 * @ingroup group_systemetic 
 * \brief Use systemetic resampling \cite Doucet2008 to generate sample count for each particle
 */
void ParticleContainer::systemetic_resampling(std::valarray<double> cum_sum, std::valarray<int>& sample_count, int sample_size){
    size_t interval_j = 0;
    size_t sample_i = 0;
    size_t N = sample_size;
    //double u_j = rand() / double(RAND_MAX) / N;
    double u_j = random_generator()->sample() / N;
    double cumsum_normalization = cum_sum[cum_sum.size()-1];

    dout << std::endl<<"systematic sampling procedue" << std::endl;
    for (size_t i=0;i<cum_sum.size();i++){dout <<  (cum_sum[i]/cumsum_normalization )<<"  ";}dout << std::endl;    

    sample_count[sample_i] = 0;
    while (sample_i < N) {
        dout << "Is " <<  u_j<<" in the interval of " << std::setw(10)<< (cum_sum[interval_j]/ cumsum_normalization) << " and " << std::setw(10)<< (cum_sum[interval_j+1]/ cumsum_normalization);
        /* invariants: */    
        assert( (cum_sum[interval_j] / cumsum_normalization) < u_j ); 
        assert( sample_i < N );        
        /* check whether u_j is in the interval [ cum_sum[interval_j], cum_sum[interval_j+1] ) */
        if ( (sample_i == N) || cum_sum[interval_j+1] / cumsum_normalization > u_j ) {
            sample_count[interval_j] += 1;        
            sample_i += 1;
            dout << "  yes, update sample count of particle " << interval_j<<" to " << sample_count[interval_j] <<std::endl;
            u_j += 1.0/double(N);    
            } 
        else {
            dout << "   no, try next interval " << std::endl;        
            //assert( sample_i < N-1 );
            interval_j += 1;
            sample_count[ interval_j ] = 0;        
            }    
        }
    interval_j=interval_j+1;
    for (;interval_j<N;interval_j++){
        sample_count[ interval_j ] = 0;    
        }
    
    dout << "systematic sampling procedue finished with total sample count " << sample_count.sum()<<std::endl<<std::endl;
    for (size_t i=0;i<sample_count.size();i++){dout << sample_count[i]<<"  ";}  dout << std::endl;
    assert(sample_count.sum()==sample_size);
    }


bool ParticleContainer::appendingStuffToFile( double x_end,  pfARG::param pfparam){
    // Record the TMRCA and weight when a heatmap is generated
    if (!pfparam.heat_bool){
        return true;
        }
           /*!
             *  \verbatim
           remove all the state prior to the minimum of
            current_printing_base and 
                previous_backbase     
                .                      backbase                     VCFfile->site()
                .                      .                            .
                .                      .     3                      .
                .                      .     x---o              6   .
                .                  2   .     |   |              x-------o
                .                  x---------o   |              |   .
                .                  |   .         |              |   .
             0  .                  |   .         x---o          |   .
             x---------o           |   .         4   |          |   .
                .      |           |   .             x----------o   .
                .      |           |   .             5              .
                .      x-----------o   .                            .
                .      1               .-------------lag------------.
                .                      .                            .
                Count::update_e_count( .                            ParticleContainer::update_state_to_data(
           x_start = previous_backbase .                              mutation data comes in here
           x_end = backbase            .
                                       .
                                       ParticleContainer::appendingStuffToFile(
                                           x_end = backbase, 
          \endverbatim
         *
         * Likelihood of the particle is updated up until state 6, but because of the lagging we are using
         * report the TMRCA up until state 2
         *           
         */ 
    if (x_end < this->current_printing_base()){
        return true;    
        }
    do {
        this->set_current_printing_base(x_end);

        if (this->current_printing_base() > 0){
            ofstream TmrcaOfstream;
            ofstream WeightOfstream;
            ofstream BLOfstream;
            TmrcaOfstream.open  ( pfparam.TMRCA_NAME.c_str() , ios::out | ios::app | ios::binary); 
            WeightOfstream.open ( pfparam.WEIGHT_NAME.c_str(), ios::out | ios::app | ios::binary); 
            BLOfstream.open     ( pfparam.BL_NAME.c_str()    , ios::out | ios::app | ios::binary);
            
            TmrcaOfstream  << this->current_printing_base();
            WeightOfstream << this->current_printing_base();
            BLOfstream     << this->current_printing_base();
            
            for ( size_t i = 0; i < this->particles.size(); i++){
                ForestState * current_state_ptr = this->particles[i];
                WeightOfstream <<"\t" << current_state_ptr->weight();
                
                while (current_state_ptr->current_base() > this->current_printing_base() && current_state_ptr->previous_state) {           
                    current_state_ptr=current_state_ptr->previous_state;            
                    }
                
                TmrcaOfstream  << "\t" << current_state_ptr->local_root()->height() / (4 * current_state_ptr->model().default_pop_size); // Normalize by 4N0
                BLOfstream     << "\t" << current_state_ptr->local_tree_length()    / (4 * current_state_ptr->model().default_pop_size); // Normalize by 4N0
                current_state_ptr=NULL;
                }
            TmrcaOfstream  << endl;            
            WeightOfstream << endl;
            BLOfstream     << endl;
            
            TmrcaOfstream.close();    
            WeightOfstream.close();                
            BLOfstream.close();
            }
        this->set_current_printing_base(this->current_printing_base() + pfparam.window);
        } while ( this->current_printing_base() < x_end);
    return true;
    }
