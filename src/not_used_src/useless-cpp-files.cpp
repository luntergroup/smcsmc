///*! 
 //* Remove previous states that are no longer used
 //*/
//void ParticleContainer::clean_old_states(double xstart){
    
    //for (size_t i = 0; i < this->particles.size(); i++){
        //ForestState* current_state = this->particles[i];
        //if ( current_state->previous_state == NULL ){
            //continue;
            //}
        //dout << "remove states before " << setw(10) << xstart ;        
        //dout << " End particle [" << i << "] lasted from " << current_state->current_base() << " to base " << current_state->next_base() << endl;        
        
        //ForestState* prior_state = current_state->previous_state;
        
        ///*!
         //* Check if the current state is pointing a previous state
         //*      if yes, check if the previous state should be removed
         //*          if yes, remove all the previous state
         //*          if no, move on the previous state, and check again
         //*/ 
         
        //while (prior_state!= NULL){

            //if ( prior_state->next_base() < xstart ){
                //assert ( current_state->current_base() >= prior_state->next_base() );
                //// reduce the pointer counter before this point!
                //prior_state->pointer_counter--;
                //current_state->previous_state = NULL;
                
                //// once the pointer counter is zero, remove the state
                //if ( prior_state->pointer_counter == (int)0 ){ 
                    //delete prior_state; 
                    //}
                //break;
                //}
            //current_state = prior_state;
            //prior_state = prior_state->previous_state;
            //}      
        //}
    //}

//ForestState::ForestState( ForestState * copied_state )
            //:Forest( copied_state ) {
	//this->init( copied_state->weight(), 
                //copied_state->site_where_weight_was_updated(), 
                //copied_state); // initialize members of ForestState
    //this->setAncestor ( copied_state->ancestor() );
	//dout << "current particle's weight is " << this->weight()<<endl;
	////copied_state->pointer_counter++;
    //}

///*! \brief Initialize members of an ForestState 
//*/
//void ForestState::init(double weight, double site ){
	//this->setParticleWeight(weight);
	//this->setSiteWhereWeightWasUpdated(site);
    //}

//void ForestState::init(double weight, double site , ForestState* previous_state){
	//this->setParticleWeight(weight);
	//this->setSiteWhereWeightWasUpdated(site);
	////this->previous_state=previous_state;
	////this->pointer_counter=(int)0;
	////this->CoaleventContainer.clear();
    ////new_forest_counter++;
    //}

//ForestState::ForestState(){
	//init();
//}

//valarray<double> ForestState::cal_marginal_likelihood_finite(Node * node){// Genealogy branch lengths are in number of generations, the mutation rate is unit of per site per generation, often in the magnitute of 10 to the power of negative 8.
	//double mutation_rate = this->model().mutation_rate();
	//valarray<double> marginal_likelihood(2);
	//dout << "subtree at " << node << " first child is " << node->first_child() <<" second child is " <<  node->second_child()<<endl;
	//if ( node->first_child() == NULL && ((node->label())>0) ){
        //marginal_likelihood[1] = node->mutation_state() ? 1.0 : 0.0;
        //marginal_likelihood[0] = node->mutation_state() ? 0.0 : 1.0;	
		//dout << "Marginal probability at " << node->label() << " is " << marginal_likelihood[0]<<"," << marginal_likelihood[1]<<endl;
		//return marginal_likelihood;
        //}
	//else{ // this is an interior node, but need to check if it is real, i.e. any of its children is a local
		//Node *left = trackLocalNode(node->first_child());
		//double t1=node->height()- left->height();
		//double ut1 = 0.5 + 0.5*exp(-t1*2*mutation_rate); // let ut1 be the probability that either end of the branch to the first child carries the same state
		//assert(ut1>=0 && ut1<=1);
		//valarray<double> y = cal_marginal_likelihood_finite(left);
		//Node *right = trackLocalNode(node->second_child());
		//double t2=node->height()- right->height();
		//double ut2 = 0.5 + 0.5*exp(-t2*2*mutation_rate); // let ut2 be the probability that either end of the branch to the second child carries the same state
        //assert(ut2>=0 && ut2<=1);
		//valarray<double> z = cal_marginal_likelihood_finite(right);		
		//marginal_likelihood[0] = (y[0]*ut1 + y[1]*(1-ut1)) * (z[0]*ut2 + z[1]*(1-ut2)) ;
		//marginal_likelihood[1] = (y[1]*ut1 + y[0]*(1-ut1)) * (z[1]*ut2 + z[0]*(1-ut2)) ;
		//dout << "Marginal probability at " << node->label() << " is " << marginal_likelihood[0]<<"," << marginal_likelihood[1]<<endl;
		
		////dout << "node is " << node<<", t1=" << t1<<", t2=" << t2<<endl;
		////dout << "prob is " << ", ut1=" << ut1<<", ut2=" << ut2<<endl;
		////dout << "prob is " << ", y[0]=" << y[0]<<", y[1]=" << y[1]<<endl;
		////dout << "prob is " << ", z[0]=" << z[0]<<", z[1]=" << z[1]<<endl;
		////dout << "marginal_likelihood[0] =(" << y[0]*ut1 <<"+" <<  y[1]*(1-ut1) <<")*(" <<  z[0]*ut2 <<"+" <<  z[1]*(1-ut2)<<")" << endl ;
        ////dout << ", marginal_likelihood  = " << marginal_likelihood[0]<<", " <<  marginal_likelihood[1]<<endl;
		//return marginal_likelihood;
        //}				
    //}


//Starevent::Starevent(){
    //this->init();
    //}
    
///*! \brief Copy the Starevent from the previous ForestState*/    
//Starevent::Starevent(Starevent* previous_Starevent){ 
    //this->pop_i_        = previous_Starevent->pop_i();
    ////this->start_height_ = previous_Starevent->start_height();
    ////this->end_height_   = previous_Starevent->end_height(); 
    //this->num_event_    = previous_Starevent->num_event();
    //this->opportunity_  = previous_Starevent->opportunity();
    //this->event_state_  = previous_Starevent->event_state();
    
    //this->change_time_i_ = previous_Starevent->change_time_i();
    ////this->counted_ = previous_Starevent->counted();
    //this->base_ = previous_Starevent->base();
    //}
    
///*! \brief Copy the Starevent from the previous ForestState*/    
//Starevent::Starevent(const Starevent & previous_Starevent){ 
    //this->pop_i_        = previous_Starevent.pop_i();
    ////this->start_height_ = previous_Starevent.start_height();
    ////this->end_height_   = previous_Starevent.end_height(); 
    //this->num_event_    = previous_Starevent.num_event();
    //this->opportunity_  = previous_Starevent.opportunity();
    //this->event_state_  = previous_Starevent.event_state();
    
    //this->change_time_i_ = previous_Starevent.change_time_i();
    ////this->counted_ = previous_Starevent.counted();
    //this->base_ = previous_Starevent.base();
    //}

//bool Coalevent::print_event(){
    //cout 
        ////<< setw(10) << this->CoaleventContainer[i]->start_height()  << " to " 
             ////<< setw(10) << this->CoaleventContainer[i]->end_height()    << ", " 
             //<< setw(13) << this->opportunity()   << " opportunity for " 
             //<< setw(2)  << this->num_event()     << " coalescent, ";
        //if ( this->event_state() == NOEVENT ){
            //cout<< " potetial coalsecent";
            //}
        //cout << endl;
    //return true;
    //}


//bool Recombevent::print_event(){
    //cout 
        ////<< setw(10) << this->start_height()  << " to " 
         ////<< setw(10) << this->end_height()    << ", " 
         //<< setw(13) << this->opportunity()   << " opportunity for " 
         //<< setw(2)  << this->num_event()     << " recombination, ";
    //if ( this->event_state() == NOEVENT ){
        //cout<< " potetial recombination";
        //}
    //cout << endl;

    //return true;    
    //}

///*!
 //* \brief Derived class of Starevent, recording the number and the time intervals of Coalescent events between two ForestState 
 //*/
//class Coalevent : public Starevent{
    ////friend class CountModel;
    ////friend class ForestState;
    ////friend class ParticleContainer;
    //public:
        //Coalevent(size_t pop_i, 
                  ////double start_time,
                  ////double end_time, 
                  //double opportunity,
                  //eventCode event_code ) 
                 //: Starevent ( pop_i, 
                    ////start_time, 
                    ////end_time, 
                    //opportunity, 
                    //event_code ){ assert (this->print_event()); };  
        //Coalevent(const Coalevent & previous_Starevent) 
                 //: Starevent (previous_Starevent) { };
        //~Coalevent(){};
    //private:
        //bool print_event();
    //};
    
///*!
 //* \brief Derived class of Starevent, recording the number and the time intervals of Recombination events between two ForestState 
 //*/    
//class Recombevent : public Starevent{
    ////friend class CountModel;
    ////friend class ForestState;
    ////friend class ParticleContainer;
    //public:
        //Recombevent(size_t pop_i, 
                    ////double start_time,
                    ////double end_time, 
                    //double opportunity,
                    //eventCode event_code ) 
                   //: Starevent ( pop_i, 
                      ////start_time, 
                      ////end_time, 
                      //opportunity, 
                      //event_code ){ assert( this->print_event()); };
        //Recombevent(const Recombevent & previous_Starevent) 
                   //: Starevent (previous_Starevent) { };
        //~Recombevent(){};
    //private:
        //bool print_event();
    //};
    

//void CountModel::compute_recomb_rate () {
    //this->resetTime();
    //double scaling_pop_size_N0 = this->population_size();
    //double recomb_opportunity = 0;
    ////double recomb_count = 0;
    //this->recomb_count_ = 0;
    
    //for (size_t time_layer_i = 0; time_layer_i<change_times_.size(); time_layer_i++){
        //for (size_t pop_j = 0 ; pop_j < this->population_number(); pop_j++ ){
            //double pop_ratio = this->population_size(pop_j) / scaling_pop_size_N0;          
            //recomb_opportunity += this->total_weighted_recomb_opportunity[time_layer_i][pop_j] / pop_ratio ;
            //this->recomb_count_ += this->total_recomb_count[time_layer_i][pop_j] / pop_ratio ;
            //}
        //// Advance to the next interval level
        //if ( current_time_idx_ == change_times_.size() - 1) break;  
        //this->increaseTime(); 
        //}
    //this->inferred_recomb_rate = this->recomb_count_ / recomb_opportunity;
    //// reset the inferred recombination rate in the current Model and CountModel?
    //}



//void CountModel::extract_and_update_count(ParticleContainer &Endparticles, double current_base, bool end_data ){
    //for (size_t i = 0; i < Endparticles.particles.size(); i++){                
        //ForestState* counting_state = Endparticles.particles[i];
        //double weight = counting_state->weight();
        //for ( size_t time_i = this->change_times_.size() - 1 ; (int)time_i >= 0 ; time_i --){
            //this->update_coal_count ( counting_state->CoaleventContainer[time_i] , time_i, weight);             
            //this->update_recomb_count ( counting_state->RecombeventContainer[time_i] , time_i, weight); 
            //this->update_migr_count ( counting_state->MigreventContainer[time_i] , time_i, weight);           
            //previous_base[time_i] = current_base - lags[time_i] > 0 ? current_base - lags[time_i] : (double)0;
            //}  
        //}
    //}



//void CountModel::resize_Starevent ( deque < Starevent *> & StareventContainer_i , int index) {
    ////cout<<"size was "<<StareventContainer_i.size();
    //size_t i = 0;
    //while (i < index && StareventContainer_i.size()>0 ){
        //StareventContainer_i[0]->pointer_counter_ --;
        //if ( StareventContainer_i[0]->pointer_counter_ == 0){
            //delete  StareventContainer_i[0];
            //}
        //StareventContainer_i.pop_front();
        //i++;
        //}
        ////cout<<" now is "<<StareventContainer_i.size()<<endl;
    //}

//void CountModel::resize_Migrevent ( deque < Migrevent *> & StareventContainer_i , int index) {
    ////cout<<"size was "<<StareventContainer_i.size();
    //size_t i = 0;
    //while (i < index && StareventContainer_i.size()>0 ){
        //StareventContainer_i[0]->pointer_counter_ --;
        //if ( StareventContainer_i[0]->pointer_counter_ == 0){
            //delete  StareventContainer_i[0];
            //}
        //StareventContainer_i.pop_front();
        //i++;
        //}
        ////cout<<" now is "<<StareventContainer_i.size()<<endl;
    //}



//void CountModel::update_migr_count ( deque < Migrevent * > & MigreventContainer_i, size_t time_i, double weight ){
    //double x_start = (double)this->previous_base[time_i];    
    //int Migrevent_index = MigreventContainer_i.size()-1;
    //while ( Migrevent_index > 0 ){
        //Migrevent* current_Migrevent = MigreventContainer_i[Migrevent_index];
        //if ( current_Migrevent->base() < x_start){
            //break;
            //}
        //if (current_Migrevent->event_state() == EVENT){
            //this->total_mig_count[time_i][current_Migrevent->pop_i()][current_Migrevent->mig_pop()] += current_Migrevent->num_event() * weight;
            //} 
        //for (size_t potential_pop = 0; potential_pop < this->total_weighted_mig_opportunity[time_i][current_Migrevent->pop_i()].size(); potential_pop++){
            //this->total_weighted_mig_opportunity[time_i][current_Migrevent->pop_i()][potential_pop] += current_Migrevent->opportunity() * weight;    
            //}    
        //Migrevent_index--;                
        //}  
    //resize_Migrevent(MigreventContainer_i, Migrevent_index);
    //}


//void CountModel::update_coal_count ( deque < Starevent *> & CoaleventContainer_i, size_t time_i, double weight ){
    //double x_start = (double)this->previous_base[time_i];
    //int Coalevent_index = CoaleventContainer_i.size()-1;
    //while ( Coalevent_index > 0 ) {        
        //Starevent * current_Coalevent = CoaleventContainer_i[Coalevent_index];
        //if ( current_Coalevent->base() < x_start){
            //break;
            //}
        //this->total_coal_count[time_i][current_Coalevent->pop_i()]                += weight * current_Coalevent->num_event();
        //this->total_weighted_coal_opportunity[time_i][current_Coalevent->pop_i()] += weight * current_Coalevent->opportunity();            
        //Coalevent_index--;
        //} 
    
    //resize_Starevent(CoaleventContainer_i, Coalevent_index);        
        
    //}    

//void CountModel::update_recomb_count ( deque < Starevent* > & RecombeventContainer_i, size_t time_i, double weight ){
    //double x_start = (double)this->previous_base[time_i];    
    //int Recombevent_index = RecombeventContainer_i.size()-1;
    //while ( Recombevent_index > 0 ){
        //Starevent* current_Recombevent = RecombeventContainer_i[Recombevent_index];
        //if ( current_Recombevent->base() < x_start){
            //break;
            //}
        //this->total_recomb_count[time_i][current_Recombevent->pop_i()]                += weight * current_Recombevent->num_event();
        //this->total_weighted_recomb_opportunity[time_i][current_Recombevent->pop_i()] += weight * current_Recombevent->opportunity();            
        //Recombevent_index--;                
        //}
    //resize_Starevent(RecombeventContainer_i, Recombevent_index); 
    //}




//useless functions
double CountModel::extract_and_update_count(ParticleContainer &Endparticles, double current_base, bool end_data ){
    //for ( size_t time_i = 0 ; time_i < this->change_times_.size(); time_i ++){
    for ( size_t time_i = this->change_times_.size() - 1 ; (int)time_i >= 0 ; time_i --){
        //cout << "at time level " << time_i << "current_base " << current_base << " this->lags[time_i] " << this->lags[time_i] <<endl;
        double x_end =  (double)this->previous_base[time_i] < ( current_base - this->lags[time_i] ) ? ( current_base - this->lags[time_i] ) : (double)this->previous_base[time_i] ;
        if (end_data){
            x_end = current_base;
            }
        for (size_t pop_j = 0 ; pop_j < this->population_number(); pop_j++ ){
            this->count_events_in_one_interval(Endparticles, time_i, pop_j, previous_base[time_i], x_end );                
            }
        previous_base[time_i] = current_base - lags[time_i] > 0 ? current_base - lags[time_i] : (double)0;
        }  
    
    double remove_particle_before_site = previous_base[0];
    for (size_t i = 0 ; i < previous_base.size() ; i++ ){
        remove_particle_before_site = remove_particle_before_site < previous_base[i] ? remove_particle_before_site : previous_base[i];
        }
    dout<< "remove_particle_before_site = "<<remove_particle_before_site<<endl;
    return remove_particle_before_site;
    }


void CountModel::count_events_in_one_interval(ParticleContainer &Endparticles, size_t time_i, size_t pop_j, double x_start, double x_end){
    if ( x_start == x_end ){ return; }    
    for (size_t i = 0; i < Endparticles.particles.size(); i++){                
        ForestState* counting_state = Endparticles.particles[i];
        double weight = counting_state->weight();
        // Skip a few state between lagging until the most updated case
        while (counting_state -> current_base() >= x_end && counting_state->previous_state){
            counting_state = counting_state->previous_state;                
        }
        
        // Start counting
        while ( counting_state -> current_base() >= x_start ) {  // Making sure there is coalescent events between the interval
            for ( size_t event_i = 0 ; event_i < counting_state->CoaleventContainer.size() ; event_i++ ){
                Coalevent * current_Coalevent = counting_state->CoaleventContainer[event_i];
                /*! Cumulate the coalescent events if the event is within the interval 
                 */
                if ( time_i == current_Coalevent->change_time_i() && pop_j == current_Coalevent->pop_i() ){
                    this->total_coal_count[time_i][current_Coalevent->pop_i()]                    += weight * current_Coalevent->num_event();
                    this->total_weighted_coal_opportunity[time_i][current_Coalevent->pop_i()]     += weight * current_Coalevent->opportunity();            
                    //current_Coalevent->set_counted(true);
                    }
                } //  < counting_state->CoaleventContainer.size() 

            for ( size_t event_i = 0 ; event_i < counting_state->RecombeventContainer.size() ; event_i++ ){
                Recombevent * current_Recombevent = counting_state->RecombeventContainer[event_i];
                /*! Cumulate the recombination events if the event is within the interval 
                 */
                if ( time_i == current_Recombevent->change_time_i() && pop_j == current_Recombevent->pop_i() ){
                    this->total_recomb_count[time_i][current_Recombevent->pop_i()]                    += weight * current_Recombevent->num_event();
                    this->total_weighted_recomb_opportunity[time_i][current_Recombevent->pop_i()]     += weight * current_Recombevent->opportunity();            
                    //current_Recombevent->set_counted(true);
                    }
                } //  < counting_state->RecombeventContainer.size() 
                
            for ( size_t j = 0; j < counting_state->MigreventContainer.size(); j++){
                Migrevent * current_Migrevent = counting_state->MigreventContainer[j];

                /*! Cumulate the recombination events if the event is within the interval 
                 */
                if (current_Migrevent->event_state() == EVENT){
                    this->total_mig_count[current_Migrevent->pop_i()][current_Migrevent->mig_pop()] += current_Migrevent->num_event() * weight;
                    } 
                for (size_t potential_pop = 0; potential_pop < this->total_weighted_mig_opportunity[current_Migrevent->pop_i()].size(); potential_pop++){
                    this->total_weighted_mig_opportunity[current_Migrevent->pop_i()][potential_pop] += current_Migrevent->opportunity() * weight;    
                    }    
                } //  < counting_state->MigreventContainer.size() 
                                
                counting_state = counting_state->previous_state;      
                if (!counting_state) break;                      
                
            }  // End of while loop: counting_state -> current_base() >= x_start
        } //  End of for loop: < Endparticles.particles.size()
        return;
    } // 


void CountModel::count_events_in_one_interval_alt(ParticleContainer &Endparticles, size_t time_i, size_t pop_j, double x_start, double x_end){
    if ( x_start == x_end ){ return; }
    for (size_t i = 0; i < Endparticles.particles.size(); i++){
        ForestState* counting_state = Endparticles.particles[i];
        double weight = counting_state->weight();
        // Skip a few state between lagging until the most updated case
        while (counting_state -> current_base() >= x_end && counting_state->previous_state){
            counting_state = counting_state->previous_state;                
        }
        
        // Start counting
        while ( counting_state -> current_base() >= x_start ) {  // Making sure there is coalescent events between the interval
            while (!counting_state->CoaleventContainer.empty()){
                Coalevent * current_Coalevent = counting_state->CoaleventContainer.back();
                /*! Cumulate the coalescent events if the event is within the interval 
                 */
                if ( time_i == current_Coalevent->change_time_i() && pop_j == current_Coalevent->pop_i() ){
                    this->total_coal_count[time_i][current_Coalevent->pop_i()]                    += weight * current_Coalevent->num_event();
                    this->total_weighted_coal_opportunity[time_i][current_Coalevent->pop_i()]     += weight * current_Coalevent->opportunity();            
                    delete current_Coalevent;
                    counting_state->CoaleventContainer.pop_back();
                    }
                else {
                    break;
                    } 
                } //  < counting_state->CoaleventContainer.size() 

            while (!counting_state->RecombeventContainer.empty()){
                Recombevent * current_Recombevent = counting_state->RecombeventContainer.back();
                /*! Cumulate the recombination events if the event is within the interval 
                 */
                if ( time_i == current_Recombevent->change_time_i() && pop_j == current_Recombevent->pop_i() ){
                    this->total_recomb_count[time_i][current_Recombevent->pop_i()]                    += weight * current_Recombevent->num_event();
                    this->total_weighted_recomb_opportunity[time_i][current_Recombevent->pop_i()]     += weight * current_Recombevent->opportunity();            
                    delete current_Recombevent;
                    counting_state->RecombeventContainer.pop_back();
                    }
                else {
                    break;
                    } 
                } //  < counting_state->RecombeventContainer.size() 
                
            for ( size_t j = 0; j < counting_state->MigreventContainer.size(); j++){
                Migrevent * current_Migrevent = counting_state->MigreventContainer[j];

                /*! Cumulate the recombination events if the event is within the interval 
                 */
                if (current_Migrevent->event_state() == EVENT){
                    this->total_mig_count[current_Migrevent->pop_i()][current_Migrevent->mig_pop()] += current_Migrevent->num_event() * weight;
                    } 
                for (size_t potential_pop = 0; potential_pop < this->total_weighted_mig_opportunity[current_Migrevent->pop_i()].size(); potential_pop++){
                    this->total_weighted_mig_opportunity[current_Migrevent->pop_i()][potential_pop] += current_Migrevent->opportunity() * weight;    
                    }    
                } //  < counting_state->MigreventContainer.size() 
                                
                counting_state = counting_state->previous_state;      
                if (!counting_state) break;                      
                
            }  // End of while loop: counting_state -> current_base() >= x_start
        } //  End of for loop: < Endparticles.particles.size()
        return;
    } // 



//void PfParam::log_end(pfTime * running_time){
    //ofstream log_file ( log_NAME.c_str(), ios::out | ios::app | ios::binary); 
    //log_file << "Initial particles building took about " << running_time->timing_[0] << " second(s) \n";
    ////log_file << "Simulation took about " << end_time - initial_state_end_time << " second(s) \n";
    //log_file << "    Update weight took " << running_time->timing_[2]<<" second(s)\n";
    //log_file << "    Resampling took " << running_time->timing_[1]<<" second(s)\n";
    //log_file.close();
    //}




//size_t CountModel::find_time_interval (double start_height, double end_height){
    //double pop_start_height = 0.0; 
    ////size_t time_i = 0;
    //auto time_i = 0;
    //for (; time_i < change_times_.size() ; time_i++){
        //double pop_end_height = (time_i == this->change_times_.size()-1) ? FLT_MAX : change_times_[time_i+1];
        //if ( log(pop_start_height) <= log(start_height) + std::numeric_limits<double>::epsilon()  && 
             //log(end_height) <= log(pop_end_height) + std::numeric_limits<double>::epsilon() ) {
            //break;
            //}
        //pop_start_height = pop_end_height;
        //} // End of for loop: < change_times.size()  
        //return time_i;
    //} // End of function: CountModel::find_time_interval( ... )



//char *convert(const std::string & s){
    //char *pc = new char[s.size()+1];
    //std::strcpy(pc, s.c_str());
    //return pc; 
    //}


//double update_cum_sum_array_find_ESS(Vparticle * pointer_to_current_particles,std::valarray<double> & weight_cum_sum,size_t Num_of_particles){
    
    //double wi_sum=0;
    //double wi_sq_sum=0;
    //weight_cum_sum=0; //Reinitialize the cum sum array 
    //for (size_t i=0; i<Num_of_particles ;i++){
    ////for (size_t i=0; i<pointer_to_current_particles->particles.size() ;i++){
        ////update the cum sum array
        //double w_i=pointer_to_current_particles->particles[i]->weight();
        //weight_cum_sum[i+1]=weight_cum_sum[i]+w_i;
        //wi_sum = wi_sum + w_i;
        ////wi_sq_sum = wi_sq_sum + pow(pointer_to_current_particles->particles[i]->weight,2);
        //wi_sq_sum = wi_sq_sum + w_i * w_i;
        ////cout<<weight_cum_sum[i]<<endl;
    //}
    ////check for the cum weight
    //dout << "### particle weights ";for (size_t i=0;i<Num_of_particles;i++){dout<<pointer_to_current_particles->particles[i]->weight()<<"  ";}dout<<std::endl<<std::endl;
    //dout << "### updated cum sum of particle weight ";for (size_t i=0;i<weight_cum_sum.size();i++){dout<<weight_cum_sum[i]<<"  ";}dout<<std::endl;
    //return wi_sum * wi_sum / wi_sq_sum;
//}

void online_reset_model_Ne(Model * model, vector <Pfcount> countContainer, vector<double>* change_times, bool online, bool print){
    if (online){
        reset_model_Ne(model, countContainer, change_times, print);
    }
    return;    
}

/*! \ingroup group_Ne */
void reset_model_Ne(Model * model, vector <Pfcount> countContainer, vector<double>* change_times, bool print){ //change countContainer to pointer
    model->resetTime(); 
    
    double total_coal_count ;
    double total_weighted_BL ;
    if (print){
    cout<<std::setw(15)<<"time" << std::setw(15)<<"count" << std::setw(15)<<"popsize" << endl;
    }
    for (size_t j=0; j < change_times->size();){
        total_coal_count=0;
        total_weighted_BL=0;
        for (size_t i=0;i<countContainer.size();i++){
            total_coal_count += countContainer[i].weightedCoalescent_[j];
            total_weighted_BL += countContainer[i].weightedBL_[j];
        }    
        double Nehat = total_weighted_BL / total_coal_count /2 ;
        
        if (print){
            //cout << "xstart = " << countContainer
            
            cout << std::setw(15) << total_weighted_BL
                 << std::setw(15) << total_coal_count
                 << std::setw(15) << Nehat<<endl;
        }
        
        dout << "Time interval " << j+1<<" " 
        <<"total coalescent = " << total_coal_count
        <<", total BL = " << total_weighted_BL
        <<", True Ne = " << model->population_size()
        <<", Ne hat= " <<  Nehat <<endl;
        
        model->addPopulationSizes(change_times->at(j), Nehat ,false,false);
        j++;
        
        if ( j<change_times->size() ){
            model->increaseTime();
        }
    }
    model->resetTime(); 
    return ;
}


//void systemetic_resampling(std::valarray<double> cum_sum, std::valarray<double> even_space, std::valarray<int> & sample_count){
    //cum_sum=cum_sum/cum_sum[cum_sum.size()-1];
    //double u=unifRand()/sample_count.size();
    //size_t interval_j=0;
    //size_t sample_i=0;
    //double uj=0;

    //dout<<std::endl<<"systematic sampling procedue" << std::endl;
    //for (size_t i=0;i<cum_sum.size();i++){dout<<cum_sum[i]<<"  ";}dout<<std::endl;
    
    //while ((interval_j<sample_count.size()) && (sample_i<sample_count.size())){//remove the second condition...
        //uj=u+even_space[sample_i];
        //dout << "Is " <<  uj<<" in the interval of " << std::setw(10)<<cum_sum[interval_j]<< " and " << std::setw(10)<< cum_sum[interval_j+1];
        //if ( (cum_sum[interval_j]<uj) && (uj<cum_sum[interval_j+1])){
            //sample_count[interval_j]++;
            //sample_i++;
            //dout << "  yes, update sample count of particle " << interval_j<<" to " << sample_count[interval_j] <<std::endl;
        //}
        //else{
            //interval_j++;
            //dout << "   no, try next interval " << std::endl;
        //}
    //}
    //dout << "systematic sampling procedue finished with total sample count " << sample_count.sum()<<std::endl<<std::endl;

//}



void check_copy_forest(Vparticle * pointer_to_current_particles,Model* model, RandomGenerator* rg,size_t Num_of_particles){
        dout << " --------------------    Initial Particles   --------------------" << std::endl;
        
        Aparticle *  new_particle = new Aparticle(model,rg);
        pointer_to_current_particles->push(new_particle);
    pointer_to_current_particles->particles[0]->particle_forest->sampleNextGenealogy();
    pointer_to_current_particles->particles[0]->particle_forest->sampleNextGenealogy();
    pointer_to_current_particles->particles[0]->particle_forest->printTree();
    pointer_to_current_particles->particles[0]->particle_forest->printTree();

    for (size_t i=1; i<Num_of_particles ;i++){
        dout << "****************** Particle # " << i << "  ******************" << std::endl;
        
        Aparticle * new_copy_particle= new Aparticle(pointer_to_current_particles->particles[0]);
        pointer_to_current_particles->push(new_copy_particle);
        pointer_to_current_particles->particles[i]->particle_forest->printTree();
        //even_space[i]=double(i)/double(Num_of_particles); // used in resampling later
    }
    for (size_t i=0; i<Num_of_particles ;i++){
        dout << "****************** Particle # " << i << "  ******************" << std::endl;
        
        pointer_to_current_particles->particles[i]->particle_forest->sampleNextGenealogy();
        //pointer_to_current_particles->particles[i]->particle_forest.printTree();
        //even_space[i]=double(i)/double(Num_of_particles); // used in resampling later
    }    

    //Vparticle* pointer_to_current_particles=&init_p_array;
    pointer_to_current_particles->makeHeap_();
    dout<<std::endl<<"Initial particles of " << Num_of_particles<<"  were built " << std::endl;        
    dout << " --------------------    Initial Particles End  --------------------" << std::endl<<std::endl;
            
    dout << "  Check the order of the particles, the particle with the minimal sequence length must be the first!!!  " << endl;
    //check particles order, the particle with the minimal sequence length must be the first!!!
    pointer_to_current_particles->check_particle_orders();
    
}


//pfARG::Vparticle * updating_current_particles(pfARG::Vparticle * pointer_to_current_particles, valarray<int> & sample_count, valarray<double> & weight_cum_sum, double curly_L, size_t Num_of_particles){
void updating_current_particles(pfARG::Vparticle * pointer_to_current_particles, valarray<int> & sample_count, valarray<double> & weight_cum_sum, double curly_L, size_t Num_of_particles){
////void updating_current_particles(pfARG::Vparticle * pointer_to_current_particles, valarray<int> & sample_count, valarray<double> & weight_cum_sum, valarray<double> & even_space, double curly_L, size_t Num_of_particles){
    //set the current site for the simulation, and sample the new genealogy

    Vparticle * vector_of_new_particles=new Vparticle(pointer_to_current_particles, sample_count);       
    //call delete and remove the particles that are not resampled....

    pointer_to_current_particles->particles.clear();
    //delete pointer_to_current_particles;// THIS IS NOT WORKING!!! introduces more errors???!!! Check for the destructor of forest!!!
    pointer_to_current_particles=NULL;          
    pointer_to_current_particles=vector_of_new_particles;
    //vector_of_new_particles=NULL;    
    //making heap and the minimal 
    pointer_to_current_particles->makeHeap_();
    //checking for heap update
    pointer_to_current_particles->check_particle_orders();
    ////cout<<vector_of_new_particles->particles.size()<<"  " << pointer_to_current_particles->particles.size()<<endl;
    
    //delete pointer_to_current_particles;
    
    //vector_of_new_particles->makeHeap_();
    //vector_of_new_particles->check_particle_orders();
    
    //return vector_of_new_particles;
}


void Aparticle::record_coalevent(const TimeInterval & current_event, double end_time){
    this->coaleventContainer.back()->set_num_coal(1);
    this->coaleventContainer.back()->set_end_height(end_time);
    if (current_event.numberOfContemporaries()>1){
        this->coaleventContainer.back()->set_num_lineages(current_event.numberOfContemporaries());
    }
    else{
        this->coaleventContainer.back()->set_num_lineages(1);
    }
    dout << "* * * * Events recorder is called start at " << this->coaleventContainer.back()->start_height()
    <<" with " << this->coaleventContainer.back()->num_lineages() <<" lineages"
    <<", with " << this->coaleventContainer.back()->num_coal()<<" coalescent at "
    <<this->coaleventContainer.back()->end_height()<<endl;
}

void Aparticle::record_recombevent(const TimeInterval & current_event, double end_time){
    this->coaleventContainer.back()->set_num_recomb(1);
    this->coaleventContainer.back()->set_end_height(end_time);
    if (current_event.numberOfContemporaries()>1){
        this->coaleventContainer.back()->set_num_lineages(current_event.numberOfContemporaries());
    }
    else{
        this->coaleventContainer.back()->set_num_lineages(1);
    }
    dout << "* * * * Events recorder is called start at " << this->coaleventContainer.back()->start_height()
    <<" with " << this->coaleventContainer.back()->num_lineages() <<" lineages"
    <<", with " << this->coaleventContainer.back()->num_recomb()<<" recombination at "
    <<this->coaleventContainer.back()->end_height()<<endl;

}

void Aparticle::record_migevent(const TimeInterval & current_event, double end_time){
    this->coaleventContainer.back()->set_num_mig(1);
    this->coaleventContainer.back()->set_end_height(end_time);
    if (current_event.numberOfContemporaries()>1){
        this->coaleventContainer.back()->set_num_lineages(current_event.numberOfContemporaries());
    }
    else{
        this->coaleventContainer.back()->set_num_lineages(1);
    }
    dout << "* * * * Events recorder is called start at " << this->coaleventContainer.back()->start_height()
    <<" with " << this->coaleventContainer.back()->num_lineages() <<" lineages"
    <<", with " << this->coaleventContainer.back()->num_mig()<<" migration at "
    <<this->coaleventContainer.back()->end_height()<<endl;

}


double Aparticle::cal_sumofiTi() const{
    double sum_of_iTi = 0;
    //cout << "sum_of_iTi= ";
    for (size_t i = 0; i < this->coaleventContainer.size() ; i++){
        sum_of_iTi = sum_of_iTi +  coaleventContainer[i]->num_lineages() * coaleventContainer[i]->length() ;
        //cout<<sum_of_iTi<<" ";
    }
    //cout<<endl;
    //if (sum_of_iTi==0){
        //cout<<endl<<"event size =" << this->events.size()<<". ";
        //for (size_t i = 0; i < this->events.size() ; i++){
        //cout<<i<<" " << events[i]->numberOfContemporaries()<<"*(" <<  events[i]->end_height()<<" - " <<  events[i]->start_height()<<")=" << events[i]->numberOfContemporaries() * ( events[i]->end_height() - events[i]->start_height() )<<", ";
        //}
        ////cout<<endl;
        //cout<<sum_of_iTi<<endl;
    //}
    return sum_of_iTi;
} 


/*! Clear coalescent events recorded when the initial particles are constructed. */
void Aparticle::clear_initial_Coalevent(){
    this->CoaleventContainer.clear();
    Coalevent* new_event = new Coalevent();
    this->CoaleventContainer.push_back(new_event);
    };
    
void pfARG::param::wright_R_script(){
    //draw TMRCA plot 
    remove(R_NAME.c_str());
    ofstream R_file;
    R_file.open (this->R_NAME.c_str(), ios::out | ios::app | ios::binary); 
    R_file<<"rm(list=ls()) \n library(stats) \n";
    R_file<<"vcf_file= \"" << this->vcf_NAME<<"\"\n vcf_data=read.table(vcf_file) \n hetro_base=vcf_data$V2[vcf_data$V7==\"PASS\"] \n";
    R_file<<"TMRCA_file= \"" << this->TMRCA_NAME<<"\"\n mydata=read.table(TMRCA_file) \n n_row=dim(mydata)[1] \n n_col=dim(mydata)[2]-1 \n ";
    R_file<<"base=mydata[,1] \n base_min=min(base) \n base_max=max(base) \n";
    //R_file<<"TMRCA=as.matrix(mydata[,-1]) \n TMRCA_mean=rowMeans(mydata[,-1]) \n TMRCA_min=min(TMRCA_mean) \n TMRCA_max=max(TMRCA_mean) \n";
    R_file<<"TMRCA=as.matrix(mydata[,-1]) \n TMRCA_mean=rowMeans(mydata[,-1]) \n TMRCA_min=min(TMRCA) \n TMRCA_max=max(TMRCA) \n";
    R_file<<"pdf(\"" << this->TMRCA_NAME<<".pdf\") \n ";
    R_file<<"plot(c(base_min,base_max),c(TMRCA_min*0.9,TMRCA_max),xlab=\"base\",ylab=\"TMRCA\",main=\"" << this->N<<" particles\",type=\"n\")\n";
//R_file<<"plot(c(base_min,base_max),c(TMRCA_min*0.9,TMRCA_max),xlab=\"base\",ylab=\"TMRCA\",main=\" << " << this->nsam<<" sample " << this->N<<" particles " << \",type=\'n\')\n";
    R_file<<"lines(base, TMRCA_mean, type=\"S\") \n ";
//#par(cin=10)
    //R_file<<"for (i in 1:(n_row-1)){ \n \t ME = 1.96 * sqrt(var(TMRCA[i,])/n_col) \n \t lines(c(base[i],base[i]),c(TMRCA_mean[i]+ME,TMRCA_mean[i]-ME),type=\"p\",pch=20,col=\"blue\")  \n \t lines(c(base[i],base[i]),c(TMRCA_mean[i]+ME,TMRCA_mean[i]-ME),type=\"l\",pch=20,col=\"blue\")  \n} \n";
    R_file<<"for (i in 1:(n_row-1)){ \n points(rep(base[i],n_col),TMRCA[i,],pch=\".\",col=\"blue\") \n} \n";
//#lines(c(base[i],base[i+1],base[i+1]),c(TMRCA_mean[i],TMRCA_mean[i],TMRCA_mean[i+1]), type="l")
    R_file<<"lines(hetro_base,rep(TMRCA_min*0.9,length(hetro_base)),type=\"p\",col=\"red\") \n dev.off() \n";
    
    //if (this->N < 150){
        R_file<<"pdf(\"" << this->TMRCA_NAME<<"-heat.pdf\") \n ";
        R_file<<"x=base \n y=TMRCA[,1] \n ";
        R_file<<"for (i in 2:(n_col)){ \n \t x=c(x,base) \n \t y=c(y,TMRCA[,i]) \n } \n";
        R_file<<"df=data.frame(x,y) \n k <- with(df,MASS:::kde2d(x,y)) \n filled.contour(k,color=heat.colors) \n dev.off() \n ";
    //}
    R_file.close();    
}






///*! 
 //* SET the model parameters, obtain number of samples from vcf file
 //*     if vcf file is empty, use default number of samples
 //*     if mutation rate is not given, use default mutation rate
 //*     if recombination rate is not given, use default recombination rate
 //*/ 
//void initialize_model(Model* model, 
                    //Param * scrm_para,
                    //pfARG::param pfARG_para,
                    //Vcf * VCFfile){
    
    //scrm_para->parse(*model);
    //size_t nsam = 2*VCFfile->nsam(); //overwrite the sample size according the data.
    
    //if (!VCFfile->withdata()){
        //nsam = pfARG_para.default_nsam;
        ////cout<<"scrm loci length is "<<model->loci_length()<<endl;
        //if ( model->loci_length() == 100000 && model->recombination_rate() == 0){
            //VCFfile->set_event_interval( pfARG_para.default_loci_length / 10 );
        //} else {
            //VCFfile->set_event_interval( model->loci_length() / 10 );
        //}
    //}    
    
    //if (model->sample_size() == 0) { // No subpopulation
        //model->addSampleSizes(0.0, vector<size_t>(1, nsam));
    //} 
    
    //else if ( model->sample_size() != 2*VCFfile->nsam() && VCFfile->withdata() ){ // Subpopulation, total number of samples do not match number of sample given by vcf file
            //throw std::invalid_argument("Total number of samples in subpopulation do not match number of samples in vcf file");    
    //}
    
    ////else if (nsam !=  model->sample_size()) { 
        ////cout  << "model->sample_size() = "<<model->sample_size()<<endl;
        ////cout <<"nsam = "<<nsam<<endl;        
        ////throw std::invalid_argument("No data given, and sum of samples in subpopulation not equal to the total sample size");    
    ////}
    
    ////check mutation rate and recombination rate    
    //if (model->recombination_rate() == 0){
        //model->set_recombination_rate(pfARG_para.default_recomb_rate, pfARG_para.default_loci_length);
    //}    

    //if (model->mutation_rate() == 0){
        //model->set_mutation_rate(pfARG_para.default_mut_rate);
    //}
    //model->finalize();
    //cout<<"model is finalized"<<endl;
    //return;
//}


