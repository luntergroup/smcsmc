//////////////////////////////////////////////////////////////////////////////////////////////////////
//                                NOT USED, FOR REMOVAL
//////////////////////////////////////////////////////////////////////////////////////////////////////
        void print_Usage_initial_particles(
            size_t ten_percent_ofN,
            size_t total_num_of_nodes,
            struct rusage *p
            );
//////////////////////////////////////////////////////////////////////////////////////////////////////
//                                NOT USED, FOR REMOVAL
//////////////////////////////////////////////////////////////////////////////////////////////////////       
        void moveDown_(size_t node);
        void moveUp_(size_t node);


//////////////////////////////////////////////////////////////////////////////////////////////////////
//                                NOT USED, FOR REMOVAL
//////////////////////////////////////////////////////////////////////////////////////////////////////


/*! \brief Particle filtering Initialization Part II
 * Create a larger initial particle pool, and resample from the pool to obtain the intial particles
 * 
 * \ingroup group_pf_init
 */ 
Vparticle::Vparticle(
    Vparticle * pointer_to_current_particles, 
    double L_hat, // 
    size_t post_initial_num_particles, 
    Model* model, RandomGenerator* rg, 
    vector <bool> data_for_init_particles,
    bool withdata
){
	this->max_weight=0;
	size_t ten_percent_ofN = post_initial_num_particles/10;
	
	struct rusage usage;
	struct rusage *p=&usage;
	
	size_t total_num_of_nodes=0;
	
// Two stages, first to recycle whichever particles that satisfy the rejection resampling step, then for whatever many left, create new particles 	
	
// stage 1
	for (size_t i=0; i<pointer_to_current_particles->particles.size();i++){
		
		if (unifRand() < (pointer_to_current_particles->particles[i]->weight() / L_hat)){// Rejection sampling step
			Aparticle * new_copy_particle= new Aparticle(*pointer_to_current_particles->particles[i]);
			new_copy_particle->previous_particle=NULL;
			this->push(new_copy_particle);
			
			total_num_of_nodes += new_copy_particle->getNodes()->size();
		
			this->print_Usage_initial_particles(ten_percent_ofN, total_num_of_nodes, p);
		}
	}
	
// stage 2
	//double mutation_rate = model->mutation_rate();
	while (this->particles.size() < post_initial_num_particles){
		
		Aparticle *  new_particle = new Aparticle(model,rg);
		new_particle->include_haplotypes_at_tips(data_for_init_particles);

		if (unifRand() < new_particle->calculate_likelihood(withdata) / L_hat){
			
			this->push(new_particle);
			
			total_num_of_nodes += new_particle->getNodes()->size();
			
			//this->print_Usage_initial_particles(ten_percent_ofN, total_num_of_nodes, p);
		}
		else{
			delete new_particle;
		}
	}
}


//void create_a_bigger_pool_of_particles(){
            //  First create a bigger list of initial particles pool with size 100 times the number of particles needed, then rescale  

        //int initial_num_particles = pfARG_para.N;      

        //Vparticle * init_pointer_to_current_particles=new Vparticle(model, rg, initial_num_particles, VCFfile->vec_of_sample_alt_bool, VCFfile->empty_file());            
        
        //double L_hat=cal_L_hat(init_pointer_to_current_particles, pfARG_para.pool_size);

        //int post_initial_num_particles = 1.2 * pfARG_para.N;
        ////assert(post_initial_num_particles>initial_num_particles);    

        //Vparticle * post_init_pointer_to_current_particles=new Vparticle(init_pointer_to_current_particles, L_hat, post_initial_num_particles, model, rg, VCFfile->vec_of_sample_alt_bool, VCFfile->empty_file());
        
        //valarray<double> initial_weight_cum_sum((post_initial_num_particles+1)); //Initialize the weight cumulated sums        
        //update_cum_sum_array(post_init_pointer_to_current_particles,initial_weight_cum_sum,post_initial_num_particles);            
        //valarray<int> initial_sample_count(post_initial_num_particles);    // if sample_count is in the while loop, this is the initializing step...
        //systemetic_resampling( initial_weight_cum_sum, initial_sample_count, pfARG_para.N);
        
        //for (size_t i=0;i<initial_sample_count.size();i++){dout<<initial_sample_count[i]<<"  ";}dout<<std::endl;
        
        ////Actual initial particles, next three lines

        //Vparticle * pointer_to_current_particles=new Vparticle(post_init_pointer_to_current_particles, initial_sample_count);             
        //for (size_t particle_i=0;particle_i<current_particles.particles.size();particle_i++){
            //current_particles.particles[particle_i]->previous_particle=NULL;
            //current_particles.particles[particle_i]->pointer_counter=0;
        //}


        ////delete init_pointer_to_current_particles;
        //delete post_init_pointer_to_current_particles; //DEBUG this should be removed as well. check...
        
        ////double printing_base = VCFfile->site();       
        ////pfARG_para.appendingStuffToFile(printing_base, pointer_to_current_particles,pfARG_para.lag);
        
        ////runningtime->stopwatch_end(0);
    //}






//void Vparticle::print_2_individual_initial_pop_est(Model * model){
	//double totalbl=0;
	//for (size_t i=0;i<this->particles.size();i++){
		//totalbl+=this->particles[i]->local_root()->height();;
	//}
	//cout << "True Ne = " << model->population_size()
    //<<", Ne hat= " <<  totalbl/this->particles.size() /2 <<endl;
//}



void Vparticle::moveDown_(size_t node)
{
  size_t left_child = node*2+1;
  size_t right_child = node*2+2;
 
  size_t replace = node;
  if (right_child < particles.size())
  {
    bool left = particles[right_child]->next_base() > particles[left_child]->next_base();
    if (left && particles[node]->next_base() > particles[left_child]->next_base())
      replace = left_child;
    else if (!left && particles[node]->next_base() > particles[right_child]->next_base())
      replace = right_child;
  }
  else if (left_child < particles.size())
  {
    if (particles[node]->next_base() > particles[left_child]->next_base())
      replace = left_child;
  }
 
  if (replace == node)
    return;
  std::swap(particles[node], particles[replace]);
  moveDown_(replace);
}
 

void Vparticle::moveUp_(size_t node)
{
  if (node == 0)
    return;
 
  size_t parent = std::floor((node-1)/2);
 
  if (particles[node]->next_base() > particles[parent]->next_base())
    return;
  std::swap(particles[node], particles[parent]);
  moveUp_(parent);
}

void Vparticle::makeHeap_()
{
    //dout << "particles.size() is " << particles.size()<<std::endl;
    //dout << "particles.size()/2 is " << particles.size()/2<<std::endl;
	//throw std::length_error("size ");
	for (int i = particles.size()/2; i >= 0; i--){
		//dout << "moving particle " << i<<std::endl;
		moveDown_(unsigned(i));
	}
}
 
/*!
 * \ingroup group_resource
 */
void Vparticle::print_Usage_initial_particles(
    size_t ten_percent_ofN,
    size_t total_num_of_nodes,
    struct rusage *p
){
	//int ret;
	int who= RUSAGE_SELF;
	if (this->particles.size() % ten_percent_ofN == 0){ 
		cout <<"Initial building " <<  this->particles.size() << " particles " ;
		//int ret=
		getrusage(who,p);
		//process(p, "usage: "); 
		cout << "average number of nodes: " << total_num_of_nodes / this->particles.size() << endl;
	}
}


//Vparticle::Vparticle(In first, In last)
//{
//  size_t n = std::distance(first, last);
//  particles.reserve(n);
//  for (In in = first; in != last; ++in)
//    particles.push_back(*in);
//  makeHeap_();
//}

void Vparticle::insert(Aparticle* new_particle)
{
  particles.push_back(new_particle);
  moveUp_(particles.size()-1);
}


void pfARG::param::appending_log_file(string log_file_input){
    ofstream log_file;
    log_file.open (log_NAME.c_str(), ios::out | ios::app | ios::binary); 
    log_file << log_file_input << "\n";
    log_file.close();
}        
        //Vparticle(Vparticle * pointer_to_current_particles, valarray<int> & sample_count);// this is used to create a new vector of new particles, that are based on the previous particles. 
        //Vparticle(Vparticle * pointer_to_current_particles, double L_hat, size_t post_initial_num_particles, Model* model, RandomGenerator* rg, vector <bool> data_for_init_particles, bool withdata);

