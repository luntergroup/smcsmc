#include "model_summary.hpp"
#include <cmath>
using namespace std;

void ModelSummary::addTree(){
    assert( !is_finalized() );
    
    tree_count_++;
    
    // simulate tree using tree_count_ as seed
    MersenneTwister *randomgenerator = new MersenneTwister( true, tree_count_ );
    //Forest* forest = new Forest( this->model, randomgenerator );
    Forest forest = Forest( this->model, randomgenerator );
    forest.buildInitialTree();
    set_current_tree_B(forest.local_root()->length_below());

    current_tree_traversed_ = false;
    
    // measuring function called on root node
    this->adjust_current_tree_measurements(forest.local_root());
    current_tree_traversed_ = true;
    //cout << "current tree traversed" << endl;
    //cout << "current tree B below " << current_tree_B_below() << endl;
    // set current tree single lineage called in processing
    this->process_current_tree_measurements();
    //cout << "current tree processed" << endl;
    //cout << "current tree B below " << current_tree_B_below() << endl;
    //cout << "current tree B within " << current_tree_B_within() << endl;
    //cout << "current tree lineage count " << current_tree_lineage_count() << endl;
    //cout << "current tree single lineage " << current_tree_single_lineage() << endl;
    //cout << "current tree B below bh " << current_tree_B_below_bh() << endl;
    
    // avg_B_below_ etc modified
    if(tree_count_==1) {  // could remove this if statement by initializing avgs to 0
      for( size_t idx=0 ; idx < this->times_.size()-1 ; idx++){
        avg_B_below_.push_back(current_tree_B_below().at(idx));
	avg_B_within_.push_back(current_tree_B_within().at(idx));
	avg_lineage_count_.push_back(current_tree_lineage_count().at(idx));
	if(current_tree_single_lineage().at(idx)) {single_lineage_count_.push_back(1);} else {single_lineage_count_.push_back(0);}
      }
      set_avg_B(current_tree_B());
      //cout << "After the first tree avg_B is " << avg_B() << endl;
      //cout << "After the first tree avg_B_below is " << avg_B_below() << endl;
      //cout << "After the first tree avg_B_within is " << avg_B_within() << endl;
      //cout << "After the first tree avg_lineage_count is " << avg_lineage_count() << endl;
      //cout << "After the first tree single_lineage_count is " << single_lineage_count() << endl;
      set_avg_B_below_bh(current_tree_B_below_bh());
    } else {
      for( size_t idx=0 ; idx < this->times_.size()-1 ; idx++){
        set_avg_B_below(idx,((tree_count_-1)/tree_count_)*avg_B_below().at(idx) +
                             (1/tree_count_)*current_tree_B_below().at(idx));
	set_avg_B_within(idx,((tree_count_-1)/tree_count_)*avg_B_within().at(idx) +
                              (1/tree_count_)*current_tree_B_within().at(idx));
	set_avg_lineage_count(idx,((tree_count_-1)/tree_count_)*avg_lineage_count().at(idx) +
                                   (1/tree_count_)*current_tree_lineage_count().at(idx));
	if(current_tree_single_lineage().at(idx)) {single_lineage_count_.at(idx)++;}
      }
      set_avg_B(((tree_count_-1)/tree_count_)*avg_B() + (1/tree_count_)*current_tree_B());
      //cout << "After the next tree avg_B is " << avg_B() << endl;
      //cout << "After the next tree avg_B_below is " << avg_B_below() << endl;
      //cout << "After the next tree avg_B_within is " << avg_B_within() << endl;
      //cout << "After the next tree avg_lineage_count is " << avg_lineage_count() << endl;
      //cout << "After the next tree single_lineage_count is " << single_lineage_count() << endl;
      //cout << "tree count is " <<  tree_count_ << endl;
      set_avg_B_below_bh(((tree_count_-1)/tree_count_)*avg_B_below_bh() + (1/tree_count_)*current_tree_B_below_bh());
    }
}

void ModelSummary::adjust_current_tree_measurements(Node* node){
    assert( !is_finalized() );
    assert( !current_tree_traversed_ );
    //cout << "adjust current tree measurements called on node" << endl;
    //cout << "node is " << node << endl;
    //cout << "node height is " << node->height() << endl;
    
    if(node->is_root()) {
	//cout << "this node is the root" << endl;
	//reset for the new tree
	for( size_t idx=0 ; idx < this->times_.size()-1 ; idx++ ) {
          set_current_tree_B_below( idx, 0 );
	  set_current_tree_B_within( idx, 0 );
	  set_current_tree_lineage_count( idx, 0 );
          set_current_tree_single_lineage( idx, false );
	}
	set_current_tree_B_below_bh( 0 );
    } else {
      //check if current_tree_B_below_bh needs to be updated
      if( node->height() < this->model->bias_height() ) { //if node.height < bh
        this->set_current_tree_B_below_bh( this->current_tree_B_below_bh() +
                               std::min( node->parent()->height() - node->height() ,
                                  this->model->bias_height() - node->height() ) );
      }
    }

    for ( size_t idx=0 ; idx < this->times_.size()-1 ; idx++ ) {
    
      //check if B_below[i] or B_within[i] need to be updated
      if( node->height() < this->times_.at(idx+1) ) { //if node.height < t[i]
	if( node->is_root() ) { //if node.parent==NULL
	  this->set_current_tree_B_within( idx , this->current_tree_B_within().at(idx) +
	                                std::min(this->times_.at(idx+1) - node->height(),
					    this->times_.at(idx+1) - this->times_.at(idx)) );
	  //above should be the same as pseudo script below
	  //B_within[i] +=min(t[i]-node.height,t[i]-t[i-1])
	  //cout << "Case 1: parent==NULL" << endl;
	  //cout << "current tree B below: " << current_tree_B_below() << endl;
	  //cout << "current tree B within: " << current_tree_B_within() << endl;
	}
	else if ( node->parent()->height() < this->times_.at(idx+1) ) {//elif parent < t[i]
	  this->set_current_tree_B_below( idx , this->current_tree_B_below().at(idx) + node->parent()->height() - node->height() );
	  this->set_current_tree_B_within( idx , this->current_tree_B_within().at(idx) +
	                                 std::min( std::min(node->parent()->height() - node->height(),
					      std::max(node->parent()->height() - this->times_.at(idx),double(0))) ,
					      this->times_.at(idx+1) - this->times_.at(idx) ));
	  //above should be the same as pseudo script below
	  //B_below[i] +=parent.height-node.height
	  //B_within[i] +=min(parent-node,parent-t[i-1],t[i]-t[i-1])
	  //cout << "Case 2: parent<t" << endl;
	  //cout << "current tree B below: " << current_tree_B_below() << endl;
	  //cout << "current tree B within: " << current_tree_B_within() << endl;
	}
	else { //else
	  this->set_current_tree_B_below( idx , this->current_tree_B_below().at(idx) +
					this->times_.at(idx+1) - node->height() );
	  this->set_current_tree_B_within( idx , this->current_tree_B_within().at(idx) +
	                                 std::min( this->times_.at(idx+1) - node->height(),
					      this->times_.at(idx+1) - this->times_.at(idx) ));
	  //above should be the same as pseudo script below
	  //B_below[i] +=t[i]-node.height
	  //B_within[i] +=min(t[i]-node.height,t[i]-t[i-1])
	  //cout << "Case 3: parent>t" << endl;
	  //cout << "current tree B below: " << current_tree_B_below() << endl;
	  //cout << "current tree B within: " << current_tree_B_within() << endl;
	}
      }
    }

    if( node->first_child() != NULL ) {this->adjust_current_tree_measurements(node->first_child());} 
    if( node->second_child() != NULL ) {this->adjust_current_tree_measurements(node->second_child());}

    //assertions on B_below, B_within, B_below_bh
    assert( this->current_tree_B_below().size() == this->times_.size()-1 );
    assert( this->current_tree_B_within().size() == this->times_.size()-1 );
}

// below calculates exp num of lineages and bool for single lineage; called AFTER having iterated through nodes
void ModelSummary::process_current_tree_measurements(){
    assert( !is_finalized() );
    assert( current_tree_traversed_ );

    for( size_t idx = 0 ; idx < this->times_.size()-1 ; idx++ ){
      set_current_tree_lineage_count(idx, this->current_tree_B_within().at(idx) /
          (this->times_.at(idx+1) - this->times_.at(idx)) );
      if( idx > 0) {
	// If this epoch has only one lineage, the bottom t has a single lineage going through it
        set_current_tree_single_lineage(idx-1, this->current_tree_lineage_count().at(idx) == 1 );
      }
    }
    set_current_tree_single_lineage(current_tree_lineage_count_.size()-1, true); //assume the most ancient time had a single lineage
}

void ModelSummary::finalize(){
    assert( !is_finalized() );
    assert( current_tree_traversed_ );

    for(size_t idx = 0 ; idx < this->times_.size()-1 ; idx++ ){
      // calculate Exp_branch_count_given_two
      exp_lineage_count_given_two_.push_back( k_calculation(avg_lineage_count().at(idx),single_lineage_count().at(idx)) );
      // calculate lags
      if(model->sample_size()==2) {
        lags_.push_back( 1 / ( model->recombination_rate()*2*this->times_.at(idx+1) ) );
      } else {
        lags_.push_back( lag_calculation(avg_B_below().at(idx),exp_lineage_count_given_two_.at(idx)) );
      }
    }

    assert(exp_lineage_count_given_two_.size()==times_.size()-1);
    assert(lags_.size()==times_.size()-1);
    finalized_ = true;
}

double ModelSummary::k_calculation(double exp_lineage_count, int single_lineage_count){

    return std::min((exp_lineage_count - double(single_lineage_count)/tree_count_)/
                    (1 - double(single_lineage_count)/tree_count_) , double(2));
}

double ModelSummary::lag_calculation(double B_below, double k){

    return 1 / ( model->recombination_rate() * B_below *
             ( (2/k)*(k-2)/(model->sample_size()-2) + .5*(model->sample_size()-k)/(model->sample_size()-2) ));

}
