#include "tree_measurements.hpp"

void TreeMeasurements::init(Forest* forest) {
    this->set_forest(forest);
}

//Will call this in smcsmc.cpp (on root node) and keep a running average of B_below, etc
//also need Prob k=1 at each t
void TreeMeasurements::adjust_branch_measurements(Node* node) {

    if(node->is_root()) {
	//reset
	for( size_t time_idx=0 ; time_idx < this->forest().model().change_times_.size() ; time_idx++ ) {
          set_B_below( time_idx , 0 );
	  set_B_within( time_idx , 0 );
	}
	set_B_below_bh( 0 );
    }

    for ( size_t time_idx=0 ; time_idx < this->forest().model().change_times_.size() ; time_idx++ ) {
      //check if B_below[i] or B_within[i] need to be updated
      if( node->height() < this->forest().model().change_times_.at(time_idx) ) { //if node.height < t[i]
	if( node->parent() == NULL) { //if node.parent==NULL
	  this->set_B_within( time_idx , this->B_within().at(time_idx) +
	                                std::min(this->forest().model().change_times_.at(time_idx) - node->height(),
					    this->forest().model().change_times_.at(time_idx) - this->forest().model().change_times_.at(time_idx-1)) );
	  //above should be the same as pseudo script below
	  //B_within[i] +=min(t[i]-node.height,t[i]-t[i-1])
	}
	else if ( node->parent()->height() < this->forest().model().change_times_.at(time_idx) ) {//elif parent < t[i]
	  this->set_B_below( time_idx , this->B_below().at(time_idx) + node->parent()->height() - node->height() );
	  this->set_B_within( time_idx , this->B_within().at(time_idx) +
	                                 std::min( std::min(node->parent()->height() - node->height(),
					      node->parent()->height() - this->forest().model().change_times_.at(time_idx-1)) ,
					      this->forest().model().change_times_.at(time_idx) - this->forest().model().change_times_.at(time_idx-1) ));
	  //above should be the same as pseudo script below
	  //B_below[i] +=parent.height-node.height
	  //B_within[i] +=min(parent-node,parent-t[i-1],t[i]-t[i-1])
	}
	else { //else
	  this->set_B_below( time_idx , this->B_below().at(time_idx) +
					this->forest().model().change_times_.at(time_idx) - node->height() );
	  this->set_B_within( time_idx , this->B_within().at(time_idx) +
	                                 std::min( this->forest().model().change_times_.at(time_idx) - node->height(),
					      this->forest().model().change_times_.at(time_idx) - this->forest().model().change_times_.at(time_idx-1) ));
	  //above should be the same as pseudo script below
	  //B_below[i] +=t[i]-node.height
	  //B_within[i] +=min(t[i]-node.height,t[i]-t[i-1])
	}
      }
    }

    //check if B_below_bh needs to be updated
    if( node->height() < this->forest().model().bias_height() ) { //if node.height < bh
      this->set_B_below_bh( this->B_below_bh() +
                               std::min( node->parent()->height() - node->height() ,
                                  this->forest().model().bias_height() - node->height() ) );
    }

    if( node->first_child() != NULL ) {this->adjust_branch_measurements(node->first_child());}
    if( node->second_child() != NULL ) {this->adjust_branch_measurements(node->second_child());}

    //assertions on B_below, B_within, B_below_bh
    assert( this->B_below().size() == this->forest().model().change_times_.size() );
    assert( this->B_within().size() == this->forest().model().change_times_.size() );
}

void TreeMeasurements::count_lineage() {
    //assertion to check B_within is appropriate

    for( size_t time_idx = 0 ; time_idx < this->forest().model().change_times_.size() ; time_idx++ ){
      this->set_lineage_count( time_idx, this->B_within().at(time_idx) /
          (this->forest().model().change_times_.at(time_idx) - this->forest().model().change_times_.at(time_idx-1)) );
      this->set_single_lineage( time_idx , (this->lineage_count(time_idx) == 1) );
    }

    //assertions to check lineage_count and single_lineage are appropriate
}

