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
    forest.buildInitialTree( false );
    set_current_tree_B(forest.local_root()->length_below());

    current_tree_traversed_ = false;

    // measuring function called on root node
    adjust_current_tree_measurements(forest.local_root());
    current_tree_traversed_ = true;

    // set current tree single lineage called in processing
    process_current_tree_measurements();

    // avg_B_below_ etc modified
    if(tree_count_==1) {
        // size the vectors containing vectors; must be a better way to initialize this...
        for( size_t idx=0 ; idx < this->times_.size()-1 ; idx++){
            avg_B_below_.push_back(current_tree_B_below().at(idx));
            avg_B_within_.push_back(current_tree_B_within().at(idx));
            avg_lineage_count_.push_back(current_tree_lineage_count().at(idx));
            if(current_tree_single_lineage().at(idx)) {single_lineage_count_.push_back(1);} else {single_lineage_count_.push_back(0);}
        }
        set_avg_B(current_tree_B());

        for( size_t idx = 0; idx < current_tree_B_below_bh().size(); idx++) {
            avg_B_below_bh_.push_back( current_tree_B_below_bh()[idx]);
        }
    } else {
        for( size_t idx=0 ; idx < times_.size()-1 ; idx++){
            set_avg_B_below(idx,((tree_count_-1)/tree_count_)*avg_B_below().at(idx) +
                            (1/tree_count_)*current_tree_B_below().at(idx));
            set_avg_B_within(idx,((tree_count_-1)/tree_count_)*avg_B_within().at(idx) +
                            (1/tree_count_)*current_tree_B_within().at(idx));
            set_avg_lineage_count(idx,((tree_count_-1)/tree_count_)*avg_lineage_count().at(idx) +
                            (1/tree_count_)*current_tree_lineage_count().at(idx));
            if(current_tree_single_lineage().at(idx)) {single_lineage_count_.at(idx)++;}
        }
        set_avg_B(((tree_count_-1)/tree_count_)*avg_B() + (1/tree_count_)*current_tree_B());

        for( size_t idx = 0; idx < avg_B_below_bh().size(); idx++) {
            set_avg_B_below_bh(idx, ((tree_count_-1)/tree_count_)*avg_B_below_bh()[idx] + (1/tree_count_)*current_tree_B_below_bh()[idx]);
        }
    }

    delete randomgenerator;

}

void ModelSummary::adjust_current_tree_measurements(Node* node){
    assert( !is_finalized() );
    assert( !current_tree_traversed_ );

    if(node->is_root()) {

        //reset for the new tree
        for( size_t idx=0 ; idx < times_.size()-1 ; idx++ ) {
            set_current_tree_B_below( idx, 0 );
            set_current_tree_B_within( idx, 0 );
            set_current_tree_lineage_count( idx, 0 );
            set_current_tree_single_lineage( idx, false );
        }
        for( size_t idx = 0; idx < avg_B_below_bh().size(); idx++) {
            set_current_tree_B_below_bh( idx, 0 );
        }
    } else {
        //check if current_tree_B_below_bh needs to be updated
        for( size_t bh_idx = 1; bh_idx < model->bias_heights().size(); bh_idx++){
            if( node->height() < model->bias_heights()[bh_idx] ) { //if node.height < bh
                set_current_tree_B_below_bh( bh_idx-1, current_tree_B_below_bh()[bh_idx-1] +
                                std::min( node->parent()->height() - node->height() ,
                                model->bias_heights()[bh_idx] - node->height() ) );
            }
        }
    }

    for ( size_t idx=0 ; idx < times_.size()-1 ; idx++ ) {

        //check if B_below[i] or B_within[i] need to be updated
        if( node->height() < times_.at(idx+1) ) { //if node.height < t[i]
            if( node->is_root() ) { //if node.parent==NULL
                set_current_tree_B_within( idx , current_tree_B_within().at(idx) +
                                    std::min(times_.at(idx+1) - node->height(),
                                    times_.at(idx+1) - times_.at(idx)) );
                //above should be the same as pseudo script below
                //B_within[i] +=min(t[i]-node.height,t[i]-t[i-1])

            } else if ( node->parent()->height() < times_.at(idx+1) ) {//elif parent < t[i]
                set_current_tree_B_below( idx , current_tree_B_below().at(idx) + node->parent()->height() - node->height() );
                set_current_tree_B_within( idx , current_tree_B_within().at(idx) +
                                    std::min( std::min(node->parent()->height() - node->height(),
                                    std::max(node->parent()->height() - times_.at(idx),double(0))) ,
                                    times_.at(idx+1) - times_.at(idx) ));
                //above should be the same as pseudo script below
                //B_below[i] +=parent.height-node.height
                //B_within[i] +=min(parent-node,parent-t[i-1],t[i]-t[i-1])

            } else {
                set_current_tree_B_below( idx , current_tree_B_below().at(idx) +
                                    times_.at(idx+1) - node->height() );
                set_current_tree_B_within( idx , current_tree_B_within().at(idx) +
                                    std::min( times_.at(idx+1) - node->height(),
                                    times_.at(idx+1) - times_.at(idx) ));
            //above should be the same as pseudo script below
            //B_below[i] +=t[i]-node.height
            //B_within[i] +=min(t[i]-node.height,t[i]-t[i-1])

            }
        }
    }

    if( node->first_child() != NULL ) {adjust_current_tree_measurements(node->first_child());}
    if( node->second_child() != NULL ) {adjust_current_tree_measurements(node->second_child());}

    //assertions on B_below, B_within
    assert( current_tree_B_below().size() == times_.size()-1 );
    assert( current_tree_B_within().size() == times_.size()-1 );
}

// below calculates exp num of lineages and bool for single lineage; called AFTER having iterated through nodes
void ModelSummary::process_current_tree_measurements(){
    assert( !is_finalized() );
    assert( current_tree_traversed_ );

    for( size_t idx = 0 ; idx < times_.size()-1 ; idx++ ){
        set_current_tree_lineage_count(idx, current_tree_B_within().at(idx) /
                                    (times_.at(idx+1) - times_.at(idx)) );
        if( idx > 0) {
            // If this epoch has only one lineage, the bottom t has a single lineage going through it
            set_current_tree_single_lineage(idx-1, current_tree_lineage_count().at(idx) == 1 );
        }
    }
    set_current_tree_single_lineage(current_tree_lineage_count_.size()-1, true); //assume the most ancient time had a single lineage
}

void ModelSummary::finalize(){
    assert( !is_finalized() );
    assert( current_tree_traversed_ );

    calculate_avg_B_within_bias_section();

    finalized_ = true;
}

void ModelSummary::calculate_avg_B_within_bias_section(){
    avg_B_within_bias_section_.empty();
    avg_B_within_bias_section_.push_back( avg_B_below_bh()[0] );
    for( size_t idx = 1; idx < avg_B_below_bh().size(); idx++){
        avg_B_within_bias_section_.push_back( avg_B_below_bh()[idx] - avg_B_below_bh()[idx-1]  );
    }
    assert( avg_B_within_bias_section().size() == model->bias_strengths().size() );
}

