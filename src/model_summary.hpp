#include <vector>

#include "scrm/src/model.h"
#include "scrm/src/node.h"
#include "scrm/src/random/mersenne_twister.h"
#include "scrm/src/forest.h"

class ModelSummary{

    public:

    ModelSummary(Model *model, double top_t){

	this->top_t_scaled = top_t * 4 *model->default_pop_size;
	this->model = model;

        for( size_t i=0 ; i < model->getNumEpochs() ; i++ ){
	    times_.push_back(model->change_times().at(i));
	}

	assert( top_t_scaled > times_.at(times_.size()) );
	times_.push_back(top_t_scaled);

	tree_count_ = 0;

	for( size_t idx=0 ; idx < this->times_.size()-1 ; idx++ ) {
          current_tree_B_below_.push_back( 0 );
	  current_tree_B_within_.push_back( 0 );
	  current_tree_lineage_count_.push_back( 0 );
          current_tree_single_lineage_.push_back( false );
	}
	set_current_tree_B_below_bh( 0 );
    };
    ~ModelSummary(){};

    std::vector<double> times_; // times_ contains all time boundaries (specifically starts with 0, ends with top_t_)

    double tree_count_;
    Model* model;
    double top_t_scaled;
    

    std::vector<double> avg_B_below_;
    std::vector<double> avg_B_within_;
    double avg_B_below_bh_;
    double avg_B_;
    std::vector<double> avg_lineage_count_;

    std::vector<int> single_lineage_count_;

    bool finalized_;
    bool current_tree_traversed_;

    std::vector<double> exp_lineage_count_given_two_;
    std::vector<double> lags_;

    // information for the current tree
    std::vector<double> current_tree_B_below_;
    std::vector<double> current_tree_B_within_;
    double current_tree_B_below_bh_;
    double current_tree_B_;
    std::vector<double> current_tree_lineage_count_;
    std::vector<bool> current_tree_single_lineage_;

    // getters
    std::vector<double> times() const {return times_;}

    std::vector<double> avg_B_below() const {return avg_B_below_;}
    std::vector<double> avg_B_within() const {return avg_B_within_;}
    double avg_B_below_bh() const {return avg_B_below_bh_;}
    double avg_B() const {return avg_B_;}
    std::vector<double> avg_lineage_count() const {return avg_lineage_count_;}

    std::vector<int> single_lineage_count() const {return single_lineage_count_;}

    bool is_finalized() const {return finalized_;}

    //need getters and setters for current_tree_*
    std::vector<double> current_tree_B_below() const {return current_tree_B_below_;}
    std::vector<double> current_tree_B_within() const {return current_tree_B_within_;}
    double current_tree_B_below_bh() const {return current_tree_B_below_bh_;}
    double current_tree_B() const {return current_tree_B_;}
    std::vector<double> current_tree_lineage_count() const {return current_tree_lineage_count_;}
    std::vector<bool> current_tree_single_lineage() const {return current_tree_single_lineage_;}
    
    /// functions for before ModelSummary is finalized
    void addTree(); //this will simulate a new Forest, take measurements, and adjust averages
    void adjust_current_tree_measurements(Node* node);
    void process_current_tree_measurements();
    void finalize();

    // setters
    void set_avg_B_below(size_t idx, double new_value ){avg_B_below_.at(idx) = new_value;}
    void set_avg_B_within(size_t idx, double new_value ){avg_B_within_.at(idx) = new_value;}
    void set_avg_B_below_bh( double new_value ){avg_B_below_bh_ = new_value;}
    void set_avg_B( double new_value ){avg_B_ = new_value;}
    void set_avg_lineage_count(size_t idx, double new_value ){avg_lineage_count_.at(idx) = new_value;}
    void set_exp_lineage_count_given_two(size_t idx, double new_value ){exp_lineage_count_given_two_.at(idx)=new_value;}

    void set_current_tree_B_below(size_t idx, double new_value ){current_tree_B_below_.at(idx) = new_value;}
    void set_current_tree_B_within(size_t idx, double new_value ){current_tree_B_within_.at(idx) = new_value;}
    void set_current_tree_B_below_bh( double new_value ){current_tree_B_below_bh_ = new_value;}
    void set_current_tree_B( double new_value ){current_tree_B_ = new_value;}
    void set_current_tree_lineage_count(size_t idx, double new_value ){current_tree_B_below_.at(idx) = new_value;}
    void set_current_tree_single_lineage(size_t idx, bool new_value ){current_tree_B_below_.at(idx) = new_value;}

    double k_calculation(double, int);
    double lag_calculation(double, double);

    /// functions for once ModelSummary is finalized
    std::vector<double> getLags() const {
	assert(is_finalized());
	return lags_;}
    
    double getBiasRatioUpper() {
	assert(is_finalized());
	return (avg_B_)/(model->bias_strength() * avg_B_below_bh_ + (avg_B_ - avg_B_below_bh_));}
    double getBiasRatioLower() {
	assert(is_finalized());
	return model->bias_strength()*avg_B_/(model->bias_strength() * avg_B_below_bh_ + (avg_B_ - avg_B_below_bh_));}

};
