#include <vector>

#include "scrm/src/model.h"
#include "scrm/src/node.h"
#include "scrm/src/random/mersenne_twister.h"
#include "scrm/src/forest.h"

class ModelSummary{

    public:

    ModelSummary(Model *model, double top_t){
        top_t_scaled = top_t * 4 * model->default_pop_size();
        this->model = model;
        finalized_ = false;

        for( size_t i=0 ; i < model->getNumEpochs() ; i++ ){
            times_.push_back(model->change_times().at(i));
        }

        assert( top_t_scaled > times_.at(times_.size()-1) );
        times_.push_back(top_t_scaled);

        tree_count_ = 0;

        for( size_t idx=0 ; idx < times_.size()-1 ; idx++ ) {
            // size the current tree epoch info vectors
            current_tree_B_below_.push_back( 0 );
            current_tree_B_within_.push_back( 0 );
            current_tree_lineage_count_.push_back( 0 );
            current_tree_single_lineage_.push_back( false );
        }
        for( size_t bias_section_idx=0; bias_section_idx < model->bias_strengths().size(); bias_section_idx++) {
            // size the current tree bias section info vector
            current_tree_B_below_bh_.push_back( 0 );
        }
    };

    ~ModelSummary(){};

    std::vector<double> times_; // times_ contains all time boundaries (specifically starts with 0, ends with top_t_)

    double tree_count_;
    Model* model;
    double top_t_scaled;

    // Need to measure below and within epochs
    std::vector<double> avg_B_below_;
    std::vector<double> avg_B_within_;
    // Need to measure below and within bias sections
    std::vector<double> avg_B_below_bh_;
    std::vector<double> avg_B_within_bias_section_;

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
    std::vector<double> current_tree_B_below_bh_;
    double current_tree_B_;
    std::vector<double> current_tree_lineage_count_;
    std::vector<bool> current_tree_single_lineage_;

    // getters
    std::vector<double> times() const {return times_;}

    std::vector<double> avg_B_below() const {return avg_B_below_;}
    std::vector<double> avg_B_within() const {return avg_B_within_;}
    std::vector<double> avg_B_below_bh() const {return avg_B_below_bh_;}
    std::vector<double> avg_B_within_bias_section() const {return avg_B_within_bias_section_;}
    double avg_B() const {return avg_B_;}
    std::vector<double> avg_lineage_count() const {return avg_lineage_count_;}

    std::vector<int> single_lineage_count() const {return single_lineage_count_;}

    bool is_finalized() const {return finalized_;}

    std::vector<double> current_tree_B_below() const {return current_tree_B_below_;}
    std::vector<double> current_tree_B_within() const {return current_tree_B_within_;}
    std::vector<double> current_tree_B_below_bh() const {return current_tree_B_below_bh_;}
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
    void set_avg_B_below_bh(size_t idx, double new_value ){avg_B_below_bh_.at(idx) = new_value;}
    void calculate_avg_B_within_bias_section();
    void set_avg_B( double new_value ){avg_B_ = new_value;}
    void set_avg_lineage_count(size_t idx, double new_value ){avg_lineage_count_.at(idx) = new_value;}
    void set_exp_lineage_count_given_two(size_t idx, double new_value ){exp_lineage_count_given_two_.at(idx)=new_value;}

    void set_current_tree_B_below(size_t idx, double new_value ){current_tree_B_below_.at(idx) = new_value;}
    void set_current_tree_B_within(size_t idx, double new_value ){current_tree_B_within_.at(idx) = new_value;}
    void set_current_tree_B_below_bh(size_t idx, double new_value ){current_tree_B_below_bh_.at(idx) = new_value;}
    void set_current_tree_B( double new_value ){current_tree_B_ = new_value;}
    void set_current_tree_lineage_count(size_t idx, double new_value ){current_tree_lineage_count_.at(idx) = new_value;}
    void set_current_tree_single_lineage(size_t idx, bool new_value ){current_tree_single_lineage_.at(idx) = new_value;}

    double k_calculation(double, int);
    double lag_calculation(double, double);

    /// functions for once ModelSummary is finalized
    std::vector<double> getLags() const {
        assert( is_finalized() );
        return lags_;
    }

    double getBiasRatio( size_t idx ) {
        assert( is_finalized() );
        assert( avg_B_within_bias_section().size()==model->bias_strengths().size() );

        double normalizing_factor = 0;
        for( size_t temp_idx=0; temp_idx < model->bias_strengths().size(); temp_idx++ ){
            normalizing_factor += avg_B_within_bias_section()[temp_idx] * model->bias_strengths()[temp_idx];
        }
        return avg_B_ * model->bias_strengths()[idx] / normalizing_factor;
    }

};
