#include <vector>

#include "scrm/src/forest.h"

class Tree_measurements {

  public:
  
  Tree_measurements( Forest *forest );
  ~Tree_measurements();

  Forest* forest_;
  void set_forest(Forest* forest) { this->forest_ = forest; }
  const Forest &forest() const { return *forest_; }

  std::vector<double> B_below_;
  std::vector<double> B_within_;
  std::vector<double> lineage_count_;
  std::vector<bool> single_lineage_;
  double B_below_bh_;

  // getters
  std::vector<double> B_below() const { return this->B_below_; }
  std::vector<double> B_within() const { return this->B_within_; }
  std::vector<double> lineage_count() const {return this->lineage_count_; }
  std::vector<bool> single_lineage() const { return this->single_lineage_; }
  double B_below_bh() const { return this->B_below_bh_; }
  // unecessary as can use B_below().at(i), but could be cleaner?
  double B_below( size_t idx ) const { return this->B_below_.at(idx); }
  double B_within( size_t idx ) const { return this->B_within_.at(idx); }
  double lineage_count( size_t idx ) const { return this->lineage_count_.at(idx); }
  bool single_lineage( size_t idx ) const { return this->single_lineage_.at(idx); }

  void init(Forest *forest);

  // adjustors
  void set_B_below( size_t time_idx , double new_value ) { B_below_.at(time_idx) = new_value; }
  void set_B_within( size_t time_idx , double new_value ) { B_within_.at(time_idx) = new_value; }
  void set_lineage_count( size_t time_idx, double new_value ) { lineage_count_.at(time_idx) = new_value; }
  void set_single_lineage( size_t time_idx, bool b) { single_lineage_.at(time_idx) = b; }
  void set_B_below_bh( double new_value ) { B_below_bh_ = new_value; }

  void adjust_branch_measurements(Node* node);
  void count_lineage();

};
