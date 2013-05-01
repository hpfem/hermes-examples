#ifndef __NIST_UTIL_H
#define __NIST_UTIL_H

#include "hermes2d.h"

using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;
using namespace Hermes::Hermes2D::RefinementSelectors;

extern const char* thresholds[4];
extern const double threshold_values[4];

class MySelector : public H1ProjBasedSelector<double>
{
public:
  MySelector(CandList cand_list) : H1ProjBasedSelector<double>(cand_list)
  {
  }
private:
  void evaluate_cands_score(Hermes::vector<Cand>& candidates, Element* e)
  {
    //calculate score of candidates
    Cand& unrefined = candidates[0];
    const int num_cands = (int)candidates.size();
    unrefined.score = 0;

    for (int i = 1; i < num_cands; i++)
    {
      Cand& cand = candidates[i];
      if(cand.error < unrefined.error)
      {
        double delta_dof = cand.dofs - unrefined.dofs;
        candidates[i].score = (log10(unrefined.error) - log10(cand.error)) / delta_dof;
      }
      else
        candidates[i].score = 0;
    }
  }

  Hermes::vector<Cand> create_candidates(Element* e, int quad_order)
  {
    Hermes::vector<Cand> candidates;

    if(this->cand_list == H2D_NONE)
    {

      // Get the current order range.
      int current_min_order, current_max_order;
      this->get_current_order_range(e, current_min_order, current_max_order);

      int order_h = H2D_GET_H_ORDER(quad_order), order_v = H2D_GET_V_ORDER(quad_order);

      if(current_max_order < std::max(order_h, order_v))
        current_max_order = std::max(order_h, order_v);

      int last_order_h = std::min(current_max_order, order_h + 1), last_order_v = std::min(current_max_order, order_v + 1);
      int last_order = H2D_MAKE_QUAD_ORDER(last_order_h, last_order_v);

      candidates.push_back(Cand(H2D_REFINEMENT_P, quad_order));
      candidates.push_back(Cand(H2D_REFINEMENT_P, last_order));
      candidates.push_back(Cand(H2D_REFINEMENT_H, quad_order, quad_order, quad_order, quad_order));
      
      return candidates;
    }
    else
    {
      return H1ProjBasedSelector<double>::create_candidates(e, quad_order);
    }
  }
};

bool adaptive_step_single_space(
  Hermes::Mixins::Loggable* logger,
  MeshSharedPtr& mesh,
  SpaceSharedPtr<double>& space, 
  MeshFunctionSharedPtr<double>& sln, 
  Selector<double>& selector,
  unsigned int order_increase,
  MeshFunctionSharedPtr<double>& ref_sln,
  Hermes::Mixins::TimeMeasurable& cpu_time,
  Solver<double>& solver,
  Views::ScalarView& sview,
  Views::OrderView & oview,
  ErrorCalculator<double>& error_calculator,
  Adapt<double>& adaptivity,
  int& as,
  double error_stop,
  double& error_reached,
  int& dof_reached,
  int& cache_searches,
  int& cache_record_found,
  int& cache_record_found_reinit,
  int& cache_record_not_found,
  double& exact_error_reached,
  MeshFunctionSharedPtr<double>& exact_sln = MeshFunctionSharedPtr<double>()
  );
#endif