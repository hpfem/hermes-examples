#ifndef __NIST_UTIL_H
#define __NIST_UTIL_H

#include "hermes2d.h"

using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;
using namespace Hermes::Hermes2D::RefinementSelectors;

extern const char* thresholds[7];
extern const double threshold_values[7];

enum hpAdaptivityStrategy
{
  noSelectionH = 0,
  noSelectionHP = 1,
  hXORpSelectionBasedOnError = 2,
  hORpSelectionBasedOnDOFs = 3,
  isoHPSelectionBasedOnDOFs = 4,
  anisoHPSelectionBasedOnDOFs = 5
};
extern const char* strategies[6];

class MySelector : public H1ProjBasedSelector<double>
{
public:
  MySelector(hpAdaptivityStrategy strategy) : H1ProjBasedSelector<double>(cand_list), strategy(strategy)
  {
  }
private:
  bool select_refinement(Element* element, int order, MeshFunction<double>* rsln, ElementToRefine& refinement)
  {
    switch(strategy)
    {
    case(noSelectionH):
      {
        refinement.split = H2D_REFINEMENT_H;
        refinement.best_refinement_polynomial_order_type[H2D_REFINEMENT_H][0] = 
          refinement.best_refinement_polynomial_order_type[H2D_REFINEMENT_H][1] = 
          refinement.best_refinement_polynomial_order_type[H2D_REFINEMENT_H][2] = 
          refinement.best_refinement_polynomial_order_type[H2D_REFINEMENT_H][3] = 
          order;
        ElementToRefine::copy_orders(refinement.refinement_polynomial_order, refinement.best_refinement_polynomial_order_type[H2D_REFINEMENT_H]);
        return true;
      }
      break;
    case(noSelectionHP):
      {
        int max_allowed_order = this->max_order;
        if(this->max_order == H2DRS_DEFAULT_ORDER)
          max_allowed_order = H2DRS_MAX_ORDER;
        int order_h = H2D_GET_H_ORDER(order), order_v = H2D_GET_V_ORDER(order);
        int increased_order_h = std::min(max_allowed_order, order_h + 1), increased_order_v = std::min(max_allowed_order, order_v + 1);
        int increased_order;
        if(element->is_triangle())
          increased_order = refinement.best_refinement_polynomial_order_type[H2D_REFINEMENT_H][0] = H2D_MAKE_QUAD_ORDER(increased_order_h, increased_order_h);
        else
          increased_order = refinement.best_refinement_polynomial_order_type[H2D_REFINEMENT_H][0] = H2D_MAKE_QUAD_ORDER(increased_order_h, increased_order_v);

        refinement.split = H2D_REFINEMENT_H;
        refinement.best_refinement_polynomial_order_type[H2D_REFINEMENT_H][0] = 
          refinement.best_refinement_polynomial_order_type[H2D_REFINEMENT_H][1] = 
          refinement.best_refinement_polynomial_order_type[H2D_REFINEMENT_H][2] = 
          refinement.best_refinement_polynomial_order_type[H2D_REFINEMENT_H][3] = 
          increased_order;
        ElementToRefine::copy_orders(refinement.refinement_polynomial_order, refinement.best_refinement_polynomial_order_type[H2D_REFINEMENT_H]);
        return true;
      }
      case(hXORpSelectionBasedOnError):
      {
        //make an uniform order in a case of a triangle
        int order_h = H2D_GET_H_ORDER(order), order_v = H2D_GET_V_ORDER(order);

        int current_min_order, current_max_order;
        this->get_current_order_range(element, current_min_order, current_max_order);

        if(current_max_order < std::max(order_h, order_v))
          current_max_order = std::max(order_h, order_v);

        int last_order_h = std::min(current_max_order, order_h + 1), last_order_v = std::min(current_max_order, order_v + 1);
        int last_order = H2D_MAKE_QUAD_ORDER(last_order_h, last_order_v);

        //build candidates.
        Hermes::vector<Cand> candidates;
        candidates.push_back(Cand(H2D_REFINEMENT_P, last_order));
        candidates.push_back(Cand(H2D_REFINEMENT_H, order, order, order, order));
        Cand* best_candidate = (candidates[0].error > candidates[1].error) ? &candidates[0] : &candidates[1];
        Cand* best_candidates_specific_type[4];
        best_candidates_specific_type[H2D_REFINEMENT_P] = &candidates[0];
        best_candidates_specific_type[H2D_REFINEMENT_H] = &candidates[1];
        best_candidates_specific_type[2] = NULL;
        best_candidates_specific_type[3] = NULL;

        //copy result to output
        refinement.split = best_candidate->split;
        ElementToRefine::copy_orders(refinement.refinement_polynomial_order, best_candidate->p);
        for(int i = 0; i < 4; i++)
          if(best_candidates_specific_type[i] != NULL)
            ElementToRefine::copy_orders(refinement.best_refinement_polynomial_order_type[i], best_candidates_specific_type[i]->p);

        ElementToRefine::copy_errors(refinement.errors, best_candidate->errors);

        //modify orders in a case of a triangle such that order_v is zero
        if(element->is_triangle())
          for(int i = 0; i < H2D_MAX_ELEMENT_SONS; i++)
            refinement.refinement_polynomial_order[i] = H2D_MAKE_QUAD_ORDER(H2D_GET_H_ORDER(refinement.refinement_polynomial_order[i]), 0);

        return true;
      }
    default:
      H1ProjBasedSelector<double>::select_refinement(element, order, rsln, refinement);
      return true;
      break;
    }
  }
  
  Hermes::vector<Cand> create_candidates(Element* e, int quad_order)
  {
    Hermes::vector<Cand> candidates;

    // Get the current order range.
    int current_min_order, current_max_order;
    this->get_current_order_range(e, current_min_order, current_max_order);

    int order_h = H2D_GET_H_ORDER(quad_order), order_v = H2D_GET_V_ORDER(quad_order);

    if(current_max_order < std::max(order_h, order_v))
      current_max_order = std::max(order_h, order_v);

    int last_order_h = std::min(current_max_order, order_h + 1), last_order_v = std::min(current_max_order, order_v + 1);
    int last_order = H2D_MAKE_QUAD_ORDER(last_order_h, last_order_v);

    switch(strategy)
    {
    case(hORpSelectionBasedOnDOFs):
      {
        candidates.push_back(Cand(H2D_REFINEMENT_P, quad_order));
      }
    case(hXORpSelectionBasedOnError):
      {
        candidates.push_back(Cand(H2D_REFINEMENT_P, last_order));
        candidates.push_back(Cand(H2D_REFINEMENT_H, quad_order, quad_order, quad_order, quad_order));
        return candidates;
      }
      break;
    case(isoHPSelectionBasedOnDOFs):
      {
        this->cand_list = H2D_HP_ISO;
        return H1ProjBasedSelector<double>::create_candidates(e, quad_order);
      }
      break;
    case(anisoHPSelectionBasedOnDOFs):
      {
        this->cand_list = H2D_HP_ANISO;
        return H1ProjBasedSelector<double>::create_candidates(e, quad_order);
      }
      break;
    }
  }

  void evaluate_cands_score(Hermes::vector<Cand>& candidates, Element* e)
  {
    switch(strategy)
    {
    case(hXORpSelectionBasedOnError):
      {
        if(candidates[0].error > candidates[1].error)
        {
          candidates[0].score = 0.0;
          candidates[1].score = 1.0;
        }
        else
        {
          candidates[1].score = 0.0;
          candidates[0].score = 1.0;
        }
      }
      break;
    default:
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
    }
  }

  int strategy;
};

char* process_arguments_main_comparison(int argc, char* argv[], Selector<double>*& selector, AdaptivityStoppingCriterion<double>*& stoppingCriterion);

bool adaptive_step_single_space(
  Hermes::Mixins::Loggable* logger,
  MeshSharedPtr& mesh,
  SpaceSharedPtr<double>& space, 
  MeshFunctionSharedPtr<double>& sln, 
  Selector<double>* selector,
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
  double& FactorizationSize,
  double& PeakMemoryUsage,
  double& Flops,
  MeshFunctionSharedPtr<double>& exact_sln = MeshFunctionSharedPtr<double>()
  );
#endif