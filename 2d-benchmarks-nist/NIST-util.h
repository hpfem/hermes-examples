#ifndef __NIST_UTIL_H
#define __NIST_UTIL_H

#include "hermes2d.h"

using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;
using namespace Hermes::Hermes2D::RefinementSelectors;

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
};

bool adaptive_step_single_space(
  MeshSharedPtr& mesh,
  SpaceSharedPtr<double>& space, 
  MeshFunctionSharedPtr<double>& sln, 
  MySelector& selector,
  MeshFunctionSharedPtr<double>& ref_sln,
  Hermes::Mixins::TimeMeasurable& cpu_time,
  Solver<double>& solver,
  Views::ScalarView& sview,
  Views::OrderView & oview,
  SimpleGraph graph_dof_est,
  SimpleGraph graph_cpu_est,
  ErrorCalculator<double>& error_calculator,
  Adapt<double>& adaptivity,
  int& as,
  double ERR_STOP,
  MeshFunctionSharedPtr<double>& exact_sln = MeshFunctionSharedPtr<double>(),
  SimpleGraph graph_dof_exact = SimpleGraph(),
  SimpleGraph graph_cpu_exact = SimpleGraph()
  );
#endif