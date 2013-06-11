#include "definitions.h"

//  This problem describes the distribution of the vector potential in
//  a 2D domain comprising a wire carrying electrical current, air, and
//  an iron which is not under voltage.
//
//  PDE: -div(1/rho grad p) - omega**2 / (rho c**2) * p = 0.
//
//  Domain: Floor plan of an existing apartment.
//
//  BC: Prescribed pressure at the source.
//      Newton matched boundary on the walls.
//
//  The following parameters can be changed:

// Number of initial uniform mesh refinements.
const int INIT_REF_NUM = 0;                       
// Initial polynomial degree of mesh elements.
const int P_INIT = 2;                             
// This is a quantitative parameter of the adapt(...) function and
// it has different meanings for various adaptive strategies.
const double THRESHOLD = 0.25;                     

// Predefined list of element refinement candidates.
const CandList CAND_LIST = H2D_HP_ANISO;
// Stopping criterion for adaptivity.
const double ERR_STOP = 1e-1;

// Problem parameters.
const double RHO = 1.25;
const double FREQ = 5e2;
const double OMEGA = 2 * M_PI * FREQ;
const double SOUND_SPEED = 353.0;
const std::complex<double> P_SOURCE(1.0, 0.0);

enum hpAdaptivityStrategy
{
  noSelectionH = 0,
  noSelectionHP = 1,
  hXORpSelectionBasedOnError = 2,
  hORpSelectionBasedOnDOFs = 3,
  isoHPSelectionBasedOnDOFs = 4,
  anisoHPSelectionBasedOnDOFs = 5
};

class MySelector : public H1ProjBasedSelector<complex>
{
public:
  MySelector(hpAdaptivityStrategy strategy) : H1ProjBasedSelector<complex>(cand_list), strategy(strategy)
  {
    if(strategy == hXORpSelectionBasedOnError)
    {
      //this->set_error_weights(1.0,1.0,1.0);
    }
  }
private:
  bool select_refinement(Element* element, int order, MeshFunction<complex>* rsln, ElementToRefine& refinement)
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
        
        this->evaluate_cands_error(candidates, element, rsln);
        
        Cand* best_candidate = (candidates[0].error < candidates[1].error) ? &candidates[0] : &candidates[1];
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
      H1ProjBasedSelector<complex>::select_refinement(element, order, rsln, refinement);
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
        return H1ProjBasedSelector<complex>::create_candidates(e, quad_order);
      }
      break;
    case(anisoHPSelectionBasedOnDOFs):
      {
        this->cand_list = H2D_HP_ANISO;
        return H1ProjBasedSelector<complex>::create_candidates(e, quad_order);
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

int main(int argc, char* argv[])
{
  // Initialize refinement selector.
  MySelector selector(hORpSelectionBasedOnDOFs);

  HermesCommonApi.set_integral_param_value(Hermes::showInternalWarnings, false);

  // Time measurement.
  Hermes::Mixins::TimeMeasurable cpu_time;
  cpu_time.tick();

  // Load the mesh.
  MeshSharedPtr mesh(new Mesh);
  MeshReaderH2D mloader;
  mloader.load("domain.mesh", mesh);

  // Error calculation & adaptivity.
  DefaultErrorCalculator<complex, HERMES_H1_NORM> errorCalculator(RelativeErrorToGlobalNorm, 1);
  // Stopping criterion for an adaptivity step.
  AdaptStoppingCriterionSingleElement<complex> stoppingCriterion(THRESHOLD);
  // Adaptivity processor class.
  Adapt<complex> adaptivity(&errorCalculator, &stoppingCriterion);

  // Perform initial mesh refinements.
  for (int i = 0; i < INIT_REF_NUM; i++) mesh->refine_all_elements();

  // Initialize boundary conditions.
  DefaultEssentialBCConst<complex> bc_essential("Source", P_SOURCE);
  EssentialBCs<complex> bcs(&bc_essential);

  // Create an H1 space with default shapeset.
  SpaceSharedPtr<complex> space(new H1Space<complex> (mesh, &bcs, P_INIT));
  int ndof = Space<complex>::get_num_dofs(space);
  //Hermes::Mixins::Loggable::Static::info("ndof = %d", ndof);

  // Initialize the weak formulation.
  CustomWeakFormAcoustics wf("Wall", RHO, SOUND_SPEED, OMEGA);

  // Initialize coarse and reference mesh solution.
  MeshFunctionSharedPtr<complex>  sln(new Solution<complex>), ref_sln(new Solution<complex>);

  // Initialize views.
  ScalarView sview("Acoustic pressure", new WinGeom(600, 0, 600, 350));
  sview.show_contours(.2);
  ScalarView eview("Error", new WinGeom(600, 377, 600, 350));
  sview.show_mesh(false);
  sview.fix_scale_width(50);
  OrderView  oview("Polynomial orders", new WinGeom(1208, 0, 600, 350));
  ScalarView ref_view("Refined elements", new WinGeom(1208, 377, 600, 350));
  ref_view.show_scale(false);
  ref_view.set_min_max_range(0, 2);

  // DOF and CPU convergence graphs initialization.
  SimpleGraph graph_dof, graph_cpu;

  //Hermes::Mixins::Loggable::Static::info("Solving on reference mesh.");
  // Time measurement.
  cpu_time.tick();

  // Perform Newton's iteration.
  Hermes::Hermes2D::NewtonSolver<complex> newton(&wf, space);
  newton.set_verbose_output(false);

  // Adaptivity loop:
  int as = 1;
  adaptivity.set_space(space);
  bool done = false;
  do
  {
    //Hermes::Mixins::Loggable::Static::info("---- Adaptivity step %d:", as);

    // Construct globally refined reference mesh and setup reference space.
    Mesh::ReferenceMeshCreator refMeshCreator(mesh);
    MeshSharedPtr ref_mesh = refMeshCreator.create_ref_mesh();

    Space<complex>::ReferenceSpaceCreator refSpaceCreator(space, ref_mesh);
    SpaceSharedPtr<complex> ref_space = refSpaceCreator.create_ref_space();

    int ndof_ref = Space<complex>::get_num_dofs(ref_space);
    wf.set_verbose_output(false);
    newton.set_space(ref_space);

    // Assemble the reference problem.
    try
    {
      newton.solve();
    }
    catch(Hermes::Exceptions::Exception e)
    {
      e.print_msg();
      throw Hermes::Exceptions::Exception("Newton's iteration failed.");
    };
    // Translate the resulting coefficient vector into the Solution<complex> sln->
    Hermes::Hermes2D::Solution<complex>::vector_to_solution(newton.get_sln_vector(), ref_space, ref_sln);

    // Project the fine mesh solution onto the coarse mesh.
    //Hermes::Mixins::Loggable::Static::info("Projecting reference solution on coarse mesh.");
    OGProjection<complex> ogProjection; ogProjection.project_global(space, ref_sln, sln);

    // Time measurement.
    cpu_time.tick();

    // View the coarse mesh solution and polynomial orders.
    MeshFunctionSharedPtr<double> acoustic_pressure(new RealFilter(sln));
    sview.show(acoustic_pressure);
    oview.show(space);

    // Calculate element errors and total error estimate.
    //Hermes::Mixins::Loggable::Static::info("Calculating error estimate.");
    errorCalculator.calculate_errors(sln, ref_sln);
    double err_est_rel = errorCalculator.get_total_error_squared() * 100;

    eview.show(errorCalculator.get_errorMeshFunction());

    // Report results.
    //Hermes::Mixins::Loggable::Static::info("ndof_coarse: %d, ndof_fine: %d, err_est_rel: %g%%", Space<complex>::get_num_dofs(space), Space<complex>::get_num_dofs(ref_space), err_est_rel);

    // Time measurement.
    cpu_time.tick();

    // Add entry to DOF and CPU convergence graphs.
    graph_dof.add_values(Space<complex>::get_num_dofs(space), err_est_rel);
    graph_dof.save("conv_dof_est.dat");
    graph_cpu.add_values(cpu_time.accumulated(), err_est_rel);
    graph_cpu.save("conv_cpu_est.dat");

    // If err_est too large, adapt the mesh.
    if (err_est_rel < ERR_STOP) done = true;
    else
    {
      //Hermes::Mixins::Loggable::Static::info("Adapting coarse mesh.");
      done = adaptivity.adapt(&selector);
      ref_view.show(adaptivity.get_refinementInfoMeshFunction());
      cpu_time.tick();
      std::cout << "Adaptivity step: " << as << ", running CPU time: " << cpu_time.accumulated_str() << std::endl;
    }

    // Increase counter.
    as++;
  }
  while (done == false);

  //Hermes::Mixins::Loggable::Static::info("Total running time: %g s", cpu_time.accumulated());

  // Show the reference solution - the final result.
  sview.set_title("Fine mesh solution magnitude");

  MeshFunctionSharedPtr<double>  ref_mag(new RealFilter(ref_sln));
  sview.show(ref_mag);

  // Output solution in VTK format.
  Linearizer lin;
  bool mode_3D = true;
  lin.save_solution_vtk(ref_mag, "sln.vtk", "Acoustic pressure", mode_3D);
  //Hermes::Mixins::Loggable::Static::info("Solution in VTK format saved to file %s.", "sln.vtk");

  // Wait for all views to be closed.
  View::wait();
  return 0;
}
