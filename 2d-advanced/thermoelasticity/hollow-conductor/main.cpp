#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "definitions.h"

// A massive hollow conductor is heated by induction and cooled by water running inside.
// We will model this problem using linear thermoelasticity equations, where the x-displacement,
// y-displacement, and the temperature will be approximated on individual meshes equipped
// with mutually independent adaptivity mechanisms. Use MULTI = true to use multimesh,
// MULTI = false for single-mesh (all solution components on the samemesh).
//
// PDE: Linear thermoelasticity.
//
// BC: u_1 = u_2 = 0 on Gamma_1
//     du_1/dn = du_2/dn = 0 elsewhere
//     temp = TEMP_INNER on Gamma_4
//     negative heat flux with HEAT_FLUX_OUTER elsewhere.

// Initial polynomial degrees in temperature mesh.
const int P_INIT_TEMP = 1;                        
// Initial polynomial degrees for displacement meshes.
const int P_INIT_DISP = 1;                        
// MULTI = true  ... use multi-mesh,
// MULTI = false ... use single-mesh.
// Note: In the single mesh option, the meshes are
// forced to be geometrically the same but the
// polynomial degrees can still vary.
const bool MULTI = true;                          
// This is a quantitative parameter of the adapt(...) function and
// it has different meanings for various adaptive strategies.
const double THRESHOLD = 0.3;                     
// Adaptive strategy:
// STRATEGY = 0 ... refine elements until sqrt(THRESHOLD) times total
//   error is processed. If more elements have similar errors, refine
//   all to keep the mesh symmetric.
// STRATEGY = 1 ... refine all elements whose error is larger
//   than THRESHOLD times maximum element error.
// STRATEGY = 2 ... refine all elements whose error is larger
//   than THRESHOLD.
const int STRATEGY = 1;                           
// Predefined list of element refinement candidates. Possible values are
// H2D_P_ISO, H2D_P_ANISO, H2D_H_ISO, H2D_H_ANISO, H2D_HP_ISO,
// H2D_HP_ANISO_H, H2D_HP_ANISO_P, H2D_HP_ANISO.
const CandList CAND_LIST = H2D_HP_ANISO;          
// Maximum allowed level of hanging nodes:
// MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
// MESH_REGULARITY = 1 ... at most one-level hanging nodes,
// MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
// Note that regular meshes are not supported, this is due to
// their notoriously bad performance.
const int MESH_REGULARITY = -1;                   
// This parameter influences the selection of
// candidates in hp-adaptivity. Default value is 1.0. 
const double CONV_EXP = 1.0;                      
// Stopping criterion for adaptivity.
const double ERR_STOP = 3.0;                      
// Adaptivity process stops when the number of degrees of freedom grows over
// this limit. This is mainly to prevent h-adaptivity to go on forever.
const int NDOF_STOP = 60000;                      
// Matrix solver: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
// SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.
MatrixSolverType matrix_solver_type = SOLVER_UMFPACK;  

// Problem parameters.
// Heat source in the material (caused by induction heating).
double HEAT_SRC = 10000.0;                        
double TEMP_INNER = 50;
double HEAT_FLUX_OUTER = -50;
// Steel: E=200 GPa.
const double E = 2e11;                            
const double nu = 0.3;
const double lambda = (E * nu) / ((1 + nu) * (1 - 2*nu));
const double mu = E / (2*(1 + nu));
const double l2m = lambda + 2*mu;
const double rho = 8000;
const double g = 9.81;
// See http://hyperphysics.phy-astr.gsu.edu/hbase/tables/thexp.html.
const double alpha = 13e-6;                       

//  Boundary markers.
const int BDY_BOTTOM = 1;
const int BDY_SIDES = 2;
const int BDY_TOP = 3;
const int BDY_HOLES = 4;

// Weak forms.
#include "definitions.cpp"

int main(int argc, char* argv[])
{
  // Time measurement.
  TimePeriod cpu_time;
  cpu_time.tick();

  // Load the mesh.
  Mesh xmesh, ymesh, tmesh;
  MeshReaderH2D mloader;
  // Master mesh.
  mloader.load("domain.mesh", &xmesh); 

  // Initialize multimesh hp-FEM.
  // Ydisp will share master mesh with xdisp.
  ymesh.copy(&xmesh);                  
  // Temp will share master mesh with xdisp.
  tmesh.copy(&xmesh);                  

  // Initialize boundary conditions.
  BCTypes bc_types_x_y;
  bc_types_x_y.add_bc_dirichlet(BDY_BOTTOM);
  bc_types_x_y.add_bc_neumann(Hermes::vector<int>(BDY_SIDES, BDY_TOP, BDY_HOLES));

  BCTypes bc_types_t;
  bc_types_t.add_bc_dirichlet(BDY_HOLES);
  bc_types_t.add_bc_neumann(Hermes::vector<int>(BDY_SIDES, BDY_TOP, BDY_BOTTOM)); 

  // Enter Dirichlet boundary values.
  BCValues bc_values_x_y;
  bc_values_x_y.add_zero(BDY_BOTTOM);

  BCValues bc_values_t;
  bc_values_t.add_const(BDY_HOLES, TEMP_INNER);

  // Create H1 spaces with default shapesets.
  H1Space<double> xdisp(&xmesh, &bc_types_x_y, &bc_values_x_y, P_INIT_DISP);
  H1Space<double> ydisp(MULTI ? &ymesh : &xmesh, &bc_types_x_y, &bc_values_x_y, P_INIT_DISP);
  H1Space<double> temp(MULTI ? &tmesh : &xmesh, &bc_types_t, &bc_values_t, P_INIT_TEMP);

  // Initialize the weak formulation.
  WeakForm wf(3);
  wf.add_matrix_form(0, 0, callback(bilinear_form_0_0));
  wf.add_matrix_form(0, 1, callback(bilinear_form_0_1), HERMES_SYM);
  wf.add_matrix_form(0, 2, callback(bilinear_form_0_2));
  wf.add_matrix_form(1, 1, callback(bilinear_form_1_1));
  wf.add_matrix_form(1, 2, callback(bilinear_form_1_2));
  wf.add_matrix_form(2, 2, callback(bilinear_form_2_2));
  wf.add_vector_form(1, callback(linear_form_1));
  wf.add_vector_form(2, callback(linear_form_2));
  wf.add_vector_form_surf(2, callback(linear_form_surf_2));

  // Initialize coarse and reference mesh solutions.
  Solution<double> xdisp_sln, ydisp_sln, temp_sln, ref_xdisp_sln, ref_ydisp_sln, ref_temp_sln;

  // Initialize refinement selector.
  H1ProjBasedSelector selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);

  // Initialize views.
  ScalarView s_view_0("Solution[xdisp]", new WinGeom(0, 0, 450, 350));
  s_view_0.show_mesh(false);
  ScalarView s_view_1("Solution[ydisp]", new WinGeom(460, 0, 450, 350));
  s_view_1.show_mesh(false);
  ScalarView s_view_2("Solution[temp]", new WinGeom(920, 0, 450, 350));
  s_view_1.show_mesh(false);
  OrderView  o_view_0("Mesh[xdisp]", new WinGeom(0, 360, 450, 350));
  OrderView  o_view_1("Mesh[ydisp]", new WinGeom(460, 360, 450, 350));
  OrderView  o_view_2("Mesh[temp]", new WinGeom(920, 360, 450, 350));

  // DOF and CPU convergence graphs.
  SimpleGraph graph_dof_est, graph_cpu_est;

  // Adaptivity loop:
  int as = 1; 
  bool done = false;
  do
  {
    info("---- Adaptivity step %d:", as);

    // Construct globally refined reference mesh and setup reference space.
    Hermes::vector<Space<double> *>* ref_spaces = Space<double>::construct_refined_spaces(Hermes::vector<Space<double> *>(&xdisp, &ydisp, &temp));

    // Assemble the reference problem.
    info("Solving on reference mesh.");
    bool is_linear = true;
    DiscreteProblem* dp = new DiscreteProblem(&wf, *ref_spaces, is_linear);
    SparseMatrix<double>* matrix = create_matrix<double>(matrix_solver_type);
    Vector<double>* rhs = create_vector<double>(matrix_solver_type);
    LinearSolver<double>* solver = create_linear_solver<double>(matrix_solver_type, matrix, rhs);
    dp->assemble(matrix, rhs);

    // Time measurement.
    cpu_time.tick();
    
    // Solve the linear system of the reference problem. If successful, obtain the solutions.
    if(solver->solve()) Solution::vector_to_solutions(solver->get_solution(), *ref_spaces, 
                                            Hermes::vector<Solution *>(&ref_xdisp_sln, &ref_ydisp_sln, &ref_temp_sln));
    else error ("Matrix solver failed.\n");
  
    // Time measurement.
    cpu_time.tick();

    // Project the fine mesh solution onto the coarse mesh.
    info("Projecting reference solution on coarse mesh.");
    OGProjection::project_global(Hermes::vector<Space<double> *>(&xdisp, &ydisp, &temp), Hermes::vector<Solution *>(&ref_xdisp_sln, &ref_ydisp_sln, &ref_temp_sln), 
                   Hermes::vector<Solution *>(&xdisp_sln, &ydisp_sln, &temp_sln), matrix_solver_type); 
   
    // View the coarse mesh solution and polynomial orders.
    s_view_0.show(&xdisp_sln); 
    o_view_0.show(&xdisp);
    s_view_1.show(&ydisp_sln); 
    o_view_1.show(&ydisp);
    s_view_2.show(&temp_sln); 
    o_view_2.show(&temp);

    // Skip visualization time.
    cpu_time.tick(HERMES_SKIP);

    // Calculate element errors.
    info("Calculating error estimate and exact error."); 
    Adapt* adaptivity = new Adapt(Hermes::vector<Space<double> *>(&xdisp, &ydisp, &temp));
    adaptivity->set_error_form(0, 0, bilinear_form_0_0<double, double>, bilinear_form_0_0<Ord, Ord>);
    adaptivity->set_error_form(0, 1, bilinear_form_0_1<double, double>, bilinear_form_0_1<Ord, Ord>);
    adaptivity->set_error_form(0, 2, bilinear_form_0_2<double, double>, bilinear_form_0_2<Ord, Ord>);
    adaptivity->set_error_form(1, 0, bilinear_form_1_0<double, double>, bilinear_form_1_0<Ord, Ord>);
    adaptivity->set_error_form(1, 1, bilinear_form_1_1<double, double>, bilinear_form_1_1<Ord, Ord>);
    adaptivity->set_error_form(1, 2, bilinear_form_1_2<double, double>, bilinear_form_1_2<Ord, Ord>);
    adaptivity->set_error_form(2, 2, bilinear_form_2_2<double, double>, bilinear_form_2_2<Ord, Ord>);

    // Calculate error estimate for each solution component and the total error estimate.
    Hermes::vector<double> err_est_rel;
    double err_est_rel_total = adaptivity->calc_err_est(Hermes::vector<Solution *>(&xdisp_sln, &ydisp_sln, &temp_sln), 
        Hermes::vector<Solution *>(&ref_xdisp_sln, &ref_ydisp_sln, &ref_temp_sln), &err_est_rel) * 100;

    // Time measurement.
    cpu_time.tick();

    // Report results.
    info("ndof_coarse[xdisp]: %d, ndof_fine[xdisp]: %d, err_est_rel[xdisp]: %g%%", 
         xdisp.Space<double>::get_num_dofs(), Space<double>::get_num_dofs((*ref_spaces)[0]), err_est_rel[0]*100);
    info("ndof_coarse[ydisp]: %d, ndof_fine[ydisp]: %d, err_est_rel[ydisp]: %g%%",
         ydisp.Space<double>::get_num_dofs(), Space<double>::get_num_dofs((*ref_spaces)[1]), err_est_rel[1]*100);
    info("ndof_coarse[temp]: %d, ndof_fine[temp]: %d, err_est_rel[temp]: %g%%",
         temp.Space<double>::get_num_dofs(), Space<double>::get_num_dofs((*ref_spaces)[2]), err_est_rel[2]*100);
    info("ndof_coarse_total: %d, ndof_fine_total: %d, err_est_rel_total: %g%%",
         Space<double>::get_num_dofs(Hermes::vector<Space<double> *>(&xdisp, &ydisp, &temp)), 
         Space<double>::get_num_dofs(*ref_spaces), err_est_rel_total);

    // Add entry to DOF and CPU convergence graphs.
    graph_dof_est.add_values(Space<double>::get_num_dofs(Hermes::vector<Space<double> *>(&xdisp, &ydisp, &temp)), err_est_rel_total);
    graph_dof_est.save("conv_dof_est.dat");
    graph_cpu_est.add_values(cpu_time.accumulated(), err_est_rel_total);
    graph_cpu_est.save("conv_cpu_est.dat");

    // If err_est too large, adapt the mesh.
    if (err_est_rel_total < ERR_STOP) 
      done = true;
    else 
    {
      info("Adapting coarse mesh.");
      done = adaptivity->adapt(Hermes::vector<RefinementSelectors::Selector *>(&selector, &selector, &selector), 
                               THRESHOLD, STRATEGY, MESH_REGULARITY);
    }
    if (Space<double>::get_num_dofs(Hermes::vector<Space<double> *>(&xdisp, &ydisp, &temp)) >= NDOF_STOP) done = true;

    // Clean up.
    delete solver;
    delete matrix;
    delete rhs;
    delete adaptivity;
    if(done == false)
      for(unsigned int i = 0; i < ref_spaces->size(); i++)
        delete (*ref_spaces)[i]->get_mesh();
    delete ref_spaces;
    delete dp;
    
    // Increase counter.
    as++;
  }
  while (done == false);

  verbose("Total running time: %g s", cpu_time.accumulated());

  // Show the reference solution - the final result.
  s_view_0.set_title("Fine mesh Solution[xdisp]");
  s_view_0.show(&ref_xdisp_sln);
  s_view_1.set_title("Fine mesh Solution[ydisp]");
  s_view_1.show(&ref_ydisp_sln);
  s_view_1.set_title("Fine mesh Solution[temp]");
  s_view_1.show(&ref_temp_sln);

  // Wait for all views to be closed.
  View::wait();
  return 0;
};

