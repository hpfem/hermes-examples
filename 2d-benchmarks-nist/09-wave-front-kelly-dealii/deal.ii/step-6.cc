/* $Id: step-6.cc 22321 2010-10-12 21:53:08Z kanschat $ */
/* Author: Wolfgang Bangerth, University of Heidelberg, 2000 */

/*    $Id: step-6.cc 22321 2010-10-12 21:53:08Z kanschat $       */
/*                                                                */
/*    Copyright (C) 2000, 2001, 2002, 2003, 2004, 2006, 2007, 2008, 2010 by the deal.II authors */
/*                                                                */
/*    This file is subject to QPL and may not be  distributed     */
/*    without copyright and license information. Please refer     */
/*    to the file deal.II/doc/license.html for the  text  and     */
/*    further information on this license.                        */

                                 // @sect3{Include files}
#include <base/quadrature_lib.h>
#include <base/function.h>
#include <base/logstream.h>
#include <lac/vector.h>
#include <lac/full_matrix.h>
#include <lac/sparse_matrix.h>
#include <lac/compressed_sparsity_pattern.h>
#include <lac/solver_cg.h>
#include <lac/sparse_direct.h>
#include <lac/precondition.h>
#include <grid/tria.h>
#include <dofs/dof_handler.h>
#include <grid/grid_generator.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/tria_boundary_lib.h>
#include <dofs/dof_accessor.h>
#include <dofs/dof_tools.h>
#include <fe/fe_values.h>
#include <numerics/vectors.h>
#include <numerics/matrices.h>
#include <numerics/data_out.h>
#include <numerics/fe_field_function.h>

#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <fe/fe_q.h>
#include <grid/grid_out.h>
#include <lac/constraint_matrix.h>
#include <grid/grid_refinement.h>
#include <numerics/error_estimator.h>
#include <numerics/vectors.h>
#include <base/table_handler.h>
#include <base/timer.h>

using namespace dealii;

template <int dim>
class SolutionBase 
{
  public:
    SolutionBase(int param);
  protected:
    double alpha, x_loc, y_loc, r_zero;
};

template <int dim>
SolutionBase<dim>::SolutionBase(int param)
{
  switch(param) 
  {
    case 0:
      alpha = 20;
      x_loc = -0.05;
      y_loc = -0.05;
      r_zero = 0.7;
      break;
    case 1:
      alpha = 1000;
      x_loc = -0.05;
      y_loc = -0.05;
      r_zero = 0.7;
      break;
    case 2:
      alpha = 1000;
      x_loc = 1.5;
      y_loc = 0.25;
      r_zero = 0.92;
      break;
    case 3:
      alpha = 50;
      x_loc = 0.5;
      y_loc = 0.5;
      r_zero = 0.25;
      break;
    default:   // The same as 0.
      alpha = 20;
      x_loc = -0.05;
      y_loc = -0.05;
      r_zero = 0.7;
      break;
  }
}

template <int dim>
class Solution : public Function<dim>,
     protected SolutionBase<dim>
{
  public:
    Solution (int param) : Function<dim>(), SolutionBase<dim>(param) {}
    
    virtual double value (const Point<dim>   &p,
        const unsigned int  component = 0) const;
    
    virtual Tensor<1,dim> gradient (const Point<dim>   &p,
            const unsigned int  component = 0) const;
};

template <>
double Solution<2>::value (const Point<2>   &p,
           const unsigned int) const
{
  return atan(this->alpha * (sqrt(pow(p[0] - this->x_loc, 2) + pow(p[1] - this->y_loc, 2)) - this->r_zero));
}

template <>
Tensor<1,2> Solution<2>::gradient (const Point<2>   &p,
               const unsigned int) const
{
  double grad[2];
  
  double x = p[0];
  double y = p[1];
  double a = pow(x - this->x_loc, 2);
  double b = pow(y - this->y_loc, 2);
  double c = sqrt(a + b);
  double d = (this->alpha*x - (this->alpha * this->x_loc));
  double e = (this->alpha*y - (this->alpha * this->y_loc));
  double f = (pow(this->alpha*c - (this->alpha * this->r_zero), 2) + 1.0);
  
  grad[0] = (d/(c * f));
  grad[1] = (e/(c * f));
  
  return Tensor<1,2>(grad);
}

template <int dim>
class RightHandSide : public Function<dim>,
          protected SolutionBase<dim>
{
  public:
    RightHandSide (int param) : Function<dim>(), SolutionBase<dim>(param) {}
    
    virtual double value (const Point<dim>   &p,
        const unsigned int  component = 0) const;
};


template <>
double RightHandSide<2>::value (const Point<2>   &p,
          const unsigned int) const
{
  double x = p[0];
  double y = p[1];
  double a = pow(x - this->x_loc, 2);
  double b = pow(y - this->y_loc, 2);
  double c = sqrt(a + b);
  double d = ((this->alpha*x - (this->alpha * this->x_loc)) * (2*x - (2 * this->x_loc)));
  double e = ((this->alpha*y - (this->alpha * this->y_loc)) * (2*y - (2 * this->y_loc)));
  double f = (pow(this->alpha*c - (this->alpha * this->r_zero), 2) + 1.0);
  double g = (this->alpha * c - (this->alpha * this->r_zero));
  
  return -( ((this->alpha/(c * f)) - (d/(2 * pow(a + b, 1.5) * f))
  - ((this->alpha * d * g)/((a + b) * pow(f, 2))) + (this->alpha/(c * f)) 
  - (e/(2 * pow(a + b, 1.5) * f))
  - ((this->alpha * e * g)/((a + b) * pow(f, 2)))));
}

enum SolverType { CG, UMFPACK };
         
template <int dim>
class LaplaceProblem
{
  public:        
    // Refinement modes:
    //
    // 0  ... global refinement
    // 1  ... number, 0.3/0.03
    // 2  ... number, 0.3/0
    // 3  ... number, 0.3/0.03, squared
    // 4  ... number, 0.3/0, squared
    // 5  ... number, sqrt(0.3)/0.03, squared
    // 6  ... number, sqrt(0.3)/0, squared
    // 7  ... fraction, 0.3/0.03
    // 8  ... fraction, 0.3/0
    // 9  ... fraction, 0.3/0.03, squared
    // 10 ... fraction, 0.3/0, squared
    // 11 ... fraction, sqrt(0.3)/0.03, squared
    // 12 ... fraction, sqrt(0.3)/0, squared
    // 13 ... optimized
    //
    LaplaceProblem (int param, SolverType solver_type, int refinement_mode);
    ~LaplaceProblem ();

    void run (double stop);

  private:
    void setup_system ();
    void assemble_system ();
    void solve ();
    void refine_grid ();
    double process_solution (const unsigned int cycle);
    void output_results (const unsigned int cycle) const;
    
    std::string itos(const unsigned int i) const // convert int to string
    {
      std::stringstream s;
      s << i;
      return s.str();
    }

    Triangulation<dim>   triangulation;

    DoFHandler<dim>      dof_handler;
    FE_Q<dim>            fe;

    ConstraintMatrix     hanging_node_constraints;

    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> system_matrix;

    Vector<double>       solution;
    Vector<double>       system_rhs;
    
    TableHandler         convergence_table;
    
    SolverType           solver_type;
    int                  refinement_mode;
    int                  param;
};

template <int dim>
LaplaceProblem<dim>::LaplaceProblem (int param, SolverType solver_type, int refinement_mode)
 : dof_handler (triangulation), fe (2), 
   solver_type (solver_type), param (param), refinement_mode (refinement_mode)
{}


template <int dim>
LaplaceProblem<dim>::~LaplaceProblem ()
{
  dof_handler.clear ();
}

template <int dim>
void LaplaceProblem<dim>::setup_system ()
{
  dof_handler.distribute_dofs (fe);

  solution.reinit (dof_handler.n_dofs());
  system_rhs.reinit (dof_handler.n_dofs());

  hanging_node_constraints.clear ();
  DoFTools::make_hanging_node_constraints (dof_handler,
             hanging_node_constraints);

  hanging_node_constraints.close ();

  CompressedSparsityPattern c_sparsity(dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern (dof_handler, c_sparsity);

  hanging_node_constraints.condense (c_sparsity);
  sparsity_pattern.copy_from(c_sparsity);
  system_matrix.reinit (sparsity_pattern);
}

template <int dim>
void LaplaceProblem<dim>::assemble_system ()
{
  const QGauss<dim>  quadrature_formula(6); // 6-point quadrature is the first one which integrates exactly polynomials of orders up to 10. This order is used in Hermes.

  FEValues<dim> fe_values (fe, quadrature_formula,
         update_values    |  update_gradients |
         update_quadrature_points  |  update_JxW_values);

  const unsigned int   dofs_per_cell = fe.dofs_per_cell;
  const unsigned int   n_q_points    = quadrature_formula.size();

  FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);
  Vector<double>       cell_rhs (dofs_per_cell);

  std::vector<unsigned int> local_dof_indices (dofs_per_cell);

  const RightHandSide<dim> right_hand_side(param);
  std::vector<double>  rhs_values (n_q_points);
  

  typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
  for (; cell!=endc; ++cell)
  {
    cell_matrix = 0;
    cell_rhs = 0;

    fe_values.reinit (cell);

    right_hand_side.value_list (fe_values.get_quadrature_points(),
                                rhs_values);
                                
    for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
      for (unsigned int i=0; i<dofs_per_cell; ++i)
      {
        for (unsigned int j=0; j<dofs_per_cell; ++j)
          cell_matrix(i,j) += (
            fe_values.shape_grad(i,q_point) *
            fe_values.shape_grad(j,q_point) *
            fe_values.JxW(q_point));

          cell_rhs(i) += (fe_values.shape_value(i,q_point) *
              rhs_values [q_point] *
              fe_values.JxW(q_point));
      }

    cell->get_dof_indices (local_dof_indices);
    for (unsigned int i=0; i<dofs_per_cell; ++i)
    {
      for (unsigned int j=0; j<dofs_per_cell; ++j)
        system_matrix.add (local_dof_indices[i],
              local_dof_indices[j],
              cell_matrix(i,j));

      system_rhs(local_dof_indices[i]) += cell_rhs(i);
    }
  }

  hanging_node_constraints.condense (system_matrix);
  hanging_node_constraints.condense (system_rhs);
  
  std::map<unsigned int,double> boundary_values;
  VectorTools::interpolate_boundary_values (dof_handler,
                                            0,
                                            Solution<dim>(param),
                                            boundary_values);
  MatrixTools::apply_boundary_values (boundary_values,
                                      system_matrix,
                                      solution,
                                      system_rhs);
  
  // Hermes DOFs : dof_handler.n_dof() - boundary_values.size();
}

template <int dim>
void LaplaceProblem<dim>::solve ()
{
  if (solver_type == CG)
  {
    SolverControl          solver_control (1000, 1e-12);
    SolverCG<>             solver (solver_control);

    PreconditionSSOR<> preconditioner;
    preconditioner.initialize(system_matrix, 1.2);

    solver.solve (system_matrix, solution, system_rhs, preconditioner);
  }
  else if (solver_type == UMFPACK)
  {
    SparseDirectUMFPACK solver;
    solution = system_rhs;
    solver.initialize<SparseMatrix<double> >(system_matrix);
    solver.solve(solution);
  }

  hanging_node_constraints.distribute (solution);
}

template <int dim>
void LaplaceProblem<dim>::refine_grid ()
{
  Assert (refinement_mod >= 0 && refinement_mode <= 13,
          StandardExceptions::ExcIndexRange(refinement_mode, 0, 13));
  
  if (refinement_mode == 0)
  {
    triangulation.refine_global (1);
    return;
  }
  
  Vector<float> estimated_error_per_cell (triangulation.n_active_cells());
  
  KellyErrorEstimator<dim>::estimate (dof_handler,
                                      QGauss<dim-1>(3),
                                      typename FunctionMap<dim>::type(),
                                      solution,
                                      estimated_error_per_cell);
  switch (refinement_mode)
  {
    case 3:
    case 4:
    case 5:
    case 6:
    case 9:
    case 10:
    case 11:
    case 12:
    {
      // Use squared element errors.
      for (int i = 0; i < estimated_error_per_cell.size(); i++)
        estimated_error_per_cell(i) *= estimated_error_per_cell(i);
      break;
    }
  }
                                      
  switch (refinement_mode) 
  {
    case 1:
    case 3:
    {              
      GridRefinement::refine_and_coarsen_fixed_number(triangulation,
                                                      estimated_error_per_cell,
                                                      0.3, 0.03);
      break;
    }
    
    case 2:
    case 4:
    {              
      GridRefinement::refine_and_coarsen_fixed_number(triangulation,
                                                      estimated_error_per_cell,
                                                      0.3, 0);
      break;
    }
    
    case 5:
    {  
      GridRefinement::refine_and_coarsen_fixed_number(triangulation,
                                                      estimated_error_per_cell,
                                                      sqrt(0.3), 0.03);
      break;
    }
    
    case 6:
    {  
      GridRefinement::refine_and_coarsen_fixed_number(triangulation,
                                                      estimated_error_per_cell,
                                                      sqrt(0.3), 0);
      break;
    }
    
    case 7:
    case 9:
    {              
      GridRefinement::refine_and_coarsen_fixed_fraction(triangulation,
                                                        estimated_error_per_cell,
                                                        0.3, 0.03);
      break;
    }
    
    case 8:
    case 10:
    {              
      GridRefinement::refine_and_coarsen_fixed_fraction(triangulation,
                                                        estimated_error_per_cell,
                                                        0.3, 0);
      break;
    }
    
    case 11:
    {
      GridRefinement::refine_and_coarsen_fixed_fraction(triangulation,
                                                        estimated_error_per_cell,
                                                        sqrt(0.3), 0.03);
      break;
    }
    
    case 12:
    {  
      GridRefinement::refine_and_coarsen_fixed_fraction(triangulation,
                                                        estimated_error_per_cell,
                                                        sqrt(0.3), 0);
      break;
    }
    
    case 13:
    {
      GridRefinement::refine_and_coarsen_optimize(triangulation, 
                                                  estimated_error_per_cell);
      break;
    }
  }
  
  triangulation.execute_coarsening_and_refinement ();
}

template <int dim>
void LaplaceProblem<dim>::output_results (const unsigned int cycle) const
{
  std::string filename = "grid-";
  filename += itos(cycle);
  filename += ".eps";

  std::ofstream output (filename.c_str());

  GridOut grid_out;
  grid_out.write_eps (triangulation, output);
}

template <int dim>
double LaplaceProblem<dim>::process_solution(const unsigned int cycle)
{
  Vector<double> difference_per_cell (triangulation.n_active_cells());
  VectorTools::integrate_difference ( dof_handler,
                                      solution,
                                      Solution<dim>(param),
                                      difference_per_cell,
                                      QGauss<dim>(13),
                                      VectorTools::H1_norm);
  const double H1_error_exact = difference_per_cell.l2_norm();
  
  const unsigned int n_dofs=dof_handler.n_dofs();
  
  std::cout << "Cycle " << cycle << ':' 
            << std::endl
            << "   Number of degrees of freedom: "
            << n_dofs
            << std::endl
            << "   H1 error w.r.t. exact soln.:  "
            << H1_error_exact
            << std::endl;

  convergence_table.add_value("cycle", cycle);
  convergence_table.add_value("H1 error", H1_error_exact);
  convergence_table.add_value("ndof", n_dofs);
  
  return H1_error_exact;
}


template <int dim>
void LaplaceProblem<dim>::run (double stop)
{ 
  double setup_time = 0, assemble_time = 0, solve_time = 0, adapt_time = 0;
  double t_start = 0, t_end;
  double accum_time;
  
  Timer timer; 
  timer.start();
  
  for (unsigned int cycle=0; true; ++cycle)
  {
    //std::cout << "Cycle " << cycle << ':' << std::endl;
    
    t_start = timer.wall_time();
    if (cycle == 0)
    {
      GridGenerator::hyper_cube (triangulation, 0, 1);
      triangulation.refine_global (2);
    }
    else
      refine_grid ();
    t_end = timer.wall_time();
    
    adapt_time += t_end - t_start;
    
    t_start = timer.wall_time();
    setup_system ();
    t_end = timer.wall_time();
    
    setup_time += t_end - t_start;
    
    t_start = timer.wall_time();
    assemble_system ();
    t_end = timer.wall_time();
    
    assemble_time += t_end - t_start;
    
    t_start = timer.wall_time();
    solve ();
    t_end = timer.wall_time();
    
    solve_time += t_end - t_start;
        
    //output_results (cycle);
    
    t_start = timer.wall_time();
    double abs_err = process_solution(cycle);
    t_end = timer.wall_time();
    
    adapt_time += t_end - t_start;
    
    accum_time = timer.wall_time();
    
    convergence_table.add_value("total time", accum_time);
    convergence_table.add_value("setup time", setup_time);
    convergence_table.add_value("assem time", assemble_time);
    convergence_table.add_value("solve time", solve_time);
    convergence_table.add_value("adapt time", adapt_time);
    
    if (abs_err < stop || dof_handler.n_dofs() > 60000)
      break;
  }
  
  timer.stop(); 

  /*  
  DataOutBase::EpsFlags eps_flags;
  eps_flags.z_scaling = 1./4.;
  eps_flags.azimut_angle = 0;
  eps_flags.turn_angle = 0;

  DataOut<dim> data_out;
  data_out.set_flags (eps_flags);

  data_out.attach_dof_handler (dof_handler);
  data_out.add_data_vector (solution, "solution");
  data_out.build_patches ();

  std::ofstream output ("final-solution.eps");
  data_out.write_eps (output);
*/
  
  GridOut grid_out;
  std::ofstream output(("final_mesh-"+itos(refinement_mode)+".eps").c_str());
  GridOutFlags::Eps<2> flags(GridOutFlags::Eps<2>::width, 800);
  grid_out.set_flags(flags);
  grid_out.write_eps (triangulation, output);
  
  convergence_table.set_precision("H1 error", 4);
  convergence_table.set_scientific("H1 error", true);
  convergence_table.set_precision("total time", 3);
  convergence_table.set_precision("setup time", 3);
  convergence_table.set_precision("assem time", 3);
  convergence_table.set_precision("solve time", 3);
  convergence_table.set_precision("adapt time", 3);
    
  std::cout << std::endl;
  convergence_table.write_text(std::cout);
  std::ofstream fout(("conv_table-"+itos(refinement_mode)+".dat").c_str());
  convergence_table.write_text(fout);
  fout.close();  
}

int main ()
{
  try
  {
    deallog.depth_console (0);

    for (int i = 1; i <= 13; i++)
    {
      LaplaceProblem<2> laplace_problem_2d(3, UMFPACK, i);
      laplace_problem_2d.run (0.1);
    }
  }

  catch (std::exception &exc)
  {
    std::cerr << std::endl << std::endl
              << "----------------------------------------------------"
              << std::endl;
                std::cerr << "Exception on processing: " << std::endl
              << exc.what() << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;

    return 1;
  }
  catch (...)
  {
    std::cerr << std::endl << std::endl
              << "----------------------------------------------------"
              << std::endl;
                std::cerr << "Unknown exception!" << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;
    return 1;
  }

  return 0;
}
