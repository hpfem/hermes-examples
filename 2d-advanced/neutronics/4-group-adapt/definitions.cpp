////// Weak formulation in axisymmetric coordinate system  ////////////////////////////////////
#include "definitions.h"

CustomWeakForm::CustomWeakForm(const MaterialPropertyMaps& matprop,
  std::vector<MeshFunctionSharedPtr<double> >& iterates,
  double init_keff, std::string bdy_vacuum)
  : DefaultWeakFormSourceIteration<double>(matprop, iterates[0]->get_mesh(), iterates, init_keff, HERMES_AXISYM_Y)
{
  for (unsigned int g = 0; g < matprop.get_G(); g++)
  {
    add_matrix_form_surf(new VacuumBoundaryCondition::Jacobian<double>(g, bdy_vacuum, HERMES_AXISYM_Y));
    add_vector_form_surf(new VacuumBoundaryCondition::Residual<double>(g, bdy_vacuum, HERMES_AXISYM_Y));
  }
}

double ErrorForm::value(int n, double *wt, Func<double> *u_ext[],
  Func<double> *u, Func<double> *v, GeomVol<double> *e,
  Func<double>* *ext) const
{
  switch (this->normType)
  {
  case HERMES_L2_NORM:
    return l2_error_form_axisym<double, double>(n, wt, u_ext, u, v, e, ext);
  case HERMES_H1_NORM:
    return h1_error_form_axisym<double, double>(n, wt, u_ext, u, v, e, ext);
  default:
    throw Hermes::Exceptions::Exception("Only the H1 and L2 norms are currently implemented.");
    return 0.0;
  }
}

Ord ErrorForm::ord(int n, double *wt, Func<Ord> *u_ext[],
  Func<Ord> *u, Func<Ord> *v, GeomVol<Ord> *e,
  Func<Ord>* *ext) const
{
  switch (this->normType)
  {
  case HERMES_L2_NORM:
    return l2_error_form_axisym<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
  case HERMES_H1_NORM:
    return h1_error_form_axisym<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
  default:
    throw Hermes::Exceptions::Exception("Only the H1 and L2 norms are currently implemented.");
    return Ord();
  }
}

// Integral over the active core.
double integrate(MeshFunction<double>* sln, std::string area)
{
  Quad2D* quad = &g_quad_2d_std;
  sln->set_quad_2d(quad);

  double integral = 0.0;
  Element* e;
  MeshSharedPtr mesh = sln->get_mesh();
  int marker = mesh->get_element_markers_conversion().get_internal_marker(area).marker;

  for_all_active_elements(e, mesh)
  {
    if (e->marker == marker)
    {
      update_limit_table(e->get_mode());
      sln->set_active_element(e);
      RefMap* ru = sln->get_refmap();
      int o = sln->get_fn_order() + ru->get_inv_ref_order();
      limit_order(o, e->get_mode());
      sln->set_quad_order(o, H2D_FN_VAL);
      const double *uval = sln->get_fn_values();
      double* x = ru->get_phys_x(o);
      double result = 0.0;
      h1_integrate_expression(x[i] * uval[i]);
      integral += result;
    }
  }

  return 2.0 * M_PI * integral;
}

// Calculate number of negative solution values.
int get_num_of_neg(MeshFunctionSharedPtr<double> sln)
{
  Quad2D* quad = &g_quad_2d_std;
  sln->set_quad_2d(quad);
  Element* e;
  const MeshSharedPtr mesh = sln->get_mesh();

  int n = 0;

  for_all_active_elements(e, mesh)
  {
    update_limit_table(e->get_mode());
    sln->set_active_element(e);
    RefMap* ru = sln->get_refmap();
    int o = sln->get_fn_order() + ru->get_inv_ref_order();
    limit_order(o, e->get_mode());
    sln->set_quad_order(o, H2D_FN_VAL);
    const double *uval = sln->get_fn_values();
    int np = quad->get_num_points(o, e->get_mode());

    for (int i = 0; i < np; i++)
      if (uval[i] < -1e-12)
        n++;
  }

  return n;
}

int power_iteration(const MaterialPropertyMaps& matprop,
  const std::vector<SpaceSharedPtr<double> >& spaces, DefaultWeakFormSourceIteration<double>* wf,
  const std::vector<MeshFunctionSharedPtr<double> >& solutions, const std::string& fission_region,
  double tol)
{
  // Sanity checks.
  if (spaces.size() != solutions.size())
    throw Hermes::Exceptions::Exception("Spaces and solutions supplied to power_iteration do not match.");

  // Number of energy groups.
  int G = spaces.size();

  // Initialize the discrete problem.
  DiscreteProblem<double> dp(wf, spaces);
  int ndof = Space<double>::get_num_dofs(spaces);

  // The following variables will store pointers to solutions obtained at each iteration and will be needed for
  // updating the eigenvalue.
  std::vector<MeshFunctionSharedPtr<double> > new_solutions;
  for (int g = 0; g < G; g++)
    new_solutions.push_back(new Solution<double>(solutions[g]->get_mesh()));

  // Initial coefficient vector for the Newton's method.
  double* coeff_vec = new double[ndof];

  // Force the Jacobian assembling in the first iteration.
  bool Jacobian_changed = true;

  bool eigen_done = false; int it = 0;
  do
  {
    memset(coeff_vec, 0.0, ndof*sizeof(double));

    NewtonSolver<double> newton(&dp);

    try
    {
      newton.solve(coeff_vec);
    }
    catch (Hermes::Exceptions::Exception e)
    {
      e.print_msg();
      throw Hermes::Exceptions::Exception("Newton's iteration failed.");
    };

    // Convert coefficients vector into a set of Solution pointers.
    Solution<double>::vector_to_solutions(newton.get_sln_vector(), spaces, new_solutions);

    // Update fission sources.
    using WeakFormsNeutronics::Multigroup::SupportClasses::SourceFilter;
    SourceFilter new_source(new_solutions, &matprop, fission_region);
    SourceFilter old_source(solutions, &matprop, fission_region);

    // Compute the eigenvalue for current iteration.
    double k_new = wf->get_keff() * (integrate(&new_source, fission_region) / integrate(&old_source, fission_region));

    Hermes::Mixins::Loggable::Static::info("      dominant eigenvalue (est): %g, rel. difference: %g", k_new, fabs((wf->get_keff() - k_new) / k_new));

    // Stopping criterion.
    if (fabs((wf->get_keff() - k_new) / k_new) < tol) eigen_done = true;

    // Update the final eigenvalue.
    wf->update_keff(k_new);

    it++;

    // Store the new eigenvector approximation in the result.
    for (int g = 0; g < G; g++)
      solutions[g]->copy(new_solutions[g]);
  } while (!eigen_done);

  return it;
}