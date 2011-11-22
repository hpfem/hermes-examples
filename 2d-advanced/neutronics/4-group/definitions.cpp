////// Weak formulation in axisymmetric coordinate system  ////////////////////////////////////

#include "definitions.h"

CustomWeakForm::CustomWeakForm( const Hermes::Hermes2D::WeakFormsNeutronics::Multigroup::MaterialProperties::Diffusion::MaterialPropertyMaps& matprop,
                                Hermes::vector<MeshFunction<double>*>& iterates,
                                double init_keff, std::string bdy_vacuum )
  :  Hermes::Hermes2D::WeakFormsNeutronics::Multigroup::CompleteWeakForms::Diffusion::DefaultWeakFormSourceIteration<double>(matprop, iterates[0]->get_mesh(), iterates, init_keff, HERMES_AXISYM_Y)
{
  for (unsigned int g = 0; g < matprop.get_G(); g++)
  {
    add_matrix_form_surf(new Hermes::Hermes2D::WeakFormsNeutronics::Multigroup::ElementaryForms::Diffusion::VacuumBoundaryCondition::Jacobian<double>(g, bdy_vacuum, HERMES_AXISYM_Y));
    add_vector_form_surf(new Hermes::Hermes2D::WeakFormsNeutronics::Multigroup::ElementaryForms::Diffusion::VacuumBoundaryCondition::Residual<double>(g, bdy_vacuum, HERMES_AXISYM_Y));
  }
}

// Integral over the active core.
double integrate(MeshFunction<double>* sln, std::string area)
{
  Quad2D* quad = &g_quad_2d_std;
  sln->set_quad_2d(quad);
  
  double integral = 0.0;
  Element* e;
  Mesh* mesh = sln->get_mesh();
  int marker = mesh->get_element_markers_conversion().get_internal_marker(area).marker;
  
  for_all_active_elements(e, mesh)
  {
    if (e->marker == marker)
    {
      update_limit_table(e->get_mode());
      sln->set_active_element(e);
      RefMap* ru = sln->get_refmap();
      int o = 20;
      limit_order(o);
      sln->set_quad_order(o, H2D_FN_VAL);
      double *uval = sln->get_fn_values();
      double* x = ru->get_phys_x(o);
      double result = 0.0;
      h1_integrate_expression(x[i] * uval[i]);
      integral += result;
    }
  }
  
  return 2.0 * M_PI * integral;
}
