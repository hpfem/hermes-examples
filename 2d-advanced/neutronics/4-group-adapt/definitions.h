////// Weak formulation in axisymmetric coordinate system  ////////////////////////////////////

#include "hermes2d.h"

/* Namespaces used */

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;
using namespace Hermes::Hermes2D::RefinementSelectors;
using namespace WeakFormsNeutronics::Multigroup::CompleteWeakForms::Diffusion;

class CustomWeakForm : public DefaultWeakFormSourceIteration<double>
{
  public:
    CustomWeakForm(const MaterialPropertyMaps& matprop,
                   Hermes::vector<MeshFunction<double>*>& iterates,
                   double init_keff, std::string bdy_vacuum);
};

// Integral over the active core.
double integrate(MeshFunction<double>* sln, Mesh* mesh, std::string area);
int get_num_of_neg(MeshFunction<double> *sln);

// Jacobian matrix (same as stiffness matrix since projections are linear).
class H1AxisymProjectionJacobian : public MatrixFormVol<double>
{
public:
  H1AxisymProjectionJacobian(int i) : MatrixFormVol<double>(i, i, HERMES_ANY, HERMES_SYM) {};

  double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v,
                Geom<double> *e, ExtData<double> *ext) const
  {
    return h1_axisym_projection_biform<double, double>(n, wt, u_ext, u, v, e, ext);
  }

  Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
          Geom<Ord> *e, ExtData<Ord> *ext) const
  {
    return h1_axisym_projection_biform<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
  }

private:
  
  template<typename Real, typename Scalar>
  static Scalar h1_axisym_projection_biform(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
                                            Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
  {
    _F_
    Scalar result = (Scalar)0;
    for (int i = 0; i < n; i++)
      result += wt[i] * e->x[i] * (u->val[i] * v->val[i] + u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]);
    return result;
  }
};

// Residual.
class H1AxisymProjectionResidual : public VectorFormVol<double>
{
public:
  H1AxisymProjectionResidual(int i, MeshFunction<double>* ext) : VectorFormVol<double>(i)
  {
    this->ext.push_back(ext);
  }

  double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v,
                Geom<double> *e, ExtData<double> *ext) const
  {
    return h1_axisym_projection_liform<double, double>(n, wt, u_ext, v, e, ext);
  }

  Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
          Geom<Ord> *e, ExtData<Ord> *ext) const
  {
    return h1_axisym_projection_liform<Ord, Ord>(n, wt, u_ext, v, e, ext);
  }

private:
  template<typename Real, typename Scalar>
  Scalar h1_axisym_projection_liform(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v,
                                     Geom<Real> *e, ExtData<Scalar> *ext) const
  {
    _F_
    Scalar result = (Scalar)0;
    for (int i = 0; i < n; i++)
      result += wt[i] * e->x[i] * ( (u_ext[this->i]->val[i] - ext->fn[0]->val[i]) * v->val[i] 
                                  + (u_ext[this->i]->dx[i]  - ext->fn[0]->dx[i])  * v->dx[i] 
                                  + (u_ext[this->i]->dy[i]  - ext->fn[0]->dy[i])  * v->dy[i]  );
    return result;
  }
};

// Matrix forms for error calculation.
class ErrorForm : public Adapt<double>::MatrixFormVolError
{
public:
  ErrorForm(ProjNormType type) : Adapt<double>::MatrixFormVolError(type) {};

  /// Error bilinear form.
  virtual double value(int n, double *wt, Func<double> *u_ext[],
                        Func<double> *u, Func<double> *v, Geom<double> *e,
                        ExtData<double> *ext) const;

  /// Error bilinear form to estimate order of a function.
  virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[],
                  Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e,
                  ExtData<Ord> *ext) const;

private:
  template<typename Real, typename Scalar>
  static Scalar l2_error_form_axisym(int n, double *wt, Func<Scalar> *u_ext[], Func<Scalar> *u,
                                     Func<Scalar> *v, Geom<Real> *e, ExtData<Scalar> *ext)
  {
    Scalar result = (Scalar)0;
    for (int i = 0; i < n; i++)
      result += wt[i] * e->x[i] * (u->val[i] * conj(v->val[i]));
    return result;
  }

  template<typename Real, typename Scalar>
  static Scalar h1_error_form_axisym(int n, double *wt, Func<Scalar> *u_ext[], Func<Scalar> *u,
                                     Func<Scalar> *v, Geom<Real> *e, ExtData<Scalar> *ext)
  {
    Scalar result = (Scalar)0;
    for (int i = 0; i < n; i++)
      result += wt[i] * e->x[i] * (u->val[i] * conj(v->val[i]) + u->dx[i] * conj(v->dx[i])
                                 + u->dy[i] * conj(v->dy[i]));
    return result;
  }  
};

/// \brief Power iteration. 
///
/// Starts from an initial guess stored in the argument 'solutions' and updates it by the final result after the iteration
/// has converged, also updating the global eigenvalue 'k_eff'.
///
/// \param[in]     hermes2d   Class encapsulating global Hermes2D functions.
/// \param[in]     spaces     Pointers to spaces on which the solutions are defined (one space for each energy group).
/// \param[in]     wf         Pointer to the weak form of the problem.
/// \param[in,out] solution   A set of Solution* pointers to solution components (neutron fluxes in each group). 
///                           Initial guess for the iteration on input, converged result on output.
/// \param[in] fission_region String specifiying the part of the solution domain where fission occurs.
/// \param[in]     tol        Relative difference between two successive eigenvalue approximations that stops the iteration.
/// \param[in,out] mat        Pointer to a matrix to which the system associated with the power iteration will be assembled.
/// \param[in,out] rhs        Pointer to a vector to which the right hand sides of the power iteration will be successively assembled.
/// \param[in]     solver     Solver for the resulting matrix problem (specified by \c mat and \c rhs).
///
/// \return  number of iterations needed for convergence within the specified tolerance.
///
int power_iteration(const MaterialPropertyMaps& matprop, 
                    const Hermes::vector<const Space<double>*>& spaces, DefaultWeakFormSourceIteration<double>* wf, 
                    const Hermes::vector<MeshFunction<double> *>& solution, const std::string& fission_region, 
                    double tol, Hermes::MatrixSolverType matrix_solver);
