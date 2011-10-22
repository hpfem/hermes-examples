#include "hermes2d.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;
using namespace Hermes::Hermes2D::RefinementSelectors;

/* Weak forms */

class CustomMatrixFormVol : public MatrixFormVol<double>
{
public:
  CustomMatrixFormVol(int i, int j, Mesh* mesh) 
      : MatrixFormVol<double>(i, j), mesh(mesh) {};

  template<typename Real, typename Scalar>
  Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
      Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const;

  virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u,
      Func<double> *v, Geom<double> *e, ExtData<double> *ext) const;

  virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
      Geom<Ord> *e, ExtData<Ord> *ext) const;

  Mesh* mesh;
};

class CustomVectorFormVol : public VectorFormVol<double>
{
public:
  CustomVectorFormVol(int i, Mesh* mesh) : VectorFormVol<double>(i), mesh(mesh) {};

  template<typename Real, typename Scalar>
  Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[],
      Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const;

  virtual double value(int n, double *wt, Func<double> *u_ext[],
      Func<double> *v, Geom<double> *e, ExtData<double> *ext) const;

  virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, 
      Geom<Ord> *e, ExtData<Ord> *ext) const;

  Mesh* mesh;
};

class CustomMatrixFormSurf : public MatrixFormSurf<double>
{
public:
  CustomMatrixFormSurf(int i, int j, std::string marker) 
      : MatrixFormSurf<double>(i, j, marker) {};

  template<typename Real, typename Scalar>
  Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
      Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const;

  virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u,
      Func<double> *v, Geom<double> *e, ExtData<double> *ext) const;

  virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
      Geom<Ord> *e, ExtData<Ord> *ext) const;
};

class CustomVectorFormSurf : public VectorFormSurf<double>
{
public:
  CustomVectorFormSurf(int i, std::string marker) 
      : VectorFormSurf<double>(i, marker) {};

  template<typename Real, typename Scalar>
  Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], 
      Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const;

  virtual double value(int n, double *wt, Func<double> *u_ext[],
      Func<double> *v, Geom<double> *e, ExtData<double> *ext) const;

  virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, 
      Geom<Ord> *e, ExtData<Ord> *ext) const;
};

class CustomWeakFormPoisson : public WeakForm<double>
{
public:
  CustomWeakFormPoisson(std::string omega_1, std::string omega_2, 
                        std::string omega_3, std::string omega_4, 
                        std::string omega_5, std::string bdy_left, 
                        std::string bdy_top, std::string bdy_right, 
                        std::string bdy_bottom, Mesh* mesh);

  Mesh* mesh;

  const std::string omega_1;
  const std::string omega_2;
  const std::string omega_3;
  const std::string omega_4;
  const std::string omega_5;

  const double p_1;
  const double p_2;
  const double p_3;
  const double p_4;
  const double p_5;

  const double q_1;
  const double q_2;
  const double q_3;
  const double q_4;
  const double q_5;

  const double f_1;
  const double f_2;
  const double f_3;
  const double f_4;
  const double f_5;

  // Boundary markers.
  const std::string bdy_left;
  const std::string bdy_top;
  const std::string bdy_right;
  const std::string bdy_bottom;

  // Boundary condition coefficients for the four sides.
  const double c_left;
  const double c_top;
  const double c_right;
  const double c_bottom;

  const double g_n_left;
  const double g_n_top;
  const double g_n_right;
  const double g_n_bottom;
};

