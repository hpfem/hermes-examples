#include "hermes2d.h"
#include "../NIST-util.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;
using namespace Hermes::Hermes2D::RefinementSelectors;

/* Weak forms */

class CustomMatrixFormVol : public MatrixFormVol<double>
{
public:
  CustomMatrixFormVol(int i, int j, MeshSharedPtr mesh) 
      : MatrixFormVol<double>(i, j), mesh(mesh) {};

  template<typename Real, typename Scalar>
  Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
      Func<Real> *v, Geom<Real> *e, Func<Scalar>* *ext) const;

  virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u,
      Func<double> *v, Geom<double> *e, Func<double>* *ext) const;

  virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
      Geom<Ord> *e, Func<Ord>* *ext) const;

  MatrixFormVol<double>* clone() const;

  MeshSharedPtr mesh;
};

class CustomVectorFormVol : public VectorFormVol<double>
{
public:
  CustomVectorFormVol(int i, MeshSharedPtr mesh) : VectorFormVol<double>(i), mesh(mesh) {};

  template<typename Real, typename Scalar>
  Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[],
      Func<Real> *v, Geom<Real> *e, Func<Scalar>* *ext) const;

  virtual double value(int n, double *wt, Func<double> *u_ext[],
      Func<double> *v, Geom<double> *e, Func<double>* *ext) const;

  virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, 
      Geom<Ord> *e, Func<Ord>* *ext) const;

  VectorFormVol<double>* clone() const;

  MeshSharedPtr mesh;
};

class CustomMatrixFormSurf : public MatrixFormSurf<double>
{
public:
  CustomMatrixFormSurf(int i, int j, std::string marker) 
    : MatrixFormSurf<double>(i, j) { this->set_area(marker); };

  template<typename Real, typename Scalar>
  Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
      Func<Real> *v, Geom<Real> *e, Func<Scalar>* *ext) const;

  virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u,
      Func<double> *v, Geom<double> *e, Func<double>* *ext) const;

  virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
      Geom<Ord> *e, Func<Ord>* *ext) const;

  MatrixFormSurf<double>* clone() const;

};

class CustomVectorFormSurf : public VectorFormSurf<double>
{
public:
  CustomVectorFormSurf(int i, std::string marker) 
      : VectorFormSurf<double>(i) { this->set_area(marker); };

  template<typename Real, typename Scalar>
  Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], 
      Func<Real> *v, Geom<Real> *e, Func<Scalar>* *ext) const;

  virtual double value(int n, double *wt, Func<double> *u_ext[],
      Func<double> *v, Geom<double> *e, Func<double>* *ext) const;

  virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, 
      Geom<Ord> *e, Func<Ord>* *ext) const;

  VectorFormSurf<double>* clone() const;
};

class CustomWeakFormPoisson : public WeakForm<double>
{
public:
  CustomWeakFormPoisson(std::string omega_1, std::string omega_2, 
                        std::string omega_3, std::string omega_4, 
                        std::string omega_5, std::string bdy_left, 
                        std::string bdy_top, std::string bdy_right, 
                        std::string bdy_bottom, MeshSharedPtr mesh);

  MeshSharedPtr mesh;

  const std::string omega_1;
  const std::string omega_2;
  const std::string omega_3;
  const std::string omega_4;
  const std::string omega_5;

  double p_1;
  double p_2;
  double p_3;
  double p_4;
  double p_5;

  double q_1;
  double q_2;
  double q_3;
  double q_4;
  double q_5;

  double f_1;
  double f_2;
  double f_3;
  double f_4;
  double f_5;

  // Boundary markers.
  const std::string bdy_left;
  const std::string bdy_top;
  const std::string bdy_right;
  const std::string bdy_bottom;

  // Boundary condition coefficients for the four sides.
  double c_left;
  double c_top;
  double c_right;
  double c_bottom;

  double g_n_left;
  double g_n_top;
  double g_n_right;
  double g_n_bottom;

  virtual WeakForm* clone() const;
};

