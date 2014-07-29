#include "hermes2d.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;

template<typename Scalar>
class volume_matrix_magnetic_harmonic_planar_linear_harmonic_laplace_1_1_1_1 : public MatrixFormVol<Scalar>
{
public:
  volume_matrix_magnetic_harmonic_planar_linear_harmonic_laplace_1_1_1_1(unsigned int i, unsigned int j, int offsetI, int offsetJ);

  virtual Scalar value(int n, double *wt, Hermes::Hermes2D::Func<Scalar> *u_ext[], Hermes::Hermes2D::Func<double> *u,
    Hermes::Hermes2D::Func<double> *v, Hermes::Hermes2D::GeomVol<double> *e, Hermes::Hermes2D::Func<Scalar> **ext) const;
  virtual Hermes::Ord ord(int n, double *wt, Hermes::Hermes2D::Func<Hermes::Ord> *u_ext[], Hermes::Hermes2D::Func<Hermes::Ord> *u,
    Hermes::Hermes2D::Func<Hermes::Ord> *v, Hermes::Hermes2D::GeomVol<Hermes::Ord> *e, Hermes::Hermes2D::Func<Hermes::Ord> **ext) const;

  volume_matrix_magnetic_harmonic_planar_linear_harmonic_laplace_1_1_1_1<Scalar>* clone() const;

private:
};

template<typename Scalar>
class volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_laplace_1_1_1_1 : public MatrixFormVol<Scalar>
{
public:
  volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_laplace_1_1_1_1(unsigned int i, unsigned int j, int offsetI, int offsetJ);

  virtual Scalar value(int n, double *wt, Hermes::Hermes2D::Func<Scalar> *u_ext[], Hermes::Hermes2D::Func<double> *u,
    Hermes::Hermes2D::Func<double> *v, Hermes::Hermes2D::GeomVol<double> *e, Hermes::Hermes2D::Func<Scalar> **ext) const;
  virtual Hermes::Ord ord(int n, double *wt, Hermes::Hermes2D::Func<Hermes::Ord> *u_ext[], Hermes::Hermes2D::Func<Hermes::Ord> *u,
    Hermes::Hermes2D::Func<Hermes::Ord> *v, Hermes::Hermes2D::GeomVol<Hermes::Ord> *e, Hermes::Hermes2D::Func<Hermes::Ord> **ext) const;

  volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_laplace_1_1_1_1<Scalar>* clone() const;

private:
};

template<typename Scalar>
class volume_matrix_magnetic_harmonic_planar_linear_harmonic_laplace_2_2_2_2 : public MatrixFormVol<Scalar>
{
public:
  volume_matrix_magnetic_harmonic_planar_linear_harmonic_laplace_2_2_2_2(unsigned int i, unsigned int j, int offsetI, int offsetJ);

  virtual Scalar value(int n, double *wt, Hermes::Hermes2D::Func<Scalar> *u_ext[], Hermes::Hermes2D::Func<double> *u,
    Hermes::Hermes2D::Func<double> *v, Hermes::Hermes2D::GeomVol<double> *e, Hermes::Hermes2D::Func<Scalar> **ext) const;
  virtual Hermes::Ord ord(int n, double *wt, Hermes::Hermes2D::Func<Hermes::Ord> *u_ext[], Hermes::Hermes2D::Func<Hermes::Ord> *u,
    Hermes::Hermes2D::Func<Hermes::Ord> *v, Hermes::Hermes2D::GeomVol<Hermes::Ord> *e, Hermes::Hermes2D::Func<Hermes::Ord> **ext) const;

  volume_matrix_magnetic_harmonic_planar_linear_harmonic_laplace_2_2_2_2<Scalar>* clone() const;

private:
};

template<typename Scalar>
class volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_laplace_2_2_2_2 : public MatrixFormVol<Scalar>
{
public:
  volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_laplace_2_2_2_2(unsigned int i, unsigned int j, int offsetI, int offsetJ);

  virtual Scalar value(int n, double *wt, Hermes::Hermes2D::Func<Scalar> *u_ext[], Hermes::Hermes2D::Func<double> *u,
    Hermes::Hermes2D::Func<double> *v, Hermes::Hermes2D::GeomVol<double> *e, Hermes::Hermes2D::Func<Scalar> **ext) const;
  virtual Hermes::Ord ord(int n, double *wt, Hermes::Hermes2D::Func<Hermes::Ord> *u_ext[], Hermes::Hermes2D::Func<Hermes::Ord> *u,
    Hermes::Hermes2D::Func<Hermes::Ord> *v, Hermes::Hermes2D::GeomVol<Hermes::Ord> *e, Hermes::Hermes2D::Func<Hermes::Ord> **ext) const;

  volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_laplace_2_2_2_2<Scalar>* clone() const;

private:
};

template<typename Scalar>
class volume_matrix_magnetic_harmonic_planar_linear_harmonic_velocity_1_1_1_1 : public MatrixFormVol<Scalar>
{
public:
  volume_matrix_magnetic_harmonic_planar_linear_harmonic_velocity_1_1_1_1(unsigned int i, unsigned int j, int offsetI, int offsetJ);

  virtual Scalar value(int n, double *wt, Hermes::Hermes2D::Func<Scalar> *u_ext[], Hermes::Hermes2D::Func<double> *u,
    Hermes::Hermes2D::Func<double> *v, Hermes::Hermes2D::GeomVol<double> *e, Hermes::Hermes2D::Func<Scalar> **ext) const;
  virtual Hermes::Ord ord(int n, double *wt, Hermes::Hermes2D::Func<Hermes::Ord> *u_ext[], Hermes::Hermes2D::Func<Hermes::Ord> *u,
    Hermes::Hermes2D::Func<Hermes::Ord> *v, Hermes::Hermes2D::GeomVol<Hermes::Ord> *e, Hermes::Hermes2D::Func<Hermes::Ord> **ext) const;

  volume_matrix_magnetic_harmonic_planar_linear_harmonic_velocity_1_1_1_1<Scalar>* clone() const;

private:
};

template<typename Scalar>
class volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_velocity_1_1_1_1 : public MatrixFormVol<Scalar>
{
public:
  volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_velocity_1_1_1_1(unsigned int i, unsigned int j, int offsetI, int offsetJ);

  virtual Scalar value(int n, double *wt, Hermes::Hermes2D::Func<Scalar> *u_ext[], Hermes::Hermes2D::Func<double> *u,
    Hermes::Hermes2D::Func<double> *v, Hermes::Hermes2D::GeomVol<double> *e, Hermes::Hermes2D::Func<Scalar> **ext) const;
  virtual Hermes::Ord ord(int n, double *wt, Hermes::Hermes2D::Func<Hermes::Ord> *u_ext[], Hermes::Hermes2D::Func<Hermes::Ord> *u,
    Hermes::Hermes2D::Func<Hermes::Ord> *v, Hermes::Hermes2D::GeomVol<Hermes::Ord> *e, Hermes::Hermes2D::Func<Hermes::Ord> **ext) const;

  volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_velocity_1_1_1_1<Scalar>* clone() const;

private:
};

template<typename Scalar>
class volume_matrix_magnetic_harmonic_planar_linear_harmonic_velocity_2_2_2_2 : public MatrixFormVol<Scalar>
{
public:
  volume_matrix_magnetic_harmonic_planar_linear_harmonic_velocity_2_2_2_2(unsigned int i, unsigned int j, int offsetI, int offsetJ);

  virtual Scalar value(int n, double *wt, Hermes::Hermes2D::Func<Scalar> *u_ext[], Hermes::Hermes2D::Func<double> *u,
    Hermes::Hermes2D::Func<double> *v, Hermes::Hermes2D::GeomVol<double> *e, Hermes::Hermes2D::Func<Scalar> **ext) const;
  virtual Hermes::Ord ord(int n, double *wt, Hermes::Hermes2D::Func<Hermes::Ord> *u_ext[], Hermes::Hermes2D::Func<Hermes::Ord> *u,
    Hermes::Hermes2D::Func<Hermes::Ord> *v, Hermes::Hermes2D::GeomVol<Hermes::Ord> *e, Hermes::Hermes2D::Func<Hermes::Ord> **ext) const;

  volume_matrix_magnetic_harmonic_planar_linear_harmonic_velocity_2_2_2_2<Scalar>* clone() const;

private:
};

template<typename Scalar>
class volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_velocity_2_2_2_2 : public MatrixFormVol<Scalar>
{
public:
  volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_velocity_2_2_2_2(unsigned int i, unsigned int j, int offsetI, int offsetJ);

  virtual Scalar value(int n, double *wt, Hermes::Hermes2D::Func<Scalar> *u_ext[], Hermes::Hermes2D::Func<double> *u,
    Hermes::Hermes2D::Func<double> *v, Hermes::Hermes2D::GeomVol<double> *e, Hermes::Hermes2D::Func<Scalar> **ext) const;
  virtual Hermes::Ord ord(int n, double *wt, Hermes::Hermes2D::Func<Hermes::Ord> *u_ext[], Hermes::Hermes2D::Func<Hermes::Ord> *u,
    Hermes::Hermes2D::Func<Hermes::Ord> *v, Hermes::Hermes2D::GeomVol<Hermes::Ord> *e, Hermes::Hermes2D::Func<Hermes::Ord> **ext) const;

  volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_velocity_2_2_2_2<Scalar>* clone() const;

private:
};

template<typename Scalar>
class volume_matrix_magnetic_harmonic_planar_linear_harmonic_gamma_1_2_1_2 : public MatrixFormVol<Scalar>
{
public:
  volume_matrix_magnetic_harmonic_planar_linear_harmonic_gamma_1_2_1_2(unsigned int i, unsigned int j, int offsetI, int offsetJ);

  virtual Scalar value(int n, double *wt, Hermes::Hermes2D::Func<Scalar> *u_ext[], Hermes::Hermes2D::Func<double> *u,
    Hermes::Hermes2D::Func<double> *v, Hermes::Hermes2D::GeomVol<double> *e, Hermes::Hermes2D::Func<Scalar> **ext) const;
  virtual Hermes::Ord ord(int n, double *wt, Hermes::Hermes2D::Func<Hermes::Ord> *u_ext[], Hermes::Hermes2D::Func<Hermes::Ord> *u,
    Hermes::Hermes2D::Func<Hermes::Ord> *v, Hermes::Hermes2D::GeomVol<Hermes::Ord> *e, Hermes::Hermes2D::Func<Hermes::Ord> **ext) const;

  volume_matrix_magnetic_harmonic_planar_linear_harmonic_gamma_1_2_1_2<Scalar>* clone() const;

private:
};

template<typename Scalar>
class volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_gamma_1_2_1_2 : public MatrixFormVol<Scalar>
{
public:
  volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_gamma_1_2_1_2(unsigned int i, unsigned int j, int offsetI, int offsetJ);

  virtual Scalar value(int n, double *wt, Hermes::Hermes2D::Func<Scalar> *u_ext[], Hermes::Hermes2D::Func<double> *u,
    Hermes::Hermes2D::Func<double> *v, Hermes::Hermes2D::GeomVol<double> *e, Hermes::Hermes2D::Func<Scalar> **ext) const;
  virtual Hermes::Ord ord(int n, double *wt, Hermes::Hermes2D::Func<Hermes::Ord> *u_ext[], Hermes::Hermes2D::Func<Hermes::Ord> *u,
    Hermes::Hermes2D::Func<Hermes::Ord> *v, Hermes::Hermes2D::GeomVol<Hermes::Ord> *e, Hermes::Hermes2D::Func<Hermes::Ord> **ext) const;

  volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_gamma_1_2_1_2<Scalar>* clone() const;

private:
};

template<typename Scalar>
class volume_matrix_magnetic_harmonic_planar_linear_harmonic_gamma_2_1_2_1 : public MatrixFormVol<Scalar>
{
public:
  volume_matrix_magnetic_harmonic_planar_linear_harmonic_gamma_2_1_2_1(unsigned int i, unsigned int j, int offsetI, int offsetJ);

  virtual Scalar value(int n, double *wt, Hermes::Hermes2D::Func<Scalar> *u_ext[], Hermes::Hermes2D::Func<double> *u,
    Hermes::Hermes2D::Func<double> *v, Hermes::Hermes2D::GeomVol<double> *e, Hermes::Hermes2D::Func<Scalar> **ext) const;
  virtual Hermes::Ord ord(int n, double *wt, Hermes::Hermes2D::Func<Hermes::Ord> *u_ext[], Hermes::Hermes2D::Func<Hermes::Ord> *u,
    Hermes::Hermes2D::Func<Hermes::Ord> *v, Hermes::Hermes2D::GeomVol<Hermes::Ord> *e, Hermes::Hermes2D::Func<Hermes::Ord> **ext) const;

  volume_matrix_magnetic_harmonic_planar_linear_harmonic_gamma_2_1_2_1<Scalar>* clone() const;

private:
};

template<typename Scalar>
class volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_gamma_2_1_2_1 : public MatrixFormVol<Scalar>
{
public:
  volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_gamma_2_1_2_1(unsigned int i, unsigned int j, int offsetI, int offsetJ);

  virtual Scalar value(int n, double *wt, Hermes::Hermes2D::Func<Scalar> *u_ext[], Hermes::Hermes2D::Func<double> *u,
    Hermes::Hermes2D::Func<double> *v, Hermes::Hermes2D::GeomVol<double> *e, Hermes::Hermes2D::Func<Scalar> **ext) const;
  virtual Hermes::Ord ord(int n, double *wt, Hermes::Hermes2D::Func<Hermes::Ord> *u_ext[], Hermes::Hermes2D::Func<Hermes::Ord> *u,
    Hermes::Hermes2D::Func<Hermes::Ord> *v, Hermes::Hermes2D::GeomVol<Hermes::Ord> *e, Hermes::Hermes2D::Func<Hermes::Ord> **ext) const;

  volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_gamma_2_1_2_1<Scalar>* clone() const;

private:
};

template<typename Scalar>
class volume_matrix_magnetic_harmonic_planar_linear_harmonic_current_1_4_1_4 : public MatrixFormVol<Scalar>
{
public:
  volume_matrix_magnetic_harmonic_planar_linear_harmonic_current_1_4_1_4(unsigned int i, unsigned int j, int offsetI, int offsetJ);

  virtual Scalar value(int n, double *wt, Hermes::Hermes2D::Func<Scalar> *u_ext[], Hermes::Hermes2D::Func<double> *u,
    Hermes::Hermes2D::Func<double> *v, Hermes::Hermes2D::GeomVol<double> *e, Hermes::Hermes2D::Func<Scalar> **ext) const;
  virtual Hermes::Ord ord(int n, double *wt, Hermes::Hermes2D::Func<Hermes::Ord> *u_ext[], Hermes::Hermes2D::Func<Hermes::Ord> *u,
    Hermes::Hermes2D::Func<Hermes::Ord> *v, Hermes::Hermes2D::GeomVol<Hermes::Ord> *e, Hermes::Hermes2D::Func<Hermes::Ord> **ext) const;

  volume_matrix_magnetic_harmonic_planar_linear_harmonic_current_1_4_1_4<Scalar>* clone() const;

private:
};

template<typename Scalar>
class volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_current_1_4_1_4 : public MatrixFormVol<Scalar>
{
public:
  volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_current_1_4_1_4(unsigned int i, unsigned int j, int offsetI, int offsetJ);

  virtual Scalar value(int n, double *wt, Hermes::Hermes2D::Func<Scalar> *u_ext[], Hermes::Hermes2D::Func<double> *u,
    Hermes::Hermes2D::Func<double> *v, Hermes::Hermes2D::GeomVol<double> *e, Hermes::Hermes2D::Func<Scalar> **ext) const;
  virtual Hermes::Ord ord(int n, double *wt, Hermes::Hermes2D::Func<Hermes::Ord> *u_ext[], Hermes::Hermes2D::Func<Hermes::Ord> *u,
    Hermes::Hermes2D::Func<Hermes::Ord> *v, Hermes::Hermes2D::GeomVol<Hermes::Ord> *e, Hermes::Hermes2D::Func<Hermes::Ord> **ext) const;

  volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_current_1_4_1_4<Scalar>* clone() const;

private:
};

template<typename Scalar>
class volume_matrix_magnetic_harmonic_planar_linear_harmonic_current_2_3_2_3 : public MatrixFormVol<Scalar>
{
public:
  volume_matrix_magnetic_harmonic_planar_linear_harmonic_current_2_3_2_3(unsigned int i, unsigned int j, int offsetI, int offsetJ);

  virtual Scalar value(int n, double *wt, Hermes::Hermes2D::Func<Scalar> *u_ext[], Hermes::Hermes2D::Func<double> *u,
    Hermes::Hermes2D::Func<double> *v, Hermes::Hermes2D::GeomVol<double> *e, Hermes::Hermes2D::Func<Scalar> **ext) const;
  virtual Hermes::Ord ord(int n, double *wt, Hermes::Hermes2D::Func<Hermes::Ord> *u_ext[], Hermes::Hermes2D::Func<Hermes::Ord> *u,
    Hermes::Hermes2D::Func<Hermes::Ord> *v, Hermes::Hermes2D::GeomVol<Hermes::Ord> *e, Hermes::Hermes2D::Func<Hermes::Ord> **ext) const;

  volume_matrix_magnetic_harmonic_planar_linear_harmonic_current_2_3_2_3<Scalar>* clone() const;

private:
};

template<typename Scalar>
class volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_current_2_3_2_3 : public MatrixFormVol<Scalar>
{
public:
  volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_current_2_3_2_3(unsigned int i, unsigned int j, int offsetI, int offsetJ);

  virtual Scalar value(int n, double *wt, Hermes::Hermes2D::Func<Scalar> *u_ext[], Hermes::Hermes2D::Func<double> *u,
    Hermes::Hermes2D::Func<double> *v, Hermes::Hermes2D::GeomVol<double> *e, Hermes::Hermes2D::Func<Scalar> **ext) const;
  virtual Hermes::Ord ord(int n, double *wt, Hermes::Hermes2D::Func<Hermes::Ord> *u_ext[], Hermes::Hermes2D::Func<Hermes::Ord> *u,
    Hermes::Hermes2D::Func<Hermes::Ord> *v, Hermes::Hermes2D::GeomVol<Hermes::Ord> *e, Hermes::Hermes2D::Func<Hermes::Ord> **ext) const;

  volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_current_2_3_2_3<Scalar>* clone() const;

private:
};

template<typename Scalar>
class volume_matrix_magnetic_harmonic_planar_linear_harmonic_current_3_1_3_1 : public MatrixFormVol<Scalar>
{
public:
  volume_matrix_magnetic_harmonic_planar_linear_harmonic_current_3_1_3_1(unsigned int i, unsigned int j, int offsetI, int offsetJ);

  virtual Scalar value(int n, double *wt, Hermes::Hermes2D::Func<Scalar> *u_ext[], Hermes::Hermes2D::Func<double> *u,
    Hermes::Hermes2D::Func<double> *v, Hermes::Hermes2D::GeomVol<double> *e, Hermes::Hermes2D::Func<Scalar> **ext) const;
  virtual Hermes::Ord ord(int n, double *wt, Hermes::Hermes2D::Func<Hermes::Ord> *u_ext[], Hermes::Hermes2D::Func<Hermes::Ord> *u,
    Hermes::Hermes2D::Func<Hermes::Ord> *v, Hermes::Hermes2D::GeomVol<Hermes::Ord> *e, Hermes::Hermes2D::Func<Hermes::Ord> **ext) const;

  volume_matrix_magnetic_harmonic_planar_linear_harmonic_current_3_1_3_1<Scalar>* clone() const;

private:
};

template<typename Scalar>
class volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_current_3_1_3_1 : public MatrixFormVol<Scalar>
{
public:
  volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_current_3_1_3_1(unsigned int i, unsigned int j, int offsetI, int offsetJ);

  virtual Scalar value(int n, double *wt, Hermes::Hermes2D::Func<Scalar> *u_ext[], Hermes::Hermes2D::Func<double> *u,
    Hermes::Hermes2D::Func<double> *v, Hermes::Hermes2D::GeomVol<double> *e, Hermes::Hermes2D::Func<Scalar> **ext) const;
  virtual Hermes::Ord ord(int n, double *wt, Hermes::Hermes2D::Func<Hermes::Ord> *u_ext[], Hermes::Hermes2D::Func<Hermes::Ord> *u,
    Hermes::Hermes2D::Func<Hermes::Ord> *v, Hermes::Hermes2D::GeomVol<Hermes::Ord> *e, Hermes::Hermes2D::Func<Hermes::Ord> **ext) const;

  volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_current_3_1_3_1<Scalar>* clone() const;

private:
};

template<typename Scalar>
class volume_matrix_magnetic_harmonic_planar_linear_harmonic_current_4_2_4_2 : public MatrixFormVol<Scalar>
{
public:
  volume_matrix_magnetic_harmonic_planar_linear_harmonic_current_4_2_4_2(unsigned int i, unsigned int j, int offsetI, int offsetJ);

  virtual Scalar value(int n, double *wt, Hermes::Hermes2D::Func<Scalar> *u_ext[], Hermes::Hermes2D::Func<double> *u,
    Hermes::Hermes2D::Func<double> *v, Hermes::Hermes2D::GeomVol<double> *e, Hermes::Hermes2D::Func<Scalar> **ext) const;
  virtual Hermes::Ord ord(int n, double *wt, Hermes::Hermes2D::Func<Hermes::Ord> *u_ext[], Hermes::Hermes2D::Func<Hermes::Ord> *u,
    Hermes::Hermes2D::Func<Hermes::Ord> *v, Hermes::Hermes2D::GeomVol<Hermes::Ord> *e, Hermes::Hermes2D::Func<Hermes::Ord> **ext) const;

  volume_matrix_magnetic_harmonic_planar_linear_harmonic_current_4_2_4_2<Scalar>* clone() const;

private:
};

template<typename Scalar>
class volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_current_4_2_4_2 : public MatrixFormVol<Scalar>
{
public:
  volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_current_4_2_4_2(unsigned int i, unsigned int j, int offsetI, int offsetJ);

  virtual Scalar value(int n, double *wt, Hermes::Hermes2D::Func<Scalar> *u_ext[], Hermes::Hermes2D::Func<double> *u,
    Hermes::Hermes2D::Func<double> *v, Hermes::Hermes2D::GeomVol<double> *e, Hermes::Hermes2D::Func<Scalar> **ext) const;
  virtual Hermes::Ord ord(int n, double *wt, Hermes::Hermes2D::Func<Hermes::Ord> *u_ext[], Hermes::Hermes2D::Func<Hermes::Ord> *u,
    Hermes::Hermes2D::Func<Hermes::Ord> *v, Hermes::Hermes2D::GeomVol<Hermes::Ord> *e, Hermes::Hermes2D::Func<Hermes::Ord> **ext) const;

  volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_current_4_2_4_2<Scalar>* clone() const;

private:
};

template<typename Scalar>
class volume_matrix_magnetic_harmonic_planar_linear_harmonic_current_3_3_3_3 : public MatrixFormVol<Scalar>
{
public:
  volume_matrix_magnetic_harmonic_planar_linear_harmonic_current_3_3_3_3(unsigned int i, unsigned int j, int offsetI, int offsetJ);

  virtual Scalar value(int n, double *wt, Hermes::Hermes2D::Func<Scalar> *u_ext[], Hermes::Hermes2D::Func<double> *u,
    Hermes::Hermes2D::Func<double> *v, Hermes::Hermes2D::GeomVol<double> *e, Hermes::Hermes2D::Func<Scalar> **ext) const;
  virtual Hermes::Ord ord(int n, double *wt, Hermes::Hermes2D::Func<Hermes::Ord> *u_ext[], Hermes::Hermes2D::Func<Hermes::Ord> *u,
    Hermes::Hermes2D::Func<Hermes::Ord> *v, Hermes::Hermes2D::GeomVol<Hermes::Ord> *e, Hermes::Hermes2D::Func<Hermes::Ord> **ext) const;

  volume_matrix_magnetic_harmonic_planar_linear_harmonic_current_3_3_3_3<Scalar>* clone() const;

private:
};

template<typename Scalar>
class volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_current_3_3_3_3 : public MatrixFormVol<Scalar>
{
public:
  volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_current_3_3_3_3(unsigned int i, unsigned int j, int offsetI, int offsetJ);

  virtual Scalar value(int n, double *wt, Hermes::Hermes2D::Func<Scalar> *u_ext[], Hermes::Hermes2D::Func<double> *u,
    Hermes::Hermes2D::Func<double> *v, Hermes::Hermes2D::GeomVol<double> *e, Hermes::Hermes2D::Func<Scalar> **ext) const;
  virtual Hermes::Ord ord(int n, double *wt, Hermes::Hermes2D::Func<Hermes::Ord> *u_ext[], Hermes::Hermes2D::Func<Hermes::Ord> *u,
    Hermes::Hermes2D::Func<Hermes::Ord> *v, Hermes::Hermes2D::GeomVol<Hermes::Ord> *e, Hermes::Hermes2D::Func<Hermes::Ord> **ext) const;

  volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_current_3_3_3_3<Scalar>* clone() const;

private:
};

template<typename Scalar>
class volume_matrix_magnetic_harmonic_planar_linear_harmonic_current_4_4_4_4 : public MatrixFormVol<Scalar>
{
public:
  volume_matrix_magnetic_harmonic_planar_linear_harmonic_current_4_4_4_4(unsigned int i, unsigned int j, int offsetI, int offsetJ);

  virtual Scalar value(int n, double *wt, Hermes::Hermes2D::Func<Scalar> *u_ext[], Hermes::Hermes2D::Func<double> *u,
    Hermes::Hermes2D::Func<double> *v, Hermes::Hermes2D::GeomVol<double> *e, Hermes::Hermes2D::Func<Scalar> **ext) const;
  virtual Hermes::Ord ord(int n, double *wt, Hermes::Hermes2D::Func<Hermes::Ord> *u_ext[], Hermes::Hermes2D::Func<Hermes::Ord> *u,
    Hermes::Hermes2D::Func<Hermes::Ord> *v, Hermes::Hermes2D::GeomVol<Hermes::Ord> *e, Hermes::Hermes2D::Func<Hermes::Ord> **ext) const;

  volume_matrix_magnetic_harmonic_planar_linear_harmonic_current_4_4_4_4<Scalar>* clone() const;

private:
};


template<typename Scalar>
class volume_vector_magnetic_harmonic_planar_linear_harmonic_rhs_1_1_remanence_1_1 : public VectorFormVol<Scalar>
{
public:
  volume_vector_magnetic_harmonic_planar_linear_harmonic_rhs_1_1_remanence_1_1(unsigned int i, unsigned int j, int offsetI, int offsetJ, int* offsetPreviousTimeExt, int* offsetCouplingExt);

  virtual Scalar value(int n, double *wt, Hermes::Hermes2D::Func<Scalar> *u_ext[], Hermes::Hermes2D::Func<double> *v,
    Hermes::Hermes2D::GeomVol<double> *e, Hermes::Hermes2D::Func<Scalar> **ext) const;
  virtual Hermes::Ord ord(int n, double *wt, Hermes::Hermes2D::Func<Hermes::Ord> *u_ext[], Hermes::Hermes2D::Func<Hermes::Ord> *v,
    Hermes::Hermes2D::GeomVol<Hermes::Ord> *e, Hermes::Hermes2D::Func<Hermes::Ord> **ext) const;

  volume_vector_magnetic_harmonic_planar_linear_harmonic_rhs_1_1_remanence_1_1<Scalar>* clone() const;

private:
  unsigned int j;
};

template<typename Scalar>
class volume_vector_magnetic_harmonic_axisymmetric_linear_harmonic_rhs_1_1_remanence_1_1 : public VectorFormVol<Scalar>
{
public:
  volume_vector_magnetic_harmonic_axisymmetric_linear_harmonic_rhs_1_1_remanence_1_1(unsigned int i, unsigned int j, int offsetI, int offsetJ, int* offsetPreviousTimeExt, int* offsetCouplingExt);

  virtual Scalar value(int n, double *wt, Hermes::Hermes2D::Func<Scalar> *u_ext[], Hermes::Hermes2D::Func<double> *v,
    Hermes::Hermes2D::GeomVol<double> *e, Hermes::Hermes2D::Func<Scalar> **ext) const;
  virtual Hermes::Ord ord(int n, double *wt, Hermes::Hermes2D::Func<Hermes::Ord> *u_ext[], Hermes::Hermes2D::Func<Hermes::Ord> *v,
    Hermes::Hermes2D::GeomVol<Hermes::Ord> *e, Hermes::Hermes2D::Func<Hermes::Ord> **ext) const;

  volume_vector_magnetic_harmonic_axisymmetric_linear_harmonic_rhs_1_1_remanence_1_1<Scalar>* clone() const;

private:
  unsigned int j;
};

template<typename Scalar>
class volume_vector_magnetic_harmonic_planar_linear_harmonic_rhs_1_1_current_1_1 : public VectorFormVol<Scalar>
{
public:
  volume_vector_magnetic_harmonic_planar_linear_harmonic_rhs_1_1_current_1_1(unsigned int i, unsigned int j, int offsetI, int offsetJ, int* offsetPreviousTimeExt, int* offsetCouplingExt);

  virtual Scalar value(int n, double *wt, Hermes::Hermes2D::Func<Scalar> *u_ext[], Hermes::Hermes2D::Func<double> *v,
    Hermes::Hermes2D::GeomVol<double> *e, Hermes::Hermes2D::Func<Scalar> **ext) const;
  virtual Hermes::Ord ord(int n, double *wt, Hermes::Hermes2D::Func<Hermes::Ord> *u_ext[], Hermes::Hermes2D::Func<Hermes::Ord> *v,
    Hermes::Hermes2D::GeomVol<Hermes::Ord> *e, Hermes::Hermes2D::Func<Hermes::Ord> **ext) const;

  volume_vector_magnetic_harmonic_planar_linear_harmonic_rhs_1_1_current_1_1<Scalar>* clone() const;

private:
  unsigned int j;
};

template<typename Scalar>
class volume_vector_magnetic_harmonic_axisymmetric_linear_harmonic_rhs_1_1_current_1_1 : public VectorFormVol<Scalar>
{
public:
  volume_vector_magnetic_harmonic_axisymmetric_linear_harmonic_rhs_1_1_current_1_1(unsigned int i, unsigned int j, int offsetI, int offsetJ, int* offsetPreviousTimeExt, int* offsetCouplingExt);

  virtual Scalar value(int n, double *wt, Hermes::Hermes2D::Func<Scalar> *u_ext[], Hermes::Hermes2D::Func<double> *v,
    Hermes::Hermes2D::GeomVol<double> *e, Hermes::Hermes2D::Func<Scalar> **ext) const;
  virtual Hermes::Ord ord(int n, double *wt, Hermes::Hermes2D::Func<Hermes::Ord> *u_ext[], Hermes::Hermes2D::Func<Hermes::Ord> *v,
    Hermes::Hermes2D::GeomVol<Hermes::Ord> *e, Hermes::Hermes2D::Func<Hermes::Ord> **ext) const;

  volume_vector_magnetic_harmonic_axisymmetric_linear_harmonic_rhs_1_1_current_1_1<Scalar>* clone() const;

private:
  unsigned int j;
};

template<typename Scalar>
class volume_vector_magnetic_harmonic_planar_linear_harmonic_rhs_2_2_current_2_2 : public VectorFormVol<Scalar>
{
public:
  volume_vector_magnetic_harmonic_planar_linear_harmonic_rhs_2_2_current_2_2(unsigned int i, unsigned int j, int offsetI, int offsetJ, int* offsetPreviousTimeExt, int* offsetCouplingExt);

  virtual Scalar value(int n, double *wt, Hermes::Hermes2D::Func<Scalar> *u_ext[], Hermes::Hermes2D::Func<double> *v,
    Hermes::Hermes2D::GeomVol<double> *e, Hermes::Hermes2D::Func<Scalar> **ext) const;
  virtual Hermes::Ord ord(int n, double *wt, Hermes::Hermes2D::Func<Hermes::Ord> *u_ext[], Hermes::Hermes2D::Func<Hermes::Ord> *v,
    Hermes::Hermes2D::GeomVol<Hermes::Ord> *e, Hermes::Hermes2D::Func<Hermes::Ord> **ext) const;

  volume_vector_magnetic_harmonic_planar_linear_harmonic_rhs_2_2_current_2_2<Scalar>* clone() const;

private:
  unsigned int j;
};

template<typename Scalar>
class volume_vector_magnetic_harmonic_axisymmetric_linear_harmonic_rhs_2_2_current_2_2 : public VectorFormVol<Scalar>
{
public:
  volume_vector_magnetic_harmonic_axisymmetric_linear_harmonic_rhs_2_2_current_2_2(unsigned int i, unsigned int j, int offsetI, int offsetJ, int* offsetPreviousTimeExt, int* offsetCouplingExt);

  virtual Scalar value(int n, double *wt, Hermes::Hermes2D::Func<Scalar> *u_ext[], Hermes::Hermes2D::Func<double> *v,
    Hermes::Hermes2D::GeomVol<double> *e, Hermes::Hermes2D::Func<Scalar> **ext) const;
  virtual Hermes::Ord ord(int n, double *wt, Hermes::Hermes2D::Func<Hermes::Ord> *u_ext[], Hermes::Hermes2D::Func<Hermes::Ord> *v,
    Hermes::Hermes2D::GeomVol<Hermes::Ord> *e, Hermes::Hermes2D::Func<Hermes::Ord> **ext) const;

  volume_vector_magnetic_harmonic_axisymmetric_linear_harmonic_rhs_2_2_current_2_2<Scalar>* clone() const;

private:
  unsigned int j;
};

template<typename Scalar>
class volume_vector_magnetic_harmonic_planar_linear_harmonic_vec_current_3_3_3_3 : public VectorFormVol<Scalar>
{
public:
  volume_vector_magnetic_harmonic_planar_linear_harmonic_vec_current_3_3_3_3(unsigned int i, unsigned int j, int offsetI, int offsetJ, int* offsetPreviousTimeExt, int* offsetCouplingExt);

  virtual Scalar value(int n, double *wt, Hermes::Hermes2D::Func<Scalar> *u_ext[], Hermes::Hermes2D::Func<double> *v,
    Hermes::Hermes2D::GeomVol<double> *e, Hermes::Hermes2D::Func<Scalar> **ext) const;
  virtual Hermes::Ord ord(int n, double *wt, Hermes::Hermes2D::Func<Hermes::Ord> *u_ext[], Hermes::Hermes2D::Func<Hermes::Ord> *v,
    Hermes::Hermes2D::GeomVol<Hermes::Ord> *e, Hermes::Hermes2D::Func<Hermes::Ord> **ext) const;

  volume_vector_magnetic_harmonic_planar_linear_harmonic_vec_current_3_3_3_3<Scalar>* clone() const;

private:
  unsigned int j;
};

template<typename Scalar>
class volume_vector_magnetic_harmonic_axisymmetric_linear_harmonic_vec_current_3_3_3_3 : public VectorFormVol<Scalar>
{
public:
  volume_vector_magnetic_harmonic_axisymmetric_linear_harmonic_vec_current_3_3_3_3(unsigned int i, unsigned int j, int offsetI, int offsetJ, int* offsetPreviousTimeExt, int* offsetCouplingExt);

  virtual Scalar value(int n, double *wt, Hermes::Hermes2D::Func<Scalar> *u_ext[], Hermes::Hermes2D::Func<double> *v,
    Hermes::Hermes2D::GeomVol<double> *e, Hermes::Hermes2D::Func<Scalar> **ext) const;
  virtual Hermes::Ord ord(int n, double *wt, Hermes::Hermes2D::Func<Hermes::Ord> *u_ext[], Hermes::Hermes2D::Func<Hermes::Ord> *v,
    Hermes::Hermes2D::GeomVol<Hermes::Ord> *e, Hermes::Hermes2D::Func<Hermes::Ord> **ext) const;

  volume_vector_magnetic_harmonic_axisymmetric_linear_harmonic_vec_current_3_3_3_3<Scalar>* clone() const;

private:
  unsigned int j;
};

template<typename Scalar>
class volume_vector_magnetic_harmonic_planar_linear_harmonic_vec_current_4_4_4_4 : public VectorFormVol<Scalar>
{
public:
  volume_vector_magnetic_harmonic_planar_linear_harmonic_vec_current_4_4_4_4(unsigned int i, unsigned int j, int offsetI, int offsetJ, int* offsetPreviousTimeExt, int* offsetCouplingExt);

  virtual Scalar value(int n, double *wt, Hermes::Hermes2D::Func<Scalar> *u_ext[], Hermes::Hermes2D::Func<double> *v,
    Hermes::Hermes2D::GeomVol<double> *e, Hermes::Hermes2D::Func<Scalar> **ext) const;
  virtual Hermes::Ord ord(int n, double *wt, Hermes::Hermes2D::Func<Hermes::Ord> *u_ext[], Hermes::Hermes2D::Func<Hermes::Ord> *v,
    Hermes::Hermes2D::GeomVol<Hermes::Ord> *e, Hermes::Hermes2D::Func<Hermes::Ord> **ext) const;

  volume_vector_magnetic_harmonic_planar_linear_harmonic_vec_current_4_4_4_4<Scalar>* clone() const;

private:
  unsigned int j;
};


template<typename Scalar>
class surface_vector_magnetic_harmonic_planar_linear_harmonic_current_1_1_1_magnetic_surface_current : public VectorFormSurf<Scalar>
{
public:
  surface_vector_magnetic_harmonic_planar_linear_harmonic_current_1_1_1_magnetic_surface_current(unsigned int i, unsigned int j, int offsetI, int offsetJ);

  virtual Scalar value(int n, double *wt, Hermes::Hermes2D::Func<Scalar> *u_ext[], Hermes::Hermes2D::Func<double> *v,
    Hermes::Hermes2D::GeomVol<double> *e, Hermes::Hermes2D::Func<Scalar> **ext) const;
  virtual Hermes::Ord ord(int n, double *wt, Hermes::Hermes2D::Func<Hermes::Ord> *u_ext[], Hermes::Hermes2D::Func<Hermes::Ord> *v,
    Hermes::Hermes2D::GeomVol<Hermes::Ord> *e, Hermes::Hermes2D::Func<Hermes::Ord> **ext) const;
  surface_vector_magnetic_harmonic_planar_linear_harmonic_current_1_1_1_magnetic_surface_current<Scalar>* clone() const;

private:
  unsigned int j;


  double Jr;
  double Ji;
};


template<typename Scalar>
class surface_vector_magnetic_harmonic_planar_linear_harmonic_current_2_2_2_magnetic_surface_current : public VectorFormSurf<Scalar>
{
public:
  surface_vector_magnetic_harmonic_planar_linear_harmonic_current_2_2_2_magnetic_surface_current(unsigned int i, unsigned int j, int offsetI, int offsetJ);

  virtual Scalar value(int n, double *wt, Hermes::Hermes2D::Func<Scalar> *u_ext[], Hermes::Hermes2D::Func<double> *v,
    Hermes::Hermes2D::GeomVol<double> *e, Hermes::Hermes2D::Func<Scalar> **ext) const;
  virtual Hermes::Ord ord(int n, double *wt, Hermes::Hermes2D::Func<Hermes::Ord> *u_ext[], Hermes::Hermes2D::Func<Hermes::Ord> *v,
    Hermes::Hermes2D::GeomVol<Hermes::Ord> *e, Hermes::Hermes2D::Func<Hermes::Ord> **ext) const;
  surface_vector_magnetic_harmonic_planar_linear_harmonic_current_2_2_2_magnetic_surface_current<Scalar>* clone() const;

private:
  unsigned int j;


  double Jr;
  double Ji;
};


template<typename Scalar>
class exact_magnetic_harmonic_planar_linear_harmonic_essential_1_1_0_magnetic_potential : public ExactSolutionScalarAgros<Scalar>
{
public:
  exact_magnetic_harmonic_planar_linear_harmonic_essential_1_1_0_magnetic_potential(Hermes::Hermes2D::MeshSharedPtr mesh);

  Scalar value(double x, double y) const;
  void derivatives(double x, double y, Scalar& dx, Scalar& dy) const;

  Hermes::Ord ord(double x, double y) const
  {
    return Hermes::Ord(Hermes::Ord::get_max_order());
  }

private:

  double Ar;
  double Ai;
};


template<typename Scalar>
class exact_magnetic_harmonic_planar_linear_harmonic_essential_2_2_0_magnetic_potential : public ExactSolutionScalarAgros<Scalar>
{
public:
  exact_magnetic_harmonic_planar_linear_harmonic_essential_2_2_0_magnetic_potential(Hermes::Hermes2D::MeshSharedPtr mesh);

  Scalar value(double x, double y) const;
  void derivatives(double x, double y, Scalar& dx, Scalar& dy) const;

  Hermes::Ord ord(double x, double y) const
  {
    return Hermes::Ord(Hermes::Ord::get_max_order());
  }

private:

  double Ar;
  double Ai;
};
