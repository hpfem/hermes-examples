#include "definitions.h"


template <typename Scalar>
volume_matrix_magnetic_harmonic_planar_linear_harmonic_laplace_1_1_1_1<Scalar>::volume_matrix_magnetic_harmonic_planar_linear_harmonic_laplace_1_1_1_1(unsigned int i, unsigned int j, int offsetI, int offsetJ)
: MatrixFormVolAgros<Scalar>(i, j, offsetI, offsetJ)
{
}


template <typename Scalar>
Scalar volume_matrix_magnetic_harmonic_planar_linear_harmonic_laplace_1_1_1_1<Scalar>::value(int n, double *wt, Hermes::Hermes2D::Func<Scalar> *u_ext[], Hermes::Hermes2D::Func<double> *u,
  Hermes::Hermes2D::Func<double> *v, Hermes::Hermes2D::Geom<double> *e, Hermes::Hermes2D::Func<Scalar> **ext) const
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
  {
    result += wt[i] * (ext[14]->val[i] * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]));
  }
  return result;
}

template <typename Scalar>
Hermes::Ord volume_matrix_magnetic_harmonic_planar_linear_harmonic_laplace_1_1_1_1<Scalar>::ord(int n, double *wt, Hermes::Hermes2D::Func<Hermes::Ord> *u_ext[], Hermes::Hermes2D::Func<Hermes::Ord> *u,
  Hermes::Hermes2D::Func<Hermes::Ord> *v, Hermes::Hermes2D::Geom<Hermes::Ord> *e, Hermes::Hermes2D::Func<Hermes::Ord> **ext) const
{
  Hermes::Ord result(0);
  for (int i = 0; i < n; i++)
  {
    result += wt[i] * (ext[14]->val[i] * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]));
  }
  return result;
}

template <typename Scalar>
volume_matrix_magnetic_harmonic_planar_linear_harmonic_laplace_1_1_1_1<Scalar>* volume_matrix_magnetic_harmonic_planar_linear_harmonic_laplace_1_1_1_1<Scalar>::clone() const
{
  //return new volume_matrix_magnetic_harmonic_planar_linear_harmonic_laplace_1_1_1_1(this->i, this->j, this->m_offsetI, this->m_offsetJ);
  return new volume_matrix_magnetic_harmonic_planar_linear_harmonic_laplace_1_1_1_1(*this);
}


template <typename Scalar>
volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_laplace_1_1_1_1<Scalar>::volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_laplace_1_1_1_1(unsigned int i, unsigned int j, int offsetI, int offsetJ)
: MatrixFormVolAgros<Scalar>(i, j, offsetI, offsetJ)
{
}


template <typename Scalar>
Scalar volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_laplace_1_1_1_1<Scalar>::value(int n, double *wt, Hermes::Hermes2D::Func<Scalar> *u_ext[], Hermes::Hermes2D::Func<double> *u,
  Hermes::Hermes2D::Func<double> *v, Hermes::Hermes2D::Geom<double> *e, Hermes::Hermes2D::Func<Scalar> **ext) const
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
  {
    result += wt[i] * (ext[14]->val[i] * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i] + u->val[i] * v->dx[i] / e->x[i]));
  }
  return result;
}

template <typename Scalar>
Hermes::Ord volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_laplace_1_1_1_1<Scalar>::ord(int n, double *wt, Hermes::Hermes2D::Func<Hermes::Ord> *u_ext[], Hermes::Hermes2D::Func<Hermes::Ord> *u,
  Hermes::Hermes2D::Func<Hermes::Ord> *v, Hermes::Hermes2D::Geom<Hermes::Ord> *e, Hermes::Hermes2D::Func<Hermes::Ord> **ext) const
{
  Hermes::Ord result(0);
  for (int i = 0; i < n; i++)
  {
    result += wt[i] * (ext[14]->val[i] * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i] + u->val[i] * v->dx[i] / e->x[i]));
  }
  return result;
}

template <typename Scalar>
volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_laplace_1_1_1_1<Scalar>* volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_laplace_1_1_1_1<Scalar>::clone() const
{
  //return new volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_laplace_1_1_1_1(this->i, this->j, this->m_offsetI, this->m_offsetJ);
  return new volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_laplace_1_1_1_1(*this);
}


template <typename Scalar>
volume_matrix_magnetic_harmonic_planar_linear_harmonic_laplace_2_2_2_2<Scalar>::volume_matrix_magnetic_harmonic_planar_linear_harmonic_laplace_2_2_2_2(unsigned int i, unsigned int j, int offsetI, int offsetJ)
: MatrixFormVolAgros<Scalar>(i, j, offsetI, offsetJ)
{
}


template <typename Scalar>
Scalar volume_matrix_magnetic_harmonic_planar_linear_harmonic_laplace_2_2_2_2<Scalar>::value(int n, double *wt, Hermes::Hermes2D::Func<Scalar> *u_ext[], Hermes::Hermes2D::Func<double> *u,
  Hermes::Hermes2D::Func<double> *v, Hermes::Hermes2D::Geom<double> *e, Hermes::Hermes2D::Func<Scalar> **ext) const
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
  {
    result += wt[i] * (ext[14]->val[i] * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]));
  }
  return result;
}

template <typename Scalar>
Hermes::Ord volume_matrix_magnetic_harmonic_planar_linear_harmonic_laplace_2_2_2_2<Scalar>::ord(int n, double *wt, Hermes::Hermes2D::Func<Hermes::Ord> *u_ext[], Hermes::Hermes2D::Func<Hermes::Ord> *u,
  Hermes::Hermes2D::Func<Hermes::Ord> *v, Hermes::Hermes2D::Geom<Hermes::Ord> *e, Hermes::Hermes2D::Func<Hermes::Ord> **ext) const
{
  Hermes::Ord result(0);
  for (int i = 0; i < n; i++)
  {
    result += wt[i] * (ext[14]->val[i] * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]));
  }
  return result;
}

template <typename Scalar>
volume_matrix_magnetic_harmonic_planar_linear_harmonic_laplace_2_2_2_2<Scalar>* volume_matrix_magnetic_harmonic_planar_linear_harmonic_laplace_2_2_2_2<Scalar>::clone() const
{
  //return new volume_matrix_magnetic_harmonic_planar_linear_harmonic_laplace_2_2_2_2(this->i, this->j, this->m_offsetI, this->m_offsetJ);
  return new volume_matrix_magnetic_harmonic_planar_linear_harmonic_laplace_2_2_2_2(*this);
}


template <typename Scalar>
volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_laplace_2_2_2_2<Scalar>::volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_laplace_2_2_2_2(unsigned int i, unsigned int j, int offsetI, int offsetJ)
: MatrixFormVolAgros<Scalar>(i, j, offsetI, offsetJ)
{
}


template <typename Scalar>
Scalar volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_laplace_2_2_2_2<Scalar>::value(int n, double *wt, Hermes::Hermes2D::Func<Scalar> *u_ext[], Hermes::Hermes2D::Func<double> *u,
  Hermes::Hermes2D::Func<double> *v, Hermes::Hermes2D::Geom<double> *e, Hermes::Hermes2D::Func<Scalar> **ext) const
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
  {
    result += wt[i] * (ext[14]->val[i] * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i] + u->val[i] * v->dx[i] / e->x[i]));
  }
  return result;
}

template <typename Scalar>
Hermes::Ord volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_laplace_2_2_2_2<Scalar>::ord(int n, double *wt, Hermes::Hermes2D::Func<Hermes::Ord> *u_ext[], Hermes::Hermes2D::Func<Hermes::Ord> *u,
  Hermes::Hermes2D::Func<Hermes::Ord> *v, Hermes::Hermes2D::Geom<Hermes::Ord> *e, Hermes::Hermes2D::Func<Hermes::Ord> **ext) const
{
  Hermes::Ord result(0);
  for (int i = 0; i < n; i++)
  {
    result += wt[i] * (ext[14]->val[i] * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i] + u->val[i] * v->dx[i] / e->x[i]));
  }
  return result;
}

template <typename Scalar>
volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_laplace_2_2_2_2<Scalar>* volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_laplace_2_2_2_2<Scalar>::clone() const
{
  //return new volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_laplace_2_2_2_2(this->i, this->j, this->m_offsetI, this->m_offsetJ);
  return new volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_laplace_2_2_2_2(*this);
}


template <typename Scalar>
volume_matrix_magnetic_harmonic_planar_linear_harmonic_velocity_1_1_1_1<Scalar>::volume_matrix_magnetic_harmonic_planar_linear_harmonic_velocity_1_1_1_1(unsigned int i, unsigned int j, int offsetI, int offsetJ)
: MatrixFormVolAgros<Scalar>(i, j, offsetI, offsetJ)
{
}


template <typename Scalar>
Scalar volume_matrix_magnetic_harmonic_planar_linear_harmonic_velocity_1_1_1_1<Scalar>::value(int n, double *wt, Hermes::Hermes2D::Func<Scalar> *u_ext[], Hermes::Hermes2D::Func<double> *u,
  Hermes::Hermes2D::Func<double> *v, Hermes::Hermes2D::Geom<double> *e, Hermes::Hermes2D::Func<Scalar> **ext) const
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
  {
    result += wt[i] * (-ext[2]->val[i] * u->val[i] * ((ext[5]->val[i] - e->y[i] * ext[7]->val[i])*v->dx[i] + (ext[6]->val[i] + e->x[i] * ext[7]->val[i])*v->dy[i]));
  }
  return result;
}

template <typename Scalar>
Hermes::Ord volume_matrix_magnetic_harmonic_planar_linear_harmonic_velocity_1_1_1_1<Scalar>::ord(int n, double *wt, Hermes::Hermes2D::Func<Hermes::Ord> *u_ext[], Hermes::Hermes2D::Func<Hermes::Ord> *u,
  Hermes::Hermes2D::Func<Hermes::Ord> *v, Hermes::Hermes2D::Geom<Hermes::Ord> *e, Hermes::Hermes2D::Func<Hermes::Ord> **ext) const
{
  Hermes::Ord result(0);
  for (int i = 0; i < n; i++)
  {
    result += wt[i] * (-ext[2]->val[i] * u->val[i] * ((ext[5]->val[i] - e->y[i] * ext[7]->val[i])*v->dx[i] + (ext[6]->val[i] + e->x[i] * ext[7]->val[i])*v->dy[i]));
  }
  return result;
}

template <typename Scalar>
volume_matrix_magnetic_harmonic_planar_linear_harmonic_velocity_1_1_1_1<Scalar>* volume_matrix_magnetic_harmonic_planar_linear_harmonic_velocity_1_1_1_1<Scalar>::clone() const
{
  //return new volume_matrix_magnetic_harmonic_planar_linear_harmonic_velocity_1_1_1_1(this->i, this->j, this->m_offsetI, this->m_offsetJ);
  return new volume_matrix_magnetic_harmonic_planar_linear_harmonic_velocity_1_1_1_1(*this);
}


template <typename Scalar>
volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_velocity_1_1_1_1<Scalar>::volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_velocity_1_1_1_1(unsigned int i, unsigned int j, int offsetI, int offsetJ)
: MatrixFormVolAgros<Scalar>(i, j, offsetI, offsetJ)
{
}


template <typename Scalar>
Scalar volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_velocity_1_1_1_1<Scalar>::value(int n, double *wt, Hermes::Hermes2D::Func<Scalar> *u_ext[], Hermes::Hermes2D::Func<double> *u,
  Hermes::Hermes2D::Func<double> *v, Hermes::Hermes2D::Geom<double> *e, Hermes::Hermes2D::Func<Scalar> **ext) const
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
  {
    result += wt[i] * (-ext[2]->val[i] * u->val[i] * ext[6]->val[i] * v->dy[i]);
  }
  return result;
}

template <typename Scalar>
Hermes::Ord volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_velocity_1_1_1_1<Scalar>::ord(int n, double *wt, Hermes::Hermes2D::Func<Hermes::Ord> *u_ext[], Hermes::Hermes2D::Func<Hermes::Ord> *u,
  Hermes::Hermes2D::Func<Hermes::Ord> *v, Hermes::Hermes2D::Geom<Hermes::Ord> *e, Hermes::Hermes2D::Func<Hermes::Ord> **ext) const
{
  Hermes::Ord result(0);
  for (int i = 0; i < n; i++)
  {
    result += wt[i] * (-ext[2]->val[i] * u->val[i] * ext[6]->val[i] * v->dy[i]);
  }
  return result;
}

template <typename Scalar>
volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_velocity_1_1_1_1<Scalar>* volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_velocity_1_1_1_1<Scalar>::clone() const
{
  //return new volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_velocity_1_1_1_1(this->i, this->j, this->m_offsetI, this->m_offsetJ);
  return new volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_velocity_1_1_1_1(*this);
}


template <typename Scalar>
volume_matrix_magnetic_harmonic_planar_linear_harmonic_velocity_2_2_2_2<Scalar>::volume_matrix_magnetic_harmonic_planar_linear_harmonic_velocity_2_2_2_2(unsigned int i, unsigned int j, int offsetI, int offsetJ)
: MatrixFormVolAgros<Scalar>(i, j, offsetI, offsetJ)
{
}


template <typename Scalar>
Scalar volume_matrix_magnetic_harmonic_planar_linear_harmonic_velocity_2_2_2_2<Scalar>::value(int n, double *wt, Hermes::Hermes2D::Func<Scalar> *u_ext[], Hermes::Hermes2D::Func<double> *u,
  Hermes::Hermes2D::Func<double> *v, Hermes::Hermes2D::Geom<double> *e, Hermes::Hermes2D::Func<Scalar> **ext) const
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
  {
    result += wt[i] * (-ext[2]->val[i] * u->val[i] * ((ext[5]->val[i] - e->y[i] * ext[7]->val[i])*v->dx[i] + (ext[6]->val[i] + e->x[i] * ext[7]->val[i])*v->dy[i]));
  }
  return result;
}

template <typename Scalar>
Hermes::Ord volume_matrix_magnetic_harmonic_planar_linear_harmonic_velocity_2_2_2_2<Scalar>::ord(int n, double *wt, Hermes::Hermes2D::Func<Hermes::Ord> *u_ext[], Hermes::Hermes2D::Func<Hermes::Ord> *u,
  Hermes::Hermes2D::Func<Hermes::Ord> *v, Hermes::Hermes2D::Geom<Hermes::Ord> *e, Hermes::Hermes2D::Func<Hermes::Ord> **ext) const
{
  Hermes::Ord result(0);
  for (int i = 0; i < n; i++)
  {
    result += wt[i] * (-ext[2]->val[i] * u->val[i] * ((ext[5]->val[i] - e->y[i] * ext[7]->val[i])*v->dx[i] + (ext[6]->val[i] + e->x[i] * ext[7]->val[i])*v->dy[i]));
  }
  return result;
}

template <typename Scalar>
volume_matrix_magnetic_harmonic_planar_linear_harmonic_velocity_2_2_2_2<Scalar>* volume_matrix_magnetic_harmonic_planar_linear_harmonic_velocity_2_2_2_2<Scalar>::clone() const
{
  //return new volume_matrix_magnetic_harmonic_planar_linear_harmonic_velocity_2_2_2_2(this->i, this->j, this->m_offsetI, this->m_offsetJ);
  return new volume_matrix_magnetic_harmonic_planar_linear_harmonic_velocity_2_2_2_2(*this);
}


template <typename Scalar>
volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_velocity_2_2_2_2<Scalar>::volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_velocity_2_2_2_2(unsigned int i, unsigned int j, int offsetI, int offsetJ)
: MatrixFormVolAgros<Scalar>(i, j, offsetI, offsetJ)
{
}


template <typename Scalar>
Scalar volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_velocity_2_2_2_2<Scalar>::value(int n, double *wt, Hermes::Hermes2D::Func<Scalar> *u_ext[], Hermes::Hermes2D::Func<double> *u,
  Hermes::Hermes2D::Func<double> *v, Hermes::Hermes2D::Geom<double> *e, Hermes::Hermes2D::Func<Scalar> **ext) const
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
  {
    result += wt[i] * (-ext[2]->val[i] * u->val[i] * ext[6]->val[i] * v->dy[i]);
  }
  return result;
}

template <typename Scalar>
Hermes::Ord volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_velocity_2_2_2_2<Scalar>::ord(int n, double *wt, Hermes::Hermes2D::Func<Hermes::Ord> *u_ext[], Hermes::Hermes2D::Func<Hermes::Ord> *u,
  Hermes::Hermes2D::Func<Hermes::Ord> *v, Hermes::Hermes2D::Geom<Hermes::Ord> *e, Hermes::Hermes2D::Func<Hermes::Ord> **ext) const
{
  Hermes::Ord result(0);
  for (int i = 0; i < n; i++)
  {
    result += wt[i] * (-ext[2]->val[i] * u->val[i] * ext[6]->val[i] * v->dy[i]);
  }
  return result;
}

template <typename Scalar>
volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_velocity_2_2_2_2<Scalar>* volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_velocity_2_2_2_2<Scalar>::clone() const
{
  //return new volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_velocity_2_2_2_2(this->i, this->j, this->m_offsetI, this->m_offsetJ);
  return new volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_velocity_2_2_2_2(*this);
}


template <typename Scalar>
volume_matrix_magnetic_harmonic_planar_linear_harmonic_gamma_1_2_1_2<Scalar>::volume_matrix_magnetic_harmonic_planar_linear_harmonic_gamma_1_2_1_2(unsigned int i, unsigned int j, int offsetI, int offsetJ)
: MatrixFormVolAgros<Scalar>(i, j, offsetI, offsetJ)
{
}


template <typename Scalar>
Scalar volume_matrix_magnetic_harmonic_planar_linear_harmonic_gamma_1_2_1_2<Scalar>::value(int n, double *wt, Hermes::Hermes2D::Func<Scalar> *u_ext[], Hermes::Hermes2D::Func<double> *u,
  Hermes::Hermes2D::Func<double> *v, Hermes::Hermes2D::Geom<double> *e, Hermes::Hermes2D::Func<Scalar> **ext) const
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
  {
    result += wt[i] * (-2 * M_PI*this->m_markerSource->fieldInfo()->frequency()*ext[2]->val[i] * u->val[i] * v->val[i]);
  }
  return result;
}

template <typename Scalar>
Hermes::Ord volume_matrix_magnetic_harmonic_planar_linear_harmonic_gamma_1_2_1_2<Scalar>::ord(int n, double *wt, Hermes::Hermes2D::Func<Hermes::Ord> *u_ext[], Hermes::Hermes2D::Func<Hermes::Ord> *u,
  Hermes::Hermes2D::Func<Hermes::Ord> *v, Hermes::Hermes2D::Geom<Hermes::Ord> *e, Hermes::Hermes2D::Func<Hermes::Ord> **ext) const
{
  Hermes::Ord result(0);
  for (int i = 0; i < n; i++)
  {
    result += wt[i] * (-2 * M_PI*this->m_markerSource->fieldInfo()->frequency()*ext[2]->val[i] * u->val[i] * v->val[i]);
  }
  return result;
}

template <typename Scalar>
volume_matrix_magnetic_harmonic_planar_linear_harmonic_gamma_1_2_1_2<Scalar>* volume_matrix_magnetic_harmonic_planar_linear_harmonic_gamma_1_2_1_2<Scalar>::clone() const
{
  //return new volume_matrix_magnetic_harmonic_planar_linear_harmonic_gamma_1_2_1_2(this->i, this->j, this->m_offsetI, this->m_offsetJ);
  return new volume_matrix_magnetic_harmonic_planar_linear_harmonic_gamma_1_2_1_2(*this);
}


template <typename Scalar>
volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_gamma_1_2_1_2<Scalar>::volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_gamma_1_2_1_2(unsigned int i, unsigned int j, int offsetI, int offsetJ)
: MatrixFormVolAgros<Scalar>(i, j, offsetI, offsetJ)
{
}


template <typename Scalar>
Scalar volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_gamma_1_2_1_2<Scalar>::value(int n, double *wt, Hermes::Hermes2D::Func<Scalar> *u_ext[], Hermes::Hermes2D::Func<double> *u,
  Hermes::Hermes2D::Func<double> *v, Hermes::Hermes2D::Geom<double> *e, Hermes::Hermes2D::Func<Scalar> **ext) const
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
  {
    result += wt[i] * (-2 * M_PI*this->m_markerSource->fieldInfo()->frequency()*ext[2]->val[i] * u->val[i] * v->val[i]);
  }
  return result;
}

template <typename Scalar>
Hermes::Ord volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_gamma_1_2_1_2<Scalar>::ord(int n, double *wt, Hermes::Hermes2D::Func<Hermes::Ord> *u_ext[], Hermes::Hermes2D::Func<Hermes::Ord> *u,
  Hermes::Hermes2D::Func<Hermes::Ord> *v, Hermes::Hermes2D::Geom<Hermes::Ord> *e, Hermes::Hermes2D::Func<Hermes::Ord> **ext) const
{
  Hermes::Ord result(0);
  for (int i = 0; i < n; i++)
  {
    result += wt[i] * (-2 * M_PI*this->m_markerSource->fieldInfo()->frequency()*ext[2]->val[i] * u->val[i] * v->val[i]);
  }
  return result;
}

template <typename Scalar>
volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_gamma_1_2_1_2<Scalar>* volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_gamma_1_2_1_2<Scalar>::clone() const
{
  //return new volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_gamma_1_2_1_2(this->i, this->j, this->m_offsetI, this->m_offsetJ);
  return new volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_gamma_1_2_1_2(*this);
}


template <typename Scalar>
volume_matrix_magnetic_harmonic_planar_linear_harmonic_gamma_2_1_2_1<Scalar>::volume_matrix_magnetic_harmonic_planar_linear_harmonic_gamma_2_1_2_1(unsigned int i, unsigned int j, int offsetI, int offsetJ)
: MatrixFormVolAgros<Scalar>(i, j, offsetI, offsetJ)
{
}


template <typename Scalar>
Scalar volume_matrix_magnetic_harmonic_planar_linear_harmonic_gamma_2_1_2_1<Scalar>::value(int n, double *wt, Hermes::Hermes2D::Func<Scalar> *u_ext[], Hermes::Hermes2D::Func<double> *u,
  Hermes::Hermes2D::Func<double> *v, Hermes::Hermes2D::Geom<double> *e, Hermes::Hermes2D::Func<Scalar> **ext) const
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
  {
    result += wt[i] * (2 * M_PI*this->m_markerSource->fieldInfo()->frequency()*ext[2]->val[i] * u->val[i] * v->val[i]);
  }
  return result;
}

template <typename Scalar>
Hermes::Ord volume_matrix_magnetic_harmonic_planar_linear_harmonic_gamma_2_1_2_1<Scalar>::ord(int n, double *wt, Hermes::Hermes2D::Func<Hermes::Ord> *u_ext[], Hermes::Hermes2D::Func<Hermes::Ord> *u,
  Hermes::Hermes2D::Func<Hermes::Ord> *v, Hermes::Hermes2D::Geom<Hermes::Ord> *e, Hermes::Hermes2D::Func<Hermes::Ord> **ext) const
{
  Hermes::Ord result(0);
  for (int i = 0; i < n; i++)
  {
    result += wt[i] * (2 * M_PI*this->m_markerSource->fieldInfo()->frequency()*ext[2]->val[i] * u->val[i] * v->val[i]);
  }
  return result;
}

template <typename Scalar>
volume_matrix_magnetic_harmonic_planar_linear_harmonic_gamma_2_1_2_1<Scalar>* volume_matrix_magnetic_harmonic_planar_linear_harmonic_gamma_2_1_2_1<Scalar>::clone() const
{
  //return new volume_matrix_magnetic_harmonic_planar_linear_harmonic_gamma_2_1_2_1(this->i, this->j, this->m_offsetI, this->m_offsetJ);
  return new volume_matrix_magnetic_harmonic_planar_linear_harmonic_gamma_2_1_2_1(*this);
}


template <typename Scalar>
volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_gamma_2_1_2_1<Scalar>::volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_gamma_2_1_2_1(unsigned int i, unsigned int j, int offsetI, int offsetJ)
: MatrixFormVolAgros<Scalar>(i, j, offsetI, offsetJ)
{
}


template <typename Scalar>
Scalar volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_gamma_2_1_2_1<Scalar>::value(int n, double *wt, Hermes::Hermes2D::Func<Scalar> *u_ext[], Hermes::Hermes2D::Func<double> *u,
  Hermes::Hermes2D::Func<double> *v, Hermes::Hermes2D::Geom<double> *e, Hermes::Hermes2D::Func<Scalar> **ext) const
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
  {
    result += wt[i] * (2 * M_PI*this->m_markerSource->fieldInfo()->frequency()*ext[2]->val[i] * u->val[i] * v->val[i]);
  }
  return result;
}

template <typename Scalar>
Hermes::Ord volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_gamma_2_1_2_1<Scalar>::ord(int n, double *wt, Hermes::Hermes2D::Func<Hermes::Ord> *u_ext[], Hermes::Hermes2D::Func<Hermes::Ord> *u,
  Hermes::Hermes2D::Func<Hermes::Ord> *v, Hermes::Hermes2D::Geom<Hermes::Ord> *e, Hermes::Hermes2D::Func<Hermes::Ord> **ext) const
{
  Hermes::Ord result(0);
  for (int i = 0; i < n; i++)
  {
    result += wt[i] * (2 * M_PI*this->m_markerSource->fieldInfo()->frequency()*ext[2]->val[i] * u->val[i] * v->val[i]);
  }
  return result;
}

template <typename Scalar>
volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_gamma_2_1_2_1<Scalar>* volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_gamma_2_1_2_1<Scalar>::clone() const
{
  //return new volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_gamma_2_1_2_1(this->i, this->j, this->m_offsetI, this->m_offsetJ);
  return new volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_gamma_2_1_2_1(*this);
}


template <typename Scalar>
volume_matrix_magnetic_harmonic_planar_linear_harmonic_current_1_4_1_4<Scalar>::volume_matrix_magnetic_harmonic_planar_linear_harmonic_current_1_4_1_4(unsigned int i, unsigned int j, int offsetI, int offsetJ)
: MatrixFormVolAgros<Scalar>(i, j, offsetI, offsetJ)
{
}


template <typename Scalar>
Scalar volume_matrix_magnetic_harmonic_planar_linear_harmonic_current_1_4_1_4<Scalar>::value(int n, double *wt, Hermes::Hermes2D::Func<Scalar> *u_ext[], Hermes::Hermes2D::Func<double> *u,
  Hermes::Hermes2D::Func<double> *v, Hermes::Hermes2D::Geom<double> *e, Hermes::Hermes2D::Func<Scalar> **ext) const
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
  {
    result += wt[i] * (-2 * M_PI*this->m_markerSource->fieldInfo()->frequency()*ext[2]->val[i] * ext[10]->val[i] * u->val[i] * v->val[i]);
  }
  return result;
}

template <typename Scalar>
Hermes::Ord volume_matrix_magnetic_harmonic_planar_linear_harmonic_current_1_4_1_4<Scalar>::ord(int n, double *wt, Hermes::Hermes2D::Func<Hermes::Ord> *u_ext[], Hermes::Hermes2D::Func<Hermes::Ord> *u,
  Hermes::Hermes2D::Func<Hermes::Ord> *v, Hermes::Hermes2D::Geom<Hermes::Ord> *e, Hermes::Hermes2D::Func<Hermes::Ord> **ext) const
{
  Hermes::Ord result(0);
  for (int i = 0; i < n; i++)
  {
    result += wt[i] * (-2 * M_PI*this->m_markerSource->fieldInfo()->frequency()*ext[2]->val[i] * ext[10]->val[i] * u->val[i] * v->val[i]);
  }
  return result;
}

template <typename Scalar>
volume_matrix_magnetic_harmonic_planar_linear_harmonic_current_1_4_1_4<Scalar>* volume_matrix_magnetic_harmonic_planar_linear_harmonic_current_1_4_1_4<Scalar>::clone() const
{
  //return new volume_matrix_magnetic_harmonic_planar_linear_harmonic_current_1_4_1_4(this->i, this->j, this->m_offsetI, this->m_offsetJ);
  return new volume_matrix_magnetic_harmonic_planar_linear_harmonic_current_1_4_1_4(*this);
}


template <typename Scalar>
volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_current_1_4_1_4<Scalar>::volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_current_1_4_1_4(unsigned int i, unsigned int j, int offsetI, int offsetJ)
: MatrixFormVolAgros<Scalar>(i, j, offsetI, offsetJ)
{
}


template <typename Scalar>
Scalar volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_current_1_4_1_4<Scalar>::value(int n, double *wt, Hermes::Hermes2D::Func<Scalar> *u_ext[], Hermes::Hermes2D::Func<double> *u,
  Hermes::Hermes2D::Func<double> *v, Hermes::Hermes2D::Geom<double> *e, Hermes::Hermes2D::Func<Scalar> **ext) const
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
  {
    result += wt[i] * (-2 * M_PI*this->m_markerSource->fieldInfo()->frequency()*ext[2]->val[i] * ext[10]->val[i] * u->val[i] * v->val[i] / (2 * M_PI*(e->x[i] + 1e-12)));
  }
  return result;
}

template <typename Scalar>
Hermes::Ord volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_current_1_4_1_4<Scalar>::ord(int n, double *wt, Hermes::Hermes2D::Func<Hermes::Ord> *u_ext[], Hermes::Hermes2D::Func<Hermes::Ord> *u,
  Hermes::Hermes2D::Func<Hermes::Ord> *v, Hermes::Hermes2D::Geom<Hermes::Ord> *e, Hermes::Hermes2D::Func<Hermes::Ord> **ext) const
{
  Hermes::Ord result(0);
  for (int i = 0; i < n; i++)
  {
    result += wt[i] * (-2 * M_PI*this->m_markerSource->fieldInfo()->frequency()*ext[2]->val[i] * ext[10]->val[i] * u->val[i] * v->val[i] / (2 * M_PI*(e->x[i] + 1e-12)));
  }
  return result;
}

template <typename Scalar>
volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_current_1_4_1_4<Scalar>* volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_current_1_4_1_4<Scalar>::clone() const
{
  //return new volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_current_1_4_1_4(this->i, this->j, this->m_offsetI, this->m_offsetJ);
  return new volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_current_1_4_1_4(*this);
}


template <typename Scalar>
volume_matrix_magnetic_harmonic_planar_linear_harmonic_current_2_3_2_3<Scalar>::volume_matrix_magnetic_harmonic_planar_linear_harmonic_current_2_3_2_3(unsigned int i, unsigned int j, int offsetI, int offsetJ)
: MatrixFormVolAgros<Scalar>(i, j, offsetI, offsetJ)
{
}


template <typename Scalar>
Scalar volume_matrix_magnetic_harmonic_planar_linear_harmonic_current_2_3_2_3<Scalar>::value(int n, double *wt, Hermes::Hermes2D::Func<Scalar> *u_ext[], Hermes::Hermes2D::Func<double> *u,
  Hermes::Hermes2D::Func<double> *v, Hermes::Hermes2D::Geom<double> *e, Hermes::Hermes2D::Func<Scalar> **ext) const
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
  {
    result += wt[i] * (2 * M_PI*this->m_markerSource->fieldInfo()->frequency()*ext[2]->val[i] * ext[10]->val[i] * u->val[i] * v->val[i]);
  }
  return result;
}

template <typename Scalar>
Hermes::Ord volume_matrix_magnetic_harmonic_planar_linear_harmonic_current_2_3_2_3<Scalar>::ord(int n, double *wt, Hermes::Hermes2D::Func<Hermes::Ord> *u_ext[], Hermes::Hermes2D::Func<Hermes::Ord> *u,
  Hermes::Hermes2D::Func<Hermes::Ord> *v, Hermes::Hermes2D::Geom<Hermes::Ord> *e, Hermes::Hermes2D::Func<Hermes::Ord> **ext) const
{
  Hermes::Ord result(0);
  for (int i = 0; i < n; i++)
  {
    result += wt[i] * (2 * M_PI*this->m_markerSource->fieldInfo()->frequency()*ext[2]->val[i] * ext[10]->val[i] * u->val[i] * v->val[i]);
  }
  return result;
}

template <typename Scalar>
volume_matrix_magnetic_harmonic_planar_linear_harmonic_current_2_3_2_3<Scalar>* volume_matrix_magnetic_harmonic_planar_linear_harmonic_current_2_3_2_3<Scalar>::clone() const
{
  //return new volume_matrix_magnetic_harmonic_planar_linear_harmonic_current_2_3_2_3(this->i, this->j, this->m_offsetI, this->m_offsetJ);
  return new volume_matrix_magnetic_harmonic_planar_linear_harmonic_current_2_3_2_3(*this);
}


template <typename Scalar>
volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_current_2_3_2_3<Scalar>::volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_current_2_3_2_3(unsigned int i, unsigned int j, int offsetI, int offsetJ)
: MatrixFormVolAgros<Scalar>(i, j, offsetI, offsetJ)
{
}


template <typename Scalar>
Scalar volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_current_2_3_2_3<Scalar>::value(int n, double *wt, Hermes::Hermes2D::Func<Scalar> *u_ext[], Hermes::Hermes2D::Func<double> *u,
  Hermes::Hermes2D::Func<double> *v, Hermes::Hermes2D::Geom<double> *e, Hermes::Hermes2D::Func<Scalar> **ext) const
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
  {
    result += wt[i] * (2 * M_PI*this->m_markerSource->fieldInfo()->frequency()*ext[2]->val[i] * ext[10]->val[i] * u->val[i] * v->val[i] / (2 * M_PI*(e->x[i] + 1e-12)));
  }
  return result;
}

template <typename Scalar>
Hermes::Ord volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_current_2_3_2_3<Scalar>::ord(int n, double *wt, Hermes::Hermes2D::Func<Hermes::Ord> *u_ext[], Hermes::Hermes2D::Func<Hermes::Ord> *u,
  Hermes::Hermes2D::Func<Hermes::Ord> *v, Hermes::Hermes2D::Geom<Hermes::Ord> *e, Hermes::Hermes2D::Func<Hermes::Ord> **ext) const
{
  Hermes::Ord result(0);
  for (int i = 0; i < n; i++)
  {
    result += wt[i] * (2 * M_PI*this->m_markerSource->fieldInfo()->frequency()*ext[2]->val[i] * ext[10]->val[i] * u->val[i] * v->val[i] / (2 * M_PI*(e->x[i] + 1e-12)));
  }
  return result;
}

template <typename Scalar>
volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_current_2_3_2_3<Scalar>* volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_current_2_3_2_3<Scalar>::clone() const
{
  //return new volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_current_2_3_2_3(this->i, this->j, this->m_offsetI, this->m_offsetJ);
  return new volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_current_2_3_2_3(*this);
}


template <typename Scalar>
volume_matrix_magnetic_harmonic_planar_linear_harmonic_current_3_1_3_1<Scalar>::volume_matrix_magnetic_harmonic_planar_linear_harmonic_current_3_1_3_1(unsigned int i, unsigned int j, int offsetI, int offsetJ)
: MatrixFormVolAgros<Scalar>(i, j, offsetI, offsetJ)
{
}


template <typename Scalar>
Scalar volume_matrix_magnetic_harmonic_planar_linear_harmonic_current_3_1_3_1<Scalar>::value(int n, double *wt, Hermes::Hermes2D::Func<Scalar> *u_ext[], Hermes::Hermes2D::Func<double> *u,
  Hermes::Hermes2D::Func<double> *v, Hermes::Hermes2D::Geom<double> *e, Hermes::Hermes2D::Func<Scalar> **ext) const
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
  {
    result += wt[i] * (-2 * M_PI*this->m_markerSource->fieldInfo()->frequency()*ext[2]->val[i] * ext[10]->val[i] * u->val[i]);
  }
  return result;
}

template <typename Scalar>
Hermes::Ord volume_matrix_magnetic_harmonic_planar_linear_harmonic_current_3_1_3_1<Scalar>::ord(int n, double *wt, Hermes::Hermes2D::Func<Hermes::Ord> *u_ext[], Hermes::Hermes2D::Func<Hermes::Ord> *u,
  Hermes::Hermes2D::Func<Hermes::Ord> *v, Hermes::Hermes2D::Geom<Hermes::Ord> *e, Hermes::Hermes2D::Func<Hermes::Ord> **ext) const
{
  Hermes::Ord result(0);
  for (int i = 0; i < n; i++)
  {
    result += wt[i] * (-2 * M_PI*this->m_markerSource->fieldInfo()->frequency()*ext[2]->val[i] * ext[10]->val[i] * u->val[i]);
  }
  return result;
}

template <typename Scalar>
volume_matrix_magnetic_harmonic_planar_linear_harmonic_current_3_1_3_1<Scalar>* volume_matrix_magnetic_harmonic_planar_linear_harmonic_current_3_1_3_1<Scalar>::clone() const
{
  //return new volume_matrix_magnetic_harmonic_planar_linear_harmonic_current_3_1_3_1(this->i, this->j, this->m_offsetI, this->m_offsetJ);
  return new volume_matrix_magnetic_harmonic_planar_linear_harmonic_current_3_1_3_1(*this);
}


template <typename Scalar>
volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_current_3_1_3_1<Scalar>::volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_current_3_1_3_1(unsigned int i, unsigned int j, int offsetI, int offsetJ)
: MatrixFormVolAgros<Scalar>(i, j, offsetI, offsetJ)
{
}


template <typename Scalar>
Scalar volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_current_3_1_3_1<Scalar>::value(int n, double *wt, Hermes::Hermes2D::Func<Scalar> *u_ext[], Hermes::Hermes2D::Func<double> *u,
  Hermes::Hermes2D::Func<double> *v, Hermes::Hermes2D::Geom<double> *e, Hermes::Hermes2D::Func<Scalar> **ext) const
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
  {
    result += wt[i] * (-2 * M_PI*this->m_markerSource->fieldInfo()->frequency()*ext[2]->val[i] * ext[10]->val[i] * u->val[i]);
  }
  return result;
}

template <typename Scalar>
Hermes::Ord volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_current_3_1_3_1<Scalar>::ord(int n, double *wt, Hermes::Hermes2D::Func<Hermes::Ord> *u_ext[], Hermes::Hermes2D::Func<Hermes::Ord> *u,
  Hermes::Hermes2D::Func<Hermes::Ord> *v, Hermes::Hermes2D::Geom<Hermes::Ord> *e, Hermes::Hermes2D::Func<Hermes::Ord> **ext) const
{
  Hermes::Ord result(0);
  for (int i = 0; i < n; i++)
  {
    result += wt[i] * (-2 * M_PI*this->m_markerSource->fieldInfo()->frequency()*ext[2]->val[i] * ext[10]->val[i] * u->val[i]);
  }
  return result;
}

template <typename Scalar>
volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_current_3_1_3_1<Scalar>* volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_current_3_1_3_1<Scalar>::clone() const
{
  //return new volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_current_3_1_3_1(this->i, this->j, this->m_offsetI, this->m_offsetJ);
  return new volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_current_3_1_3_1(*this);
}


template <typename Scalar>
volume_matrix_magnetic_harmonic_planar_linear_harmonic_current_4_2_4_2<Scalar>::volume_matrix_magnetic_harmonic_planar_linear_harmonic_current_4_2_4_2(unsigned int i, unsigned int j, int offsetI, int offsetJ)
: MatrixFormVolAgros<Scalar>(i, j, offsetI, offsetJ)
{
}


template <typename Scalar>
Scalar volume_matrix_magnetic_harmonic_planar_linear_harmonic_current_4_2_4_2<Scalar>::value(int n, double *wt, Hermes::Hermes2D::Func<Scalar> *u_ext[], Hermes::Hermes2D::Func<double> *u,
  Hermes::Hermes2D::Func<double> *v, Hermes::Hermes2D::Geom<double> *e, Hermes::Hermes2D::Func<Scalar> **ext) const
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
  {
    result += wt[i] * (2 * M_PI*this->m_markerSource->fieldInfo()->frequency()*ext[2]->val[i] * ext[10]->val[i] * u->val[i]);
  }
  return result;
}

template <typename Scalar>
Hermes::Ord volume_matrix_magnetic_harmonic_planar_linear_harmonic_current_4_2_4_2<Scalar>::ord(int n, double *wt, Hermes::Hermes2D::Func<Hermes::Ord> *u_ext[], Hermes::Hermes2D::Func<Hermes::Ord> *u,
  Hermes::Hermes2D::Func<Hermes::Ord> *v, Hermes::Hermes2D::Geom<Hermes::Ord> *e, Hermes::Hermes2D::Func<Hermes::Ord> **ext) const
{
  Hermes::Ord result(0);
  for (int i = 0; i < n; i++)
  {
    result += wt[i] * (2 * M_PI*this->m_markerSource->fieldInfo()->frequency()*ext[2]->val[i] * ext[10]->val[i] * u->val[i]);
  }
  return result;
}

template <typename Scalar>
volume_matrix_magnetic_harmonic_planar_linear_harmonic_current_4_2_4_2<Scalar>* volume_matrix_magnetic_harmonic_planar_linear_harmonic_current_4_2_4_2<Scalar>::clone() const
{
  //return new volume_matrix_magnetic_harmonic_planar_linear_harmonic_current_4_2_4_2(this->i, this->j, this->m_offsetI, this->m_offsetJ);
  return new volume_matrix_magnetic_harmonic_planar_linear_harmonic_current_4_2_4_2(*this);
}


template <typename Scalar>
volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_current_4_2_4_2<Scalar>::volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_current_4_2_4_2(unsigned int i, unsigned int j, int offsetI, int offsetJ)
: MatrixFormVolAgros<Scalar>(i, j, offsetI, offsetJ)
{
}


template <typename Scalar>
Scalar volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_current_4_2_4_2<Scalar>::value(int n, double *wt, Hermes::Hermes2D::Func<Scalar> *u_ext[], Hermes::Hermes2D::Func<double> *u,
  Hermes::Hermes2D::Func<double> *v, Hermes::Hermes2D::Geom<double> *e, Hermes::Hermes2D::Func<Scalar> **ext) const
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
  {
    result += wt[i] * (2 * M_PI*this->m_markerSource->fieldInfo()->frequency()*ext[2]->val[i] * ext[10]->val[i] * u->val[i]);
  }
  return result;
}

template <typename Scalar>
Hermes::Ord volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_current_4_2_4_2<Scalar>::ord(int n, double *wt, Hermes::Hermes2D::Func<Hermes::Ord> *u_ext[], Hermes::Hermes2D::Func<Hermes::Ord> *u,
  Hermes::Hermes2D::Func<Hermes::Ord> *v, Hermes::Hermes2D::Geom<Hermes::Ord> *e, Hermes::Hermes2D::Func<Hermes::Ord> **ext) const
{
  Hermes::Ord result(0);
  for (int i = 0; i < n; i++)
  {
    result += wt[i] * (2 * M_PI*this->m_markerSource->fieldInfo()->frequency()*ext[2]->val[i] * ext[10]->val[i] * u->val[i]);
  }
  return result;
}

template <typename Scalar>
volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_current_4_2_4_2<Scalar>* volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_current_4_2_4_2<Scalar>::clone() const
{
  //return new volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_current_4_2_4_2(this->i, this->j, this->m_offsetI, this->m_offsetJ);
  return new volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_current_4_2_4_2(*this);
}


template <typename Scalar>
volume_matrix_magnetic_harmonic_planar_linear_harmonic_current_3_3_3_3<Scalar>::volume_matrix_magnetic_harmonic_planar_linear_harmonic_current_3_3_3_3(unsigned int i, unsigned int j, int offsetI, int offsetJ)
: MatrixFormVolAgros<Scalar>(i, j, offsetI, offsetJ)
{
}


template <typename Scalar>
Scalar volume_matrix_magnetic_harmonic_planar_linear_harmonic_current_3_3_3_3<Scalar>::value(int n, double *wt, Hermes::Hermes2D::Func<Scalar> *u_ext[], Hermes::Hermes2D::Func<double> *u,
  Hermes::Hermes2D::Func<double> *v, Hermes::Hermes2D::Geom<double> *e, Hermes::Hermes2D::Func<Scalar> **ext) const
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
  {
    result += wt[i] * (tern(ext[10]->val[i]>0, -2 * M_PI*this->m_markerSource->fieldInfo()->frequency()*ext[2]->val[i] * u->val[i], 1));
  }
  return result;
}

template <typename Scalar>
Hermes::Ord volume_matrix_magnetic_harmonic_planar_linear_harmonic_current_3_3_3_3<Scalar>::ord(int n, double *wt, Hermes::Hermes2D::Func<Hermes::Ord> *u_ext[], Hermes::Hermes2D::Func<Hermes::Ord> *u,
  Hermes::Hermes2D::Func<Hermes::Ord> *v, Hermes::Hermes2D::Geom<Hermes::Ord> *e, Hermes::Hermes2D::Func<Hermes::Ord> **ext) const
{
  Hermes::Ord result(0);
  for (int i = 0; i < n; i++)
  {
    result += wt[i] * (tern(ext[10]->val[i]>0, -2 * M_PI*this->m_markerSource->fieldInfo()->frequency()*ext[2]->val[i] * u->val[i], 1));
  }
  return result;
}

template <typename Scalar>
volume_matrix_magnetic_harmonic_planar_linear_harmonic_current_3_3_3_3<Scalar>* volume_matrix_magnetic_harmonic_planar_linear_harmonic_current_3_3_3_3<Scalar>::clone() const
{
  //return new volume_matrix_magnetic_harmonic_planar_linear_harmonic_current_3_3_3_3(this->i, this->j, this->m_offsetI, this->m_offsetJ);
  return new volume_matrix_magnetic_harmonic_planar_linear_harmonic_current_3_3_3_3(*this);
}


template <typename Scalar>
volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_current_3_3_3_3<Scalar>::volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_current_3_3_3_3(unsigned int i, unsigned int j, int offsetI, int offsetJ)
: MatrixFormVolAgros<Scalar>(i, j, offsetI, offsetJ)
{
}


template <typename Scalar>
Scalar volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_current_3_3_3_3<Scalar>::value(int n, double *wt, Hermes::Hermes2D::Func<Scalar> *u_ext[], Hermes::Hermes2D::Func<double> *u,
  Hermes::Hermes2D::Func<double> *v, Hermes::Hermes2D::Geom<double> *e, Hermes::Hermes2D::Func<Scalar> **ext) const
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
  {
    result += wt[i] * (tern(ext[10]->val[i]>0, -2 * M_PI*this->m_markerSource->fieldInfo()->frequency()*ext[2]->val[i] * u->val[i] / (2 * M_PI*(e->x[i] + 1e-12)), 1));
  }
  return result;
}

template <typename Scalar>
Hermes::Ord volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_current_3_3_3_3<Scalar>::ord(int n, double *wt, Hermes::Hermes2D::Func<Hermes::Ord> *u_ext[], Hermes::Hermes2D::Func<Hermes::Ord> *u,
  Hermes::Hermes2D::Func<Hermes::Ord> *v, Hermes::Hermes2D::Geom<Hermes::Ord> *e, Hermes::Hermes2D::Func<Hermes::Ord> **ext) const
{
  Hermes::Ord result(0);
  for (int i = 0; i < n; i++)
  {
    result += wt[i] * (tern(ext[10]->val[i]>0, -2 * M_PI*this->m_markerSource->fieldInfo()->frequency()*ext[2]->val[i] * u->val[i] / (2 * M_PI*(e->x[i] + 1e-12)), 1));
  }
  return result;
}

template <typename Scalar>
volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_current_3_3_3_3<Scalar>* volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_current_3_3_3_3<Scalar>::clone() const
{
  //return new volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_current_3_3_3_3(this->i, this->j, this->m_offsetI, this->m_offsetJ);
  return new volume_matrix_magnetic_harmonic_axisymmetric_linear_harmonic_current_3_3_3_3(*this);
}


template <typename Scalar>
volume_matrix_magnetic_harmonic_planar_linear_harmonic_current_4_4_4_4<Scalar>::volume_matrix_magnetic_harmonic_planar_linear_harmonic_current_4_4_4_4(unsigned int i, unsigned int j, int offsetI, int offsetJ)
: MatrixFormVolAgros<Scalar>(i, j, offsetI, offsetJ)
{
}


template <typename Scalar>
Scalar volume_matrix_magnetic_harmonic_planar_linear_harmonic_current_4_4_4_4<Scalar>::value(int n, double *wt, Hermes::Hermes2D::Func<Scalar> *u_ext[], Hermes::Hermes2D::Func<double> *u,
  Hermes::Hermes2D::Func<double> *v, Hermes::Hermes2D::Geom<double> *e, Hermes::Hermes2D::Func<Scalar> **ext) const
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
  {
    result += wt[i] * (tern(ext[10]->val[i]>0, 2 * M_PI*this->m_markerSource->fieldInfo()->frequency()*ext[2]->val[i] * u->val[i], 1));
  }
  return result;
}

template <typename Scalar>
Hermes::Ord volume_matrix_magnetic_harmonic_planar_linear_harmonic_current_4_4_4_4<Scalar>::ord(int n, double *wt, Hermes::Hermes2D::Func<Hermes::Ord> *u_ext[], Hermes::Hermes2D::Func<Hermes::Ord> *u,
  Hermes::Hermes2D::Func<Hermes::Ord> *v, Hermes::Hermes2D::Geom<Hermes::Ord> *e, Hermes::Hermes2D::Func<Hermes::Ord> **ext) const
{
  Hermes::Ord result(0);
  for (int i = 0; i < n; i++)
  {
    result += wt[i] * (tern(ext[10]->val[i]>0, 2 * M_PI*this->m_markerSource->fieldInfo()->frequency()*ext[2]->val[i] * u->val[i], 1));
  }
  return result;
}

template <typename Scalar>
volume_matrix_magnetic_harmonic_planar_linear_harmonic_current_4_4_4_4<Scalar>* volume_matrix_magnetic_harmonic_planar_linear_harmonic_current_4_4_4_4<Scalar>::clone() const
{
  //return new volume_matrix_magnetic_harmonic_planar_linear_harmonic_current_4_4_4_4(this->i, this->j, this->m_offsetI, this->m_offsetJ);
  return new volume_matrix_magnetic_harmonic_planar_linear_harmonic_current_4_4_4_4(*this);
}


template <typename Scalar>
volume_vector_magnetic_harmonic_planar_linear_harmonic_rhs_1_1_remanence_1_1<Scalar>::volume_vector_magnetic_harmonic_planar_linear_harmonic_rhs_1_1_remanence_1_1(unsigned int i, unsigned int j, int offsetI, int offsetJ, int* offsetPreviousTimeExt, int* offsetCouplingExt)
: VectorFormVolAgros<Scalar>(i, offsetI, offsetJ, offsetPreviousTimeExt, offsetCouplingExt), j(j)
{
}

template <typename Scalar>
Scalar volume_vector_magnetic_harmonic_planar_linear_harmonic_rhs_1_1_remanence_1_1<Scalar>::value(int n, double *wt, Hermes::Hermes2D::Func<Scalar> *u_ext[], Hermes::Hermes2D::Func<double> *v,
  Hermes::Hermes2D::Geom<double> *e, Hermes::Hermes2D::Func<Scalar> **ext) const
{
  Scalar result = 0;

  for (int i = 0; i < n; i++)
  {
    result += wt[i] * (ext[3]->val[i] * ext[14]->val[i] * (-ext[16]->val[i] * v->dx[i] + ext[17]->val[i] * v->dy[i]));
  }
  return result;
}

template <typename Scalar>
Hermes::Ord volume_vector_magnetic_harmonic_planar_linear_harmonic_rhs_1_1_remanence_1_1<Scalar>::ord(int n, double *wt, Hermes::Hermes2D::Func<Hermes::Ord> *u_ext[], Hermes::Hermes2D::Func<Hermes::Ord> *v,
  Hermes::Hermes2D::Geom<Hermes::Ord> *e, Hermes::Hermes2D::Func<Hermes::Ord> **ext) const
{
  Hermes::Ord result(0);
  for (int i = 0; i < n; i++)
  {
    result += wt[i] * (ext[3]->val[i] * ext[14]->val[i] * (-ext[16]->val[i] * v->dx[i] + ext[17]->val[i] * v->dy[i]));
  }
  return result;
}

template <typename Scalar>
volume_vector_magnetic_harmonic_planar_linear_harmonic_rhs_1_1_remanence_1_1<Scalar>* volume_vector_magnetic_harmonic_planar_linear_harmonic_rhs_1_1_remanence_1_1<Scalar>::clone() const
{
  //return new volume_vector_magnetic_harmonic_planar_linear_harmonic_rhs_1_1_remanence_1_1(this->i, this->j, this->m_offsetI, this->m_offsetJ, this->m_offsetPreviousTimeExt, this->m_offsetCouplingExt);
  return new volume_vector_magnetic_harmonic_planar_linear_harmonic_rhs_1_1_remanence_1_1(*this);
}


template <typename Scalar>
volume_vector_magnetic_harmonic_axisymmetric_linear_harmonic_rhs_1_1_remanence_1_1<Scalar>::volume_vector_magnetic_harmonic_axisymmetric_linear_harmonic_rhs_1_1_remanence_1_1(unsigned int i, unsigned int j, int offsetI, int offsetJ, int* offsetPreviousTimeExt, int* offsetCouplingExt)
: VectorFormVolAgros<Scalar>(i, offsetI, offsetJ, offsetPreviousTimeExt, offsetCouplingExt), j(j)
{
}

template <typename Scalar>
Scalar volume_vector_magnetic_harmonic_axisymmetric_linear_harmonic_rhs_1_1_remanence_1_1<Scalar>::value(int n, double *wt, Hermes::Hermes2D::Func<Scalar> *u_ext[], Hermes::Hermes2D::Func<double> *v,
  Hermes::Hermes2D::Geom<double> *e, Hermes::Hermes2D::Func<Scalar> **ext) const
{
  Scalar result = 0;

  for (int i = 0; i < n; i++)
  {
    result += wt[i] * (-ext[3]->val[i] * ext[14]->val[i] * (-ext[16]->val[i] * v->dx[i] + ext[17]->val[i] * v->dy[i]));
  }
  return result;
}

template <typename Scalar>
Hermes::Ord volume_vector_magnetic_harmonic_axisymmetric_linear_harmonic_rhs_1_1_remanence_1_1<Scalar>::ord(int n, double *wt, Hermes::Hermes2D::Func<Hermes::Ord> *u_ext[], Hermes::Hermes2D::Func<Hermes::Ord> *v,
  Hermes::Hermes2D::Geom<Hermes::Ord> *e, Hermes::Hermes2D::Func<Hermes::Ord> **ext) const
{
  Hermes::Ord result(0);
  for (int i = 0; i < n; i++)
  {
    result += wt[i] * (-ext[3]->val[i] * ext[14]->val[i] * (-ext[16]->val[i] * v->dx[i] + ext[17]->val[i] * v->dy[i]));
  }
  return result;
}

template <typename Scalar>
volume_vector_magnetic_harmonic_axisymmetric_linear_harmonic_rhs_1_1_remanence_1_1<Scalar>* volume_vector_magnetic_harmonic_axisymmetric_linear_harmonic_rhs_1_1_remanence_1_1<Scalar>::clone() const
{
  //return new volume_vector_magnetic_harmonic_axisymmetric_linear_harmonic_rhs_1_1_remanence_1_1(this->i, this->j, this->m_offsetI, this->m_offsetJ, this->m_offsetPreviousTimeExt, this->m_offsetCouplingExt);
  return new volume_vector_magnetic_harmonic_axisymmetric_linear_harmonic_rhs_1_1_remanence_1_1(*this);
}


template <typename Scalar>
volume_vector_magnetic_harmonic_planar_linear_harmonic_rhs_1_1_current_1_1<Scalar>::volume_vector_magnetic_harmonic_planar_linear_harmonic_rhs_1_1_current_1_1(unsigned int i, unsigned int j, int offsetI, int offsetJ, int* offsetPreviousTimeExt, int* offsetCouplingExt)
: VectorFormVolAgros<Scalar>(i, offsetI, offsetJ, offsetPreviousTimeExt, offsetCouplingExt), j(j)
{
}

template <typename Scalar>
Scalar volume_vector_magnetic_harmonic_planar_linear_harmonic_rhs_1_1_current_1_1<Scalar>::value(int n, double *wt, Hermes::Hermes2D::Func<Scalar> *u_ext[], Hermes::Hermes2D::Func<double> *v,
  Hermes::Hermes2D::Geom<double> *e, Hermes::Hermes2D::Func<Scalar> **ext) const
{
  Scalar result = 0;

  for (int i = 0; i < n; i++)
  {
    result += wt[i] * (ext[8]->val[i] * v->val[i]);
  }
  return result;
}

template <typename Scalar>
Hermes::Ord volume_vector_magnetic_harmonic_planar_linear_harmonic_rhs_1_1_current_1_1<Scalar>::ord(int n, double *wt, Hermes::Hermes2D::Func<Hermes::Ord> *u_ext[], Hermes::Hermes2D::Func<Hermes::Ord> *v,
  Hermes::Hermes2D::Geom<Hermes::Ord> *e, Hermes::Hermes2D::Func<Hermes::Ord> **ext) const
{
  Hermes::Ord result(0);
  for (int i = 0; i < n; i++)
  {
    result += wt[i] * (ext[8]->val[i] * v->val[i]);
  }
  return result;
}

template <typename Scalar>
volume_vector_magnetic_harmonic_planar_linear_harmonic_rhs_1_1_current_1_1<Scalar>* volume_vector_magnetic_harmonic_planar_linear_harmonic_rhs_1_1_current_1_1<Scalar>::clone() const
{
  //return new volume_vector_magnetic_harmonic_planar_linear_harmonic_rhs_1_1_current_1_1(this->i, this->j, this->m_offsetI, this->m_offsetJ, this->m_offsetPreviousTimeExt, this->m_offsetCouplingExt);
  return new volume_vector_magnetic_harmonic_planar_linear_harmonic_rhs_1_1_current_1_1(*this);
}


template <typename Scalar>
volume_vector_magnetic_harmonic_axisymmetric_linear_harmonic_rhs_1_1_current_1_1<Scalar>::volume_vector_magnetic_harmonic_axisymmetric_linear_harmonic_rhs_1_1_current_1_1(unsigned int i, unsigned int j, int offsetI, int offsetJ, int* offsetPreviousTimeExt, int* offsetCouplingExt)
: VectorFormVolAgros<Scalar>(i, offsetI, offsetJ, offsetPreviousTimeExt, offsetCouplingExt), j(j)
{
}

template <typename Scalar>
Scalar volume_vector_magnetic_harmonic_axisymmetric_linear_harmonic_rhs_1_1_current_1_1<Scalar>::value(int n, double *wt, Hermes::Hermes2D::Func<Scalar> *u_ext[], Hermes::Hermes2D::Func<double> *v,
  Hermes::Hermes2D::Geom<double> *e, Hermes::Hermes2D::Func<Scalar> **ext) const
{
  Scalar result = 0;

  for (int i = 0; i < n; i++)
  {
    result += wt[i] * (ext[8]->val[i] * v->val[i]);
  }
  return result;
}

template <typename Scalar>
Hermes::Ord volume_vector_magnetic_harmonic_axisymmetric_linear_harmonic_rhs_1_1_current_1_1<Scalar>::ord(int n, double *wt, Hermes::Hermes2D::Func<Hermes::Ord> *u_ext[], Hermes::Hermes2D::Func<Hermes::Ord> *v,
  Hermes::Hermes2D::Geom<Hermes::Ord> *e, Hermes::Hermes2D::Func<Hermes::Ord> **ext) const
{
  Hermes::Ord result(0);
  for (int i = 0; i < n; i++)
  {
    result += wt[i] * (ext[8]->val[i] * v->val[i]);
  }
  return result;
}

template <typename Scalar>
volume_vector_magnetic_harmonic_axisymmetric_linear_harmonic_rhs_1_1_current_1_1<Scalar>* volume_vector_magnetic_harmonic_axisymmetric_linear_harmonic_rhs_1_1_current_1_1<Scalar>::clone() const
{
  //return new volume_vector_magnetic_harmonic_axisymmetric_linear_harmonic_rhs_1_1_current_1_1(this->i, this->j, this->m_offsetI, this->m_offsetJ, this->m_offsetPreviousTimeExt, this->m_offsetCouplingExt);
  return new volume_vector_magnetic_harmonic_axisymmetric_linear_harmonic_rhs_1_1_current_1_1(*this);
}


template <typename Scalar>
volume_vector_magnetic_harmonic_planar_linear_harmonic_rhs_2_2_current_2_2<Scalar>::volume_vector_magnetic_harmonic_planar_linear_harmonic_rhs_2_2_current_2_2(unsigned int i, unsigned int j, int offsetI, int offsetJ, int* offsetPreviousTimeExt, int* offsetCouplingExt)
: VectorFormVolAgros<Scalar>(i, offsetI, offsetJ, offsetPreviousTimeExt, offsetCouplingExt), j(j)
{
}

template <typename Scalar>
Scalar volume_vector_magnetic_harmonic_planar_linear_harmonic_rhs_2_2_current_2_2<Scalar>::value(int n, double *wt, Hermes::Hermes2D::Func<Scalar> *u_ext[], Hermes::Hermes2D::Func<double> *v,
  Hermes::Hermes2D::Geom<double> *e, Hermes::Hermes2D::Func<Scalar> **ext) const
{
  Scalar result = 0;

  for (int i = 0; i < n; i++)
  {
    result += wt[i] * (ext[9]->val[i] * v->val[i]);
  }
  return result;
}

template <typename Scalar>
Hermes::Ord volume_vector_magnetic_harmonic_planar_linear_harmonic_rhs_2_2_current_2_2<Scalar>::ord(int n, double *wt, Hermes::Hermes2D::Func<Hermes::Ord> *u_ext[], Hermes::Hermes2D::Func<Hermes::Ord> *v,
  Hermes::Hermes2D::Geom<Hermes::Ord> *e, Hermes::Hermes2D::Func<Hermes::Ord> **ext) const
{
  Hermes::Ord result(0);
  for (int i = 0; i < n; i++)
  {
    result += wt[i] * (ext[9]->val[i] * v->val[i]);
  }
  return result;
}

template <typename Scalar>
volume_vector_magnetic_harmonic_planar_linear_harmonic_rhs_2_2_current_2_2<Scalar>* volume_vector_magnetic_harmonic_planar_linear_harmonic_rhs_2_2_current_2_2<Scalar>::clone() const
{
  //return new volume_vector_magnetic_harmonic_planar_linear_harmonic_rhs_2_2_current_2_2(this->i, this->j, this->m_offsetI, this->m_offsetJ, this->m_offsetPreviousTimeExt, this->m_offsetCouplingExt);
  return new volume_vector_magnetic_harmonic_planar_linear_harmonic_rhs_2_2_current_2_2(*this);
}


template <typename Scalar>
volume_vector_magnetic_harmonic_axisymmetric_linear_harmonic_rhs_2_2_current_2_2<Scalar>::volume_vector_magnetic_harmonic_axisymmetric_linear_harmonic_rhs_2_2_current_2_2(unsigned int i, unsigned int j, int offsetI, int offsetJ, int* offsetPreviousTimeExt, int* offsetCouplingExt)
: VectorFormVolAgros<Scalar>(i, offsetI, offsetJ, offsetPreviousTimeExt, offsetCouplingExt), j(j)
{
}

template <typename Scalar>
Scalar volume_vector_magnetic_harmonic_axisymmetric_linear_harmonic_rhs_2_2_current_2_2<Scalar>::value(int n, double *wt, Hermes::Hermes2D::Func<Scalar> *u_ext[], Hermes::Hermes2D::Func<double> *v,
  Hermes::Hermes2D::Geom<double> *e, Hermes::Hermes2D::Func<Scalar> **ext) const
{
  Scalar result = 0;

  for (int i = 0; i < n; i++)
  {
    result += wt[i] * (ext[9]->val[i] * v->val[i]);
  }
  return result;
}

template <typename Scalar>
Hermes::Ord volume_vector_magnetic_harmonic_axisymmetric_linear_harmonic_rhs_2_2_current_2_2<Scalar>::ord(int n, double *wt, Hermes::Hermes2D::Func<Hermes::Ord> *u_ext[], Hermes::Hermes2D::Func<Hermes::Ord> *v,
  Hermes::Hermes2D::Geom<Hermes::Ord> *e, Hermes::Hermes2D::Func<Hermes::Ord> **ext) const
{
  Hermes::Ord result(0);
  for (int i = 0; i < n; i++)
  {
    result += wt[i] * (ext[9]->val[i] * v->val[i]);
  }
  return result;
}

template <typename Scalar>
volume_vector_magnetic_harmonic_axisymmetric_linear_harmonic_rhs_2_2_current_2_2<Scalar>* volume_vector_magnetic_harmonic_axisymmetric_linear_harmonic_rhs_2_2_current_2_2<Scalar>::clone() const
{
  //return new volume_vector_magnetic_harmonic_axisymmetric_linear_harmonic_rhs_2_2_current_2_2(this->i, this->j, this->m_offsetI, this->m_offsetJ, this->m_offsetPreviousTimeExt, this->m_offsetCouplingExt);
  return new volume_vector_magnetic_harmonic_axisymmetric_linear_harmonic_rhs_2_2_current_2_2(*this);
}


template <typename Scalar>
volume_vector_magnetic_harmonic_planar_linear_harmonic_vec_current_3_3_3_3<Scalar>::volume_vector_magnetic_harmonic_planar_linear_harmonic_vec_current_3_3_3_3(unsigned int i, unsigned int j, int offsetI, int offsetJ, int* offsetPreviousTimeExt, int* offsetCouplingExt)
: VectorFormVolAgros<Scalar>(i, offsetI, offsetJ, offsetPreviousTimeExt, offsetCouplingExt), j(j)
{
}

template <typename Scalar>
Scalar volume_vector_magnetic_harmonic_planar_linear_harmonic_vec_current_3_3_3_3<Scalar>::value(int n, double *wt, Hermes::Hermes2D::Func<Scalar> *u_ext[], Hermes::Hermes2D::Func<double> *v,
  Hermes::Hermes2D::Geom<double> *e, Hermes::Hermes2D::Func<Scalar> **ext) const
{
  Scalar result = 0;

  for (int i = 0; i < n; i++)
  {
    result += wt[i] * (ext[10]->val[i] * (ext[12]->val[i] / this->markerVolume() - ext[9]->val[i]));
  }
  return result;
}

template <typename Scalar>
Hermes::Ord volume_vector_magnetic_harmonic_planar_linear_harmonic_vec_current_3_3_3_3<Scalar>::ord(int n, double *wt, Hermes::Hermes2D::Func<Hermes::Ord> *u_ext[], Hermes::Hermes2D::Func<Hermes::Ord> *v,
  Hermes::Hermes2D::Geom<Hermes::Ord> *e, Hermes::Hermes2D::Func<Hermes::Ord> **ext) const
{
  Hermes::Ord result(0);
  for (int i = 0; i < n; i++)
  {
    result += wt[i] * (ext[10]->val[i] * (ext[12]->val[i] / this->markerVolume() - ext[9]->val[i]));
  }
  return result;
}

template <typename Scalar>
volume_vector_magnetic_harmonic_planar_linear_harmonic_vec_current_3_3_3_3<Scalar>* volume_vector_magnetic_harmonic_planar_linear_harmonic_vec_current_3_3_3_3<Scalar>::clone() const
{
  //return new volume_vector_magnetic_harmonic_planar_linear_harmonic_vec_current_3_3_3_3(this->i, this->j, this->m_offsetI, this->m_offsetJ, this->m_offsetPreviousTimeExt, this->m_offsetCouplingExt);
  return new volume_vector_magnetic_harmonic_planar_linear_harmonic_vec_current_3_3_3_3(*this);
}


template <typename Scalar>
volume_vector_magnetic_harmonic_axisymmetric_linear_harmonic_vec_current_3_3_3_3<Scalar>::volume_vector_magnetic_harmonic_axisymmetric_linear_harmonic_vec_current_3_3_3_3(unsigned int i, unsigned int j, int offsetI, int offsetJ, int* offsetPreviousTimeExt, int* offsetCouplingExt)
: VectorFormVolAgros<Scalar>(i, offsetI, offsetJ, offsetPreviousTimeExt, offsetCouplingExt), j(j)
{
}

template <typename Scalar>
Scalar volume_vector_magnetic_harmonic_axisymmetric_linear_harmonic_vec_current_3_3_3_3<Scalar>::value(int n, double *wt, Hermes::Hermes2D::Func<Scalar> *u_ext[], Hermes::Hermes2D::Func<double> *v,
  Hermes::Hermes2D::Geom<double> *e, Hermes::Hermes2D::Func<Scalar> **ext) const
{
  Scalar result = 0;

  for (int i = 0; i < n; i++)
  {
    result += wt[i] * (ext[10]->val[i] * (ext[12]->val[i] / this->markerVolume() - ext[9]->val[i]));
  }
  return result;
}

template <typename Scalar>
Hermes::Ord volume_vector_magnetic_harmonic_axisymmetric_linear_harmonic_vec_current_3_3_3_3<Scalar>::ord(int n, double *wt, Hermes::Hermes2D::Func<Hermes::Ord> *u_ext[], Hermes::Hermes2D::Func<Hermes::Ord> *v,
  Hermes::Hermes2D::Geom<Hermes::Ord> *e, Hermes::Hermes2D::Func<Hermes::Ord> **ext) const
{
  Hermes::Ord result(0);
  for (int i = 0; i < n; i++)
  {
    result += wt[i] * (ext[10]->val[i] * (ext[12]->val[i] / this->markerVolume() - ext[9]->val[i]));
  }
  return result;
}

template <typename Scalar>
volume_vector_magnetic_harmonic_axisymmetric_linear_harmonic_vec_current_3_3_3_3<Scalar>* volume_vector_magnetic_harmonic_axisymmetric_linear_harmonic_vec_current_3_3_3_3<Scalar>::clone() const
{
  //return new volume_vector_magnetic_harmonic_axisymmetric_linear_harmonic_vec_current_3_3_3_3(this->i, this->j, this->m_offsetI, this->m_offsetJ, this->m_offsetPreviousTimeExt, this->m_offsetCouplingExt);
  return new volume_vector_magnetic_harmonic_axisymmetric_linear_harmonic_vec_current_3_3_3_3(*this);
}


template <typename Scalar>
volume_vector_magnetic_harmonic_planar_linear_harmonic_vec_current_4_4_4_4<Scalar>::volume_vector_magnetic_harmonic_planar_linear_harmonic_vec_current_4_4_4_4(unsigned int i, unsigned int j, int offsetI, int offsetJ, int* offsetPreviousTimeExt, int* offsetCouplingExt)
: VectorFormVolAgros<Scalar>(i, offsetI, offsetJ, offsetPreviousTimeExt, offsetCouplingExt), j(j)
{
}

template <typename Scalar>
Scalar volume_vector_magnetic_harmonic_planar_linear_harmonic_vec_current_4_4_4_4<Scalar>::value(int n, double *wt, Hermes::Hermes2D::Func<Scalar> *u_ext[], Hermes::Hermes2D::Func<double> *v,
  Hermes::Hermes2D::Geom<double> *e, Hermes::Hermes2D::Func<Scalar> **ext) const
{
  Scalar result = 0;

  for (int i = 0; i < n; i++)
  {
    result += wt[i] * (ext[10]->val[i] * (ext[11]->val[i] / this->markerVolume() - ext[8]->val[i]));
  }
  return result;
}

template <typename Scalar>
Hermes::Ord volume_vector_magnetic_harmonic_planar_linear_harmonic_vec_current_4_4_4_4<Scalar>::ord(int n, double *wt, Hermes::Hermes2D::Func<Hermes::Ord> *u_ext[], Hermes::Hermes2D::Func<Hermes::Ord> *v,
  Hermes::Hermes2D::Geom<Hermes::Ord> *e, Hermes::Hermes2D::Func<Hermes::Ord> **ext) const
{
  Hermes::Ord result(0);
  for (int i = 0; i < n; i++)
  {
    result += wt[i] * (ext[10]->val[i] * (ext[11]->val[i] / this->markerVolume() - ext[8]->val[i]));
  }
  return result;
}

template <typename Scalar>
volume_vector_magnetic_harmonic_planar_linear_harmonic_vec_current_4_4_4_4<Scalar>* volume_vector_magnetic_harmonic_planar_linear_harmonic_vec_current_4_4_4_4<Scalar>::clone() const
{
  //return new volume_vector_magnetic_harmonic_planar_linear_harmonic_vec_current_4_4_4_4(this->i, this->j, this->m_offsetI, this->m_offsetJ, this->m_offsetPreviousTimeExt, this->m_offsetCouplingExt);
  return new volume_vector_magnetic_harmonic_planar_linear_harmonic_vec_current_4_4_4_4(*this);
}




template <typename Scalar>
surface_vector_magnetic_harmonic_planar_linear_harmonic_current_1_1_1_magnetic_surface_current<Scalar>::surface_vector_magnetic_harmonic_planar_linear_harmonic_current_1_1_1_magnetic_surface_current(unsigned int i, unsigned int j, int offsetI, int offsetJ)

: VectorFormSurfAgros<Scalar>(i, offsetI, offsetJ), j(j)
{
}

template <typename Scalar>
Scalar surface_vector_magnetic_harmonic_planar_linear_harmonic_current_1_1_1_magnetic_surface_current<Scalar>::value(int n, double *wt, Hermes::Hermes2D::Func<Scalar> *u_ext[], Hermes::Hermes2D::Func<double> *v,
  Hermes::Hermes2D::Geom<double> *e, Hermes::Hermes2D::Func<Scalar> **ext) const
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
  {
    result += wt[i] * (-Jr*v->val[i]);
  }
  return result;
}

template <typename Scalar>
Hermes::Ord surface_vector_magnetic_harmonic_planar_linear_harmonic_current_1_1_1_magnetic_surface_current<Scalar>::ord(int n, double *wt, Hermes::Hermes2D::Func<Hermes::Ord> *u_ext[], Hermes::Hermes2D::Func<Hermes::Ord> *v,
  Hermes::Hermes2D::Geom<Hermes::Ord> *e, Hermes::Hermes2D::Func<Hermes::Ord> **ext) const
{
  Hermes::Ord result(0);
  for (int i = 0; i < n; i++)
  {
    result += wt[i] * (-Jr*v->val[i]);
  }
  return result;

}

template <typename Scalar>
surface_vector_magnetic_harmonic_planar_linear_harmonic_current_1_1_1_magnetic_surface_current<Scalar>* surface_vector_magnetic_harmonic_planar_linear_harmonic_current_1_1_1_magnetic_surface_current<Scalar>::clone() const
{
  //return new surface_vector_magnetic_harmonic_planar_linear_harmonic_current_1_1_1_magnetic_surface_current(this->i, this->j, this->m_offsetI, this->m_offsetJ);
  return new surface_vector_magnetic_harmonic_planar_linear_harmonic_current_1_1_1_magnetic_surface_current(*this);
}

template <typename Scalar>
surface_vector_magnetic_harmonic_axisymmetric_linear_harmonic_current_1_1_1_magnetic_surface_current<Scalar>::surface_vector_magnetic_harmonic_axisymmetric_linear_harmonic_current_1_1_1_magnetic_surface_current(unsigned int i, unsigned int j, int offsetI, int offsetJ)
: VectorFormSurfAgros<Scalar>(i, offsetI, offsetJ), j(j)
{
}

template <typename Scalar>
Scalar surface_vector_magnetic_harmonic_axisymmetric_linear_harmonic_current_1_1_1_magnetic_surface_current<Scalar>::value(int n, double *wt, Hermes::Hermes2D::Func<Scalar> *u_ext[], Hermes::Hermes2D::Func<double> *v,
  Hermes::Hermes2D::Geom<double> *e, Hermes::Hermes2D::Func<Scalar> **ext) const
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
  {
    result += wt[i] * (-Jr*v->val[i]);
  }
  return result;
}

template <typename Scalar>
Hermes::Ord surface_vector_magnetic_harmonic_axisymmetric_linear_harmonic_current_1_1_1_magnetic_surface_current<Scalar>::ord(int n, double *wt, Hermes::Hermes2D::Func<Hermes::Ord> *u_ext[], Hermes::Hermes2D::Func<Hermes::Ord> *v,
  Hermes::Hermes2D::Geom<Hermes::Ord> *e, Hermes::Hermes2D::Func<Hermes::Ord> **ext) const
{
  Hermes::Ord result(0);
  for (int i = 0; i < n; i++)
  {
    result += wt[i] * (-Jr*v->val[i]);
  }
  return result;

}

template <typename Scalar>
surface_vector_magnetic_harmonic_axisymmetric_linear_harmonic_current_1_1_1_magnetic_surface_current<Scalar>* surface_vector_magnetic_harmonic_axisymmetric_linear_harmonic_current_1_1_1_magnetic_surface_current<Scalar>::clone() const
{
  //return new surface_vector_magnetic_harmonic_axisymmetric_linear_harmonic_current_1_1_1_magnetic_surface_current(this->i, this->j, this->m_offsetI, this->m_offsetJ);
  return new surface_vector_magnetic_harmonic_axisymmetric_linear_harmonic_current_1_1_1_magnetic_surface_current(*this);
}

template <typename Scalar>
surface_vector_magnetic_harmonic_planar_linear_harmonic_current_2_2_2_magnetic_surface_current<Scalar>::surface_vector_magnetic_harmonic_planar_linear_harmonic_current_2_2_2_magnetic_surface_current(unsigned int i, unsigned int j, int offsetI, int offsetJ)

: VectorFormSurfAgros<Scalar>(i, offsetI, offsetJ), j(j)
{
}

template <typename Scalar>
Scalar surface_vector_magnetic_harmonic_planar_linear_harmonic_current_2_2_2_magnetic_surface_current<Scalar>::value(int n, double *wt, Hermes::Hermes2D::Func<Scalar> *u_ext[], Hermes::Hermes2D::Func<double> *v,
  Hermes::Hermes2D::Geom<double> *e, Hermes::Hermes2D::Func<Scalar> **ext) const
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
  {
    result += wt[i] * (-Ji*v->val[i]);
  }
  return result;
}

template <typename Scalar>
Hermes::Ord surface_vector_magnetic_harmonic_planar_linear_harmonic_current_2_2_2_magnetic_surface_current<Scalar>::ord(int n, double *wt, Hermes::Hermes2D::Func<Hermes::Ord> *u_ext[], Hermes::Hermes2D::Func<Hermes::Ord> *v,
  Hermes::Hermes2D::Geom<Hermes::Ord> *e, Hermes::Hermes2D::Func<Hermes::Ord> **ext) const
{
  Hermes::Ord result(0);
  for (int i = 0; i < n; i++)
  {
    result += wt[i] * (-Ji*v->val[i]);
  }
  return result;

}

template <typename Scalar>
surface_vector_magnetic_harmonic_planar_linear_harmonic_current_2_2_2_magnetic_surface_current<Scalar>* surface_vector_magnetic_harmonic_planar_linear_harmonic_current_2_2_2_magnetic_surface_current<Scalar>::clone() const
{
  //return new surface_vector_magnetic_harmonic_planar_linear_harmonic_current_2_2_2_magnetic_surface_current(this->i, this->j, this->m_offsetI, this->m_offsetJ);
  return new surface_vector_magnetic_harmonic_planar_linear_harmonic_current_2_2_2_magnetic_surface_current(*this);
}

template <typename Scalar>
exact_magnetic_harmonic_planar_linear_harmonic_essential_1_1_0_magnetic_potential<Scalar>::exact_magnetic_harmonic_planar_linear_harmonic_essential_1_1_0_magnetic_potential(Hermes::Hermes2D::MeshSharedPtr mesh)
: ExactSolutionScalarAgros<Scalar>(mesh)
{
}

template <typename Scalar>
Scalar exact_magnetic_harmonic_planar_linear_harmonic_essential_1_1_0_magnetic_potential<Scalar>::value(double x, double y) const
{
  Scalar result = Ar;
  return result;
}

template <typename Scalar>
void exact_magnetic_harmonic_planar_linear_harmonic_essential_1_1_0_magnetic_potential<Scalar>::derivatives(double x, double y, Scalar& dx, Scalar& dy) const
{

}

template <typename Scalar>
exact_magnetic_harmonic_planar_linear_harmonic_essential_2_2_0_magnetic_potential<Scalar>::exact_magnetic_harmonic_planar_linear_harmonic_essential_2_2_0_magnetic_potential(Hermes::Hermes2D::MeshSharedPtr mesh)
: ExactSolutionScalarAgros<Scalar>(mesh)
{
}

template <typename Scalar>
Scalar exact_magnetic_harmonic_planar_linear_harmonic_essential_2_2_0_magnetic_potential<Scalar>::value(double x, double y) const
{
  Scalar result = Ai;
  return result;
}

template <typename Scalar>
void exact_magnetic_harmonic_planar_linear_harmonic_essential_2_2_0_magnetic_potential<Scalar>::derivatives(double x, double y, Scalar& dx, Scalar& dy) const
{

}


template class volume_matrix_magnetic_harmonic_planar_linear_harmonic_laplace_1_1_1_1<double>;
template class volume_matrix_magnetic_harmonic_planar_linear_harmonic_laplace_2_2_2_2<double>;
template class volume_matrix_magnetic_harmonic_planar_linear_harmonic_velocity_1_1_1_1<double>;
template class volume_matrix_magnetic_harmonic_planar_linear_harmonic_velocity_2_2_2_2<double>;
template class volume_matrix_magnetic_harmonic_planar_linear_harmonic_gamma_1_2_1_2<double>;
template class volume_matrix_magnetic_harmonic_planar_linear_harmonic_gamma_2_1_2_1<double>;
template class volume_matrix_magnetic_harmonic_planar_linear_harmonic_current_1_4_1_4<double>;
template class volume_matrix_magnetic_harmonic_planar_linear_harmonic_current_2_3_2_3<double>;
template class volume_matrix_magnetic_harmonic_planar_linear_harmonic_current_3_1_3_1<double>;
template class volume_matrix_magnetic_harmonic_planar_linear_harmonic_current_4_2_4_2<double>;
template class volume_matrix_magnetic_harmonic_planar_linear_harmonic_current_3_3_3_3<double>;
template class volume_matrix_magnetic_harmonic_planar_linear_harmonic_current_4_4_4_4<double>;
template class volume_vector_magnetic_harmonic_planar_linear_harmonic_rhs_1_1_remanence_1_1<double>;
template class volume_vector_magnetic_harmonic_planar_linear_harmonic_rhs_1_1_current_1_1<double>;
template class volume_vector_magnetic_harmonic_planar_linear_harmonic_rhs_2_2_current_2_2<double>;
template class volume_vector_magnetic_harmonic_planar_linear_harmonic_vec_current_3_3_3_3<double>;
template class volume_vector_magnetic_harmonic_planar_linear_harmonic_vec_current_4_4_4_4<double>;
template class exact_magnetic_harmonic_planar_linear_harmonic_essential_1_1_0_magnetic_potential<double>;
template class exact_magnetic_harmonic_planar_linear_harmonic_essential_2_2_0_magnetic_potential<double>;
template class surface_vector_magnetic_harmonic_planar_linear_harmonic_current_1_1_1_magnetic_surface_current<double>;
template class surface_vector_magnetic_harmonic_planar_linear_harmonic_current_2_2_2_magnetic_surface_current<double>;
