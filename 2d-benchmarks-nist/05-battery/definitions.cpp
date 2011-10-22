#include "definitions.h"

template<typename Real, typename Scalar>
Scalar CustomMatrixFormVol::matrix_form(int n, double *wt, Func<Scalar> *u_ext[],     
    Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const
{
      Scalar result = Scalar(0);
      double p = 0, q = 0;
      // integration order calculation.
      if(e->elem_marker == -9999) p = q = 1;
      else {
        if(e->elem_marker == mesh->get_element_markers_conversion().get_internal_marker("e1").marker) 
        {
           p = static_cast<CustomWeakFormPoisson*>(wf)->p_1;
           q = static_cast<CustomWeakFormPoisson*>(wf)->q_1;
        }
        if(e->elem_marker == mesh->get_element_markers_conversion().get_internal_marker("e2").marker) 
        {
           p = static_cast<CustomWeakFormPoisson*>(wf)->p_2;
           q = static_cast<CustomWeakFormPoisson*>(wf)->q_2;
        }
        if(e->elem_marker == mesh->get_element_markers_conversion().get_internal_marker("e3").marker) 
        {
           p = static_cast<CustomWeakFormPoisson*>(wf)->p_3;
           q = static_cast<CustomWeakFormPoisson*>(wf)->q_3;
        }
        if(e->elem_marker == mesh->get_element_markers_conversion().get_internal_marker("e4").marker) 
        {
           p = static_cast<CustomWeakFormPoisson*>(wf)->p_4;
           q = static_cast<CustomWeakFormPoisson*>(wf)->q_4;
        }
        if(e->elem_marker == mesh->get_element_markers_conversion().get_internal_marker("e5").marker) 
        {
           p = static_cast<CustomWeakFormPoisson*>(wf)->p_5;
           q = static_cast<CustomWeakFormPoisson*>(wf)->q_5;
        }
      }

      for (int i = 0; i < n; i++)
        result += wt[i] * (p * u->dx[i] * v->dx[i] + q * u->dy[i] * v->dy[i]);
      return result;
}

double CustomMatrixFormVol::value(int n, double *wt, Func<double> *u_ext[],
    Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<double> *ext) const 
{
  return matrix_form<double, double>(n, wt, u_ext, u, v, e, ext);
}

Ord CustomMatrixFormVol::ord(int n, double *wt, Func<Ord> *u_ext[],
    Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const 
{
  return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
}

template<typename Real, typename Scalar>
Scalar CustomVectorFormVol::vector_form(int n, double *wt, Func<Scalar> *u_ext[],
    Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const
{
      double f = 0;
      Scalar result = Scalar(0);
      if(e->elem_marker == -9999)
        f = 1;
      else {
        if(e->elem_marker == mesh->get_element_markers_conversion().get_internal_marker("e1").marker)
            f = static_cast<CustomWeakFormPoisson*>(wf)->f_1;
        if(e->elem_marker == mesh->get_element_markers_conversion().get_internal_marker("e2").marker)
            f = static_cast<CustomWeakFormPoisson*>(wf)->f_2;
        if(e->elem_marker == mesh->get_element_markers_conversion().get_internal_marker("e3").marker)
            f = static_cast<CustomWeakFormPoisson*>(wf)->f_3;
        if(e->elem_marker == mesh->get_element_markers_conversion().get_internal_marker("e4").marker)
            f = static_cast<CustomWeakFormPoisson*>(wf)->f_4;
        if(e->elem_marker == mesh->get_element_markers_conversion().get_internal_marker("e5").marker)
            f = static_cast<CustomWeakFormPoisson*>(wf)->f_5;
       }

      double p = 0, q = 0;
      // integration order calculation.
      if(e->elem_marker == -9999) p = q = 1;
      else {
        if(e->elem_marker == mesh->get_element_markers_conversion().get_internal_marker("e1").marker) 
        {
           p = static_cast<CustomWeakFormPoisson*>(wf)->p_1;
           q = static_cast<CustomWeakFormPoisson*>(wf)->q_1;
        }
        if(e->elem_marker == mesh->get_element_markers_conversion().get_internal_marker("e2").marker) 
        {
           p = static_cast<CustomWeakFormPoisson*>(wf)->p_2;
           q = static_cast<CustomWeakFormPoisson*>(wf)->q_2;
        }
        if(e->elem_marker == mesh->get_element_markers_conversion().get_internal_marker("e3").marker) 
        {
           p = static_cast<CustomWeakFormPoisson*>(wf)->p_3;
           q = static_cast<CustomWeakFormPoisson*>(wf)->q_3;
        }
        if(e->elem_marker == mesh->get_element_markers_conversion().get_internal_marker("e4").marker) 
        {
           p = static_cast<CustomWeakFormPoisson*>(wf)->p_4;
           q = static_cast<CustomWeakFormPoisson*>(wf)->q_4;
        }
        if(e->elem_marker == mesh->get_element_markers_conversion().get_internal_marker("e5").marker) 
        {
           p = static_cast<CustomWeakFormPoisson*>(wf)->p_5;
           q = static_cast<CustomWeakFormPoisson*>(wf)->q_5;
        }
      }

      for (int i = 0; i < n; i++)
        result += wt[i] * (p * u_ext[0]->dx[i] * v->dx[i] + q * u_ext[0]->dy[i] * v->dy[i]);

      result = result - f * int_v<Real>(n, wt, v);
      return result;
}

double CustomVectorFormVol::value(int n, double *wt, Func<double> *u_ext[],
    Func<double> *v, Geom<double> *e, ExtData<double> *ext) const 
{
  return vector_form<double, double>(n, wt, u_ext, v, e, ext);
}

Ord CustomVectorFormVol::ord(int n, double *wt, Func<Ord> *u_ext[],
    Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const 
{
  return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
}

template<typename Real, typename Scalar>
Scalar CustomMatrixFormSurf::matrix_form(int n, double *wt, Func<Scalar> *u_ext[],     
    Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const
{
      Scalar result = Scalar(0);
      for (int i = 0; i < n; i++) 
      {
        Real x = e->x[i];
        Real y = e->y[i];
        double p = 0, q = 0, c = 1.;
        if(this->areas[0] == static_cast<CustomWeakFormPoisson*>(wf)->bdy_left) 
        {
          if (x == 0.0) 
          {
            if ((y >= 0.0 && y <= 0.8)||(y >= 23.2 && y <= 24.0)) 
            {
              p = static_cast<CustomWeakFormPoisson*>(wf)->p_1; 
              q = static_cast<CustomWeakFormPoisson*>(wf)->q_1;
            }
            if ((y >= 1.6 && y <= 3.6)||(y >= 18.8 && y <= 21.2)) 
            {
              p = static_cast<CustomWeakFormPoisson*>(wf)->p_2; 
              q = static_cast<CustomWeakFormPoisson*>(wf)->q_2;
            }
            if (y >= 3.6 && y <= 18.8) 
            {
              p = static_cast<CustomWeakFormPoisson*>(wf)->p_3; 
              q = static_cast<CustomWeakFormPoisson*>(wf)->q_3;
            }
            if ((y >= 0.8 && y <= 1.6)||(y >= 21.2 && y <= 23.2)) 
            {
              p = static_cast<CustomWeakFormPoisson*>(wf)->p_5; 
              q = static_cast<CustomWeakFormPoisson*>(wf)->q_5;
            }
          }
          c = static_cast<CustomWeakFormPoisson*>(wf)->c_left;
        }
        if(this->areas[0] == static_cast<CustomWeakFormPoisson*>(wf)->bdy_right) 
        {
          p = static_cast<CustomWeakFormPoisson*>(wf)->p_1; 
          q = static_cast<CustomWeakFormPoisson*>(wf)->q_1;
          c = static_cast<CustomWeakFormPoisson*>(wf)->c_right;
        }

        if(this->areas[0] == static_cast<CustomWeakFormPoisson*>(wf)->bdy_bottom) 
        {
          p = static_cast<CustomWeakFormPoisson*>(wf)->p_1; 
          q = static_cast<CustomWeakFormPoisson*>(wf)->q_1;
          c = static_cast<CustomWeakFormPoisson*>(wf)->c_bottom;
        }

        if(this->areas[0] == static_cast<CustomWeakFormPoisson*>(wf)->bdy_top) 
        {
          p = static_cast<CustomWeakFormPoisson*>(wf)->p_1; 
          q = static_cast<CustomWeakFormPoisson*>(wf)->q_1;
          c = static_cast<CustomWeakFormPoisson*>(wf)->c_top;

        }
        result += wt[i] * (//p * u->dx[i] * v->val[i] - q * u->dy[i] * v->val[i] 
                  + c * u->val[i] * v->val[i]);
      }
      return result;
}

double CustomMatrixFormSurf::value(int n, double *wt, Func<double> *u_ext[],
    Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<double> *ext) const 
{
  return matrix_form<double, double>(n, wt, u_ext, u, v, e, ext);
}

Ord CustomMatrixFormSurf::ord(int n, double *wt, Func<Ord> *u_ext[],
    Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const 
{
  return Ord(4.0);
}

template<typename Real, typename Scalar>
Scalar CustomVectorFormSurf::vector_form(int n, double *wt, Func<Scalar> *u_ext[],
    Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const
{
      Scalar result = Scalar(0);
      double g = 1.0;

      if(this->areas[0] == static_cast<CustomWeakFormPoisson*>(wf)->bdy_left) 
      {
        g = static_cast<CustomWeakFormPoisson*>(wf)->g_n_left; 
      }
      if(this->areas[0] == static_cast<CustomWeakFormPoisson*>(wf)->bdy_right) 
      {
        g = static_cast<CustomWeakFormPoisson*>(wf)->g_n_right;
      }

      if(this->areas[0] == static_cast<CustomWeakFormPoisson*>(wf)->bdy_bottom) 
      {
        g = static_cast<CustomWeakFormPoisson*>(wf)->g_n_bottom;
      }

      if(this->areas[0] == static_cast<CustomWeakFormPoisson*>(wf)->bdy_top) 
      {
        g = static_cast<CustomWeakFormPoisson*>(wf)->g_n_top;
      }

      for (int i = 0; i < n; i++) 
      {
        Real x = e->x[i];
        Real y = e->y[i];
        double p = 0, q = 0, c = 1.;
        if(this->areas[0] == static_cast<CustomWeakFormPoisson*>(wf)->bdy_left) 
        {
          if (x == 0.0) 
          {
            if ((y >= 0.0 && y <= 0.8)||(y >= 23.2 && y <= 24.0)) 
            {
              p = static_cast<CustomWeakFormPoisson*>(wf)->p_1; 
              q = static_cast<CustomWeakFormPoisson*>(wf)->q_1;
            }
            if ((y >= 1.6 && y <= 3.6)||(y >= 18.8 && y <= 21.2)) 
            {
              p = static_cast<CustomWeakFormPoisson*>(wf)->p_2; 
              q = static_cast<CustomWeakFormPoisson*>(wf)->q_2;
            }
            if (y >= 3.6 && y <= 18.8) 
            {
              p = static_cast<CustomWeakFormPoisson*>(wf)->p_3; 
              q = static_cast<CustomWeakFormPoisson*>(wf)->q_3;
            }
            if ((y >= 0.8 && y <= 1.6)||(y >= 21.2 && y <= 23.2)) 
            {
              p = static_cast<CustomWeakFormPoisson*>(wf)->p_5; 
              q = static_cast<CustomWeakFormPoisson*>(wf)->q_5;
            }
          }
          c = static_cast<CustomWeakFormPoisson*>(wf)->c_left;
        }
        if(this->areas[0] == static_cast<CustomWeakFormPoisson*>(wf)->bdy_right) 
        {
          p = static_cast<CustomWeakFormPoisson*>(wf)->p_1; 
          q = static_cast<CustomWeakFormPoisson*>(wf)->q_1;
          c = static_cast<CustomWeakFormPoisson*>(wf)->c_right;
        }

        if(this->areas[0] == static_cast<CustomWeakFormPoisson*>(wf)->bdy_bottom) 
        {
          p = static_cast<CustomWeakFormPoisson*>(wf)->p_1; 
          q = static_cast<CustomWeakFormPoisson*>(wf)->q_1;
          c = static_cast<CustomWeakFormPoisson*>(wf)->c_bottom;
        }

        if(this->areas[0] == static_cast<CustomWeakFormPoisson*>(wf)->bdy_top) 
        {
          p = static_cast<CustomWeakFormPoisson*>(wf)->p_1; 
          q = static_cast<CustomWeakFormPoisson*>(wf)->q_1;
          c = static_cast<CustomWeakFormPoisson*>(wf)->c_top;
        }
        result += wt[i] * (//p * u_ext[0]->dx[i] * v->val[i] - q * u_ext[0]->dy[i] * v->val[i] 
                  + c * u_ext[0]->val[i] * v->val[i]);
      }

      result = result - g * int_v<Real>(n, wt, v);

      return result;
}

double CustomVectorFormSurf::value(int n, double *wt, Func<double> *u_ext[],
    Func<double> *v, Geom<double> *e, ExtData<double> *ext) const 
{
  return vector_form<double, double>(n, wt, u_ext, v, e, ext);
}

Ord CustomVectorFormSurf::ord(int n, double *wt, Func<Ord> *u_ext[],
    Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const 
{
  return Ord(4.0);
}

CustomWeakFormPoisson::CustomWeakFormPoisson(std::string omega_1, std::string omega_2, 
                                             std::string omega_3, std::string omega_4, 
                                             std::string omega_5, std::string bdy_left, 
                                             std::string bdy_top, std::string bdy_right, 
                                             std::string bdy_bottom, Mesh* mesh) : WeakForm<double>(1),
  
  omega_1(omega_1), omega_2(omega_2), omega_3(omega_3), 
  omega_4(omega_4), omega_5(omega_5), mesh(mesh),

  p_1(25.0),
  p_2(7.0),
  p_3(5.0),
  p_4(0.2),
  p_5(0.05),

  q_1(25.0),
  q_2(0.8),
  q_3(0.0001),
  q_4(0.2),
  q_5(0.05),

  f_1(0.0),
  f_2(1.0),
  f_3(1.0),
  f_4(0.0),
  f_5(0.0),

  bdy_left(bdy_left), 
  bdy_top(bdy_top), 
  bdy_right(bdy_right), 
  bdy_bottom(bdy_bottom),

  c_left(0.0),
  c_top(1.0),
  c_right(2.0),
  c_bottom(3.0),

  g_n_left(0.0),
  g_n_top(3.0),
  g_n_right(2.0),
  g_n_bottom(1.0)

{
    add_matrix_form(new CustomMatrixFormVol(0, 0, mesh));
    add_vector_form(new CustomVectorFormVol(0, mesh));

    add_matrix_form_surf(new CustomMatrixFormSurf(0, 0, bdy_bottom));
    add_matrix_form_surf(new CustomMatrixFormSurf(0, 0, bdy_right));
    add_matrix_form_surf(new CustomMatrixFormSurf(0, 0, bdy_top));

    add_vector_form_surf(new CustomVectorFormSurf(0, bdy_bottom));
    add_vector_form_surf(new CustomVectorFormSurf(0, bdy_top));
    add_vector_form_surf(new CustomVectorFormSurf(0, bdy_left));
    add_vector_form_surf(new CustomVectorFormSurf(0, bdy_right));

}
