#include "definitions.h"

double CustomWeakForm::Jacobian::value(int n, double* wt, 
                                       Func< double >* u_ext[], Func< double >* u, Func< double >* v, 
                                       Geom< double >* e, ExtData< double >* ext) const
{
  double result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]);
  
  return result;
}

Ord CustomWeakForm::Jacobian::ord(int n, double* wt, 
                                  Func< Ord >* u_ext[], Func< Ord >* u, Func< Ord >* v, 
                                  Geom< Ord >* e, ExtData< Ord >* ext) const
{ 
  return u->dx[0] * v->dx[0] + u->dy[0] * v->dy[0];
}

double CustomWeakForm::Residual::value(int n, double* wt, Func< double >* u_ext[], Func< double >* v,
                                       Geom< double >* e, ExtData< double >* ext) const
{
  double result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * ( u_ext[0]->dx[i] * v->dx[i] + u_ext[0]->dy[i] * v->dy[i] + rhs->value(e->x[i], e->y[i]) * v->val[i] );
  
  return result;
}

Ord CustomWeakForm::Residual::ord(int n, double* wt, Func< Ord >* u_ext[], Func< Ord >* v, 
                                  Geom< Ord >* e, ExtData< Ord >* ext) const
{ 
  return u_ext[0]->dx[0] * v->dx[0] + u_ext[0]->dy[0] * v->dy[0] + rhs->value(e->x[0], e->y[0]) * v->val[0];
}

double CustomRightHandSide::value(double x, double y) const
{
  double a = pow(x - x_loc, 2);
  double b = pow(y - y_loc, 2);
  double c = sqrt(a + b);
  double d = ((alpha*x - (alpha * x_loc)) * (2*x - (2 * x_loc)));
  double e = ((alpha*y - (alpha * y_loc)) * (2*y - (2 * y_loc)));
  double f = (pow(alpha*c - (alpha * r_zero), 2) + 1.0);
  double g = (alpha * c - (alpha * r_zero));
  
  return +( ((alpha/(c * f)) - (d/(2 * pow(a + b, 1.5) * f))
          - ((alpha * d * g)/((a + b) * pow(f, 2))) + (alpha/(c * f)) 
          - (e/(2 * pow(a + b, 1.5) * f))
          - ((alpha * e * g)/((a + b) * pow(f, 2)))));
}

void CustomExactSolution::derivatives(double x, double y, double& dx, double& dy) const
{
  double a = pow(x - x_loc, 2);
  double b = pow(y - y_loc, 2);
  double c = sqrt(a + b);
  double d = (alpha*x - (alpha * x_loc));
  double e = (alpha*y - (alpha * y_loc));
  double f = (pow(alpha*c - (alpha * r_zero), 2) + 1.0);
  
  dx = (d/(c * f));
  dy = (e/(c * f));
}

double ResidualErrorForm::value(int n, double* wt, 
                                Func< double >* u_ext[], Func< double >* u, 
                                Geom< double >* e, ExtData< double >* ext) const
{
#ifdef H2D_SECOND_DERIVATIVES_ENABLED
  double result = 0.;

  for (int i = 0; i < n; i++)
    result += wt[i] * Hermes::sqr( -rhs->value(e->x[i], e->y[i]) + u->laplace[i] );

  return result * Hermes::sqr(e->diam);
#else
  error("Define H2D_SECOND_DERIVATIVES_ENABLED in hermes2d_common_defs.h"
        "if you want to use second derivatives in weak forms.");
#endif
}

Ord ResidualErrorForm::ord(int n, double* wt, 
                           Func< Ord >* u_ext[], Func< Ord >* u, 
                           Geom< Ord >* e, ExtData< Ord >* ext) const
{
#ifdef H2D_SECOND_DERIVATIVES_ENABLED
  return sqr( -rhs->value(e->x[0], e->y[0]) + u->laplace[0] );
#else
  error("Define H2D_SECOND_DERIVATIVES_ENABLED in hermes2d_common_defs.h"
        "if you want to use second derivatives in weak forms.");
#endif
}

void ConvergenceTable::save(const char* filename) const
{
  if (num_columns() == 0) error("No data columns defined.");
  for (unsigned int i = 1; i < num_columns(); i++)
    if (columns[i].data.size() != num_rows())
      error("Incompatible column sizes."); 
  
  FILE* f = fopen(filename, "w");
  if (f == NULL) error("Error writing to %s.", filename);
    
  for (unsigned int i = 0; i < num_columns(); i++)
    fprintf(f, columns[i].label.c_str());
  fprintf(f, "\n");
  
  for (int j = 0; j < num_rows(); j++)
  {
    for (unsigned int i = 0; i < num_columns(); i++)
    {
      if (columns[i].data[j].data_type == Column::Entry::INT)
        fprintf(f, columns[i].format.c_str(), columns[i].data[j].ivalue);
      else if (columns[i].data[j].data_type == Column::Entry::DOUBLE)
        fprintf(f, columns[i].format.c_str(), columns[i].data[j].dvalue);
    }
    fprintf(f, "\n");
  }
  
  fclose(f);
  verbose("Convergence table saved to file '%s'.", filename);
}

std::string itos(const unsigned int i)
{
  std::stringstream s;
  s << i;
  return s.str();
}