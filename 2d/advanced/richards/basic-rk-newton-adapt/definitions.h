#include "hermes2d.h"

/* Namespaces used */

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;
using namespace Hermes::Hermes2D::RefinementSelectors;

// K (Gardner).
double K(double h);

// dK/dh (Gardner).
double dKdh(double h);

// ddK/dhh (Gardner).
double ddKdhh(double h);

// C (Gardner).
double C(double h);

// dC/dh (Gardner).
double dCdh(double h);

// ddC/dhh (Gardner).
double ddCdhh(double h);

/* Custom non-constant Dirichlet condition */

class CustomEssentialBCNonConst : public EssentialBoundaryCondition<double> 
{
public:
  CustomEssentialBCNonConst(Hermes::vector<std::string>(markers))       
        : EssentialBoundaryCondition<double>(markers) {};

  virtual EssentialBCValueType get_value_type() const;

  virtual double value(double x, double y, double n_x, double n_y, 
                       double t_x, double t_y) const;
};

/* Weak forms */

class CustomWeakFormRichardsRK : public WeakForm<double>
{
public:
  CustomWeakFormRichardsRK();

private:

  class CustomJacobianFormVol : public MatrixFormVol<double>
  {
  public:
    CustomJacobianFormVol(int i, int j) 
          : MatrixFormVol<double>(i, j) 
    {
    }

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
                         Func<double> *v, Geom<double> *e, ExtData<double> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, 
                    Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const;

    virtual MatrixFormVol<double>* clone();
  };

  class CustomResidualFormVol : public VectorFormVol<double>
  {
  public:
    CustomResidualFormVol(int i)
          : VectorFormVol<double>(i) 
    {
    }

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e,
                         ExtData<double> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, 
                    ExtData<Ord> *ext) const;

    virtual VectorFormVol<double>* clone();
  };
};

