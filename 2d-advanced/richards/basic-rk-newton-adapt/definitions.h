#include "hermes2d.h"
#include "../constitutive.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::WeakFormsH1;
using namespace Hermes::Hermes2D::Views;
using namespace Hermes::Hermes2D::RefinementSelectors;

/* Custom non-constant Dirichlet condition */

class CustomEssentialBCNonConst : public EssentialBoundaryCondition <double>
{
public:
  CustomEssentialBCNonConst(std::vector<std::string>(markers))
    : EssentialBoundaryCondition<double>(markers)
  {
  }

  virtual EssentialBCValueType get_value_type() const;

  virtual double value(double x, double y) const;
};

/* Weak forms */

class CustomWeakFormRichardsRK : public WeakForm <double>
{
public:
  CustomWeakFormRichardsRK(ConstitutiveRelations* constitutive);

private:

  class CustomJacobianFormVol : public MatrixFormVol <double>
  {
  public:
    CustomJacobianFormVol(int i, int j, ConstitutiveRelations* constitutive)
      : MatrixFormVol<double>(i, j), constitutive(constitutive)
    {
    }

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u,
      Func<double> *v, GeomVol<double> *e, Func<double>* *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u,
      Func<Ord> *v, GeomVol<Ord> *e, Func<Ord>* *ext) const;

    virtual MatrixFormVol<double>* clone() const;
    ConstitutiveRelations* constitutive;
  };

  class CustomResidualFormVol : public VectorFormVol <double>
  {
  public:
    CustomResidualFormVol(int i, ConstitutiveRelations* constitutive)
      : VectorFormVol<double>(i), constitutive(constitutive)
    {
    }

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, GeomVol<double> *e,
      Func<double>* *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, GeomVol<Ord> *e,
      Func<Ord>* *ext) const;

    virtual VectorFormVol<double>* clone() const;
    ConstitutiveRelations* constitutive;
  };
  ConstitutiveRelations* constitutive;

  WeakForm<double>* clone() const
  {
    CustomWeakFormRichardsRK* wf = new CustomWeakFormRichardsRK(*this);
    wf->constitutive = this->constitutive;
    return wf;
  }
};
