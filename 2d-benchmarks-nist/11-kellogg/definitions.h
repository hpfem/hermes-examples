#include "hermes2d.h"
#include "../NIST-util.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;

/* Exact solution */

class CustomExactSolution : public ExactSolutionScalar<double>
{
public:
  CustomExactSolution(MeshSharedPtr mesh, double sigma, double tau, double rho)
      : ExactSolutionScalar<double>(mesh), sigma(sigma), tau(tau), rho(rho) {};

  virtual double value(double x, double y) const;

  virtual void derivatives(double x, double y, double& dx, double& dy) const;

  virtual Ord ord (double x, double y) const; 

  MeshFunction<double>* clone() const { return new CustomExactSolution(mesh, sigma, tau, rho); }

  double sigma;
  double tau;
  double rho;
};

/* Weak forms */

class CustomWeakFormPoisson : public WeakForm<double>
{
public:
  CustomWeakFormPoisson(std::string area_1, double r, std::string area_2);
};
