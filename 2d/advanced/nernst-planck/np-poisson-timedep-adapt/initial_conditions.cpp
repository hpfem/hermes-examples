
#include "header.h"


class InitialSolutionVoltage : public ExactSolutionScalar<double>
{
public:
  InitialSolutionVoltage(Mesh* mesh) : ExactSolutionScalar<double>(mesh) {};

  virtual double value (double x, double y) const {
    return 0.0;
  };

  virtual void derivatives (double x, double y, double& dx, double& dy) const {
    dx = 0;
    dy = 0;
  };

  virtual Ord ord(Ord x, Ord y) const {
    return Ord(0);
  }
};

class InitialSolutionConcentration : public ExactSolutionScalar<double>
{
public:
  InitialSolutionConcentration(Mesh* mesh, double C0) : ExactSolutionScalar<double>(mesh), C0(C0) {};

  virtual double value (double x, double y) const {
    return C0;
  };

  virtual void derivatives (double x, double y, double& dx, double& dy) const {
    dx = 0;
    dy = 0;
  };

  virtual Ord ord(Ord x, Ord y) const {
    return Ord(0);
  }

  double C0;
};
