class InitialSolutionEulerDensity : public ExactSolutionScalar<double>
{
public:
  InitialSolutionEulerDensity(Mesh* mesh, double constant) : ExactSolutionScalar<double>(mesh), constant(constant) {};

  virtual double value (double x, double y) const {
    return constant; 
  };

  virtual void derivatives (double x, double y, double& dx, double& dy) const {
    dx = 0;
    dy = 0;
  };

  virtual Ord ord(Ord x, Ord y) const {
    return Ord(constant);
  }

  // Value.
  double constant;
};

class InitialSolutionEulerDensityVelX : public ExactSolutionScalar<double>
{
public:
  InitialSolutionEulerDensityVelX(Mesh* mesh, double constant) : ExactSolutionScalar<double>(mesh), constant(constant) {};

  virtual double value (double x, double y) const {
    return constant; 
  };

  virtual void derivatives (double x, double y, double& dx, double& dy) const {
    dx = 0;
    dy = 0;
  };

  virtual Ord ord(Ord x, Ord y) const {
    return Ord(constant);
  }

  // Value.
  double constant;
};

class InitialSolutionEulerDensityVelY : public ExactSolutionScalar<double>
{
public:
  InitialSolutionEulerDensityVelY(Mesh* mesh, double constant) : ExactSolutionScalar<double>(mesh), constant(constant) {};

  virtual double value (double x, double y) const {
    return constant; 
  };

  virtual void derivatives (double x, double y, double& dx, double& dy) const {
    dx = 0;
    dy = 0;
  };

  virtual Ord ord(Ord x, Ord y) const {
    return Ord(constant);
  }

  // Value.
  double constant;
};

class InitialSolutionEulerDensityEnergy : public ExactSolutionScalar<double>
{
public:
  InitialSolutionEulerDensityEnergy(Mesh* mesh, double constant) : ExactSolutionScalar<double>(mesh), constant(constant) {};

  virtual double value (double x, double y) const {
    return constant; 
  };

  virtual void derivatives (double x, double y, double& dx, double& dy) const {
    dx = 0;
    dy = 0;
  };

  virtual Ord ord(Ord x, Ord y) const {
    return Ord(constant);
  }

  // Value.
  double constant;
};

class InitialSolutionConcentration : public ExactSolutionScalar<double>
{
public:
  InitialSolutionConcentration(Mesh* mesh, double constant) : ExactSolutionScalar<double>(mesh), constant(constant) {};

  virtual double value (double x, double y) const {
    return constant; 
  };

  virtual void derivatives (double x, double y, double& dx, double& dy) const {
    dx = 0;
    dy = 0;
  };

  virtual Ord ord(Ord x, Ord y) const {
    return Ord(constant);
  }

  // Value.
  double constant;
};

/* Essential boundary condition for the coupled problem. */

class ConcentrationTimedepEssentialBC : public EssentialBoundaryCondition<double> {
public:
  ConcentrationTimedepEssentialBC(std::string marker, double constant, double startup_time) 
           : EssentialBoundaryCondition<double>(Hermes::vector<std::string>()), startup_time(startup_time), constant(constant)
  {
    markers.push_back(marker);
  }

  ~ConcentrationTimedepEssentialBC() {};

  inline EssentialBCValueType get_value_type() const { 
    return EssentialBoundaryCondition<double>::BC_FUNCTION; 
  }

  virtual double value(double x, double y, double n_x, double n_y, double t_x, double t_y) const
  {
    if(this->get_current_time() < startup_time)
      return 0.0;
    else
      if(this->get_current_time() < 2 * startup_time)
        return ((this->get_current_time() - startup_time) / startup_time) * constant;
      else
        return constant;
  }

  double startup_time;
  double constant;
};
