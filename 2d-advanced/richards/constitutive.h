class ConstitutiveRelations
{
public:
  ConstitutiveRelations(double alpha, double theta_s, double theta_r, double k_s) : alpha(alpha), theta_s(theta_s), theta_r(theta_r), k_s(k_s)
  {}

  virtual double K(double h) = 0;
  virtual double dKdh(double h) = 0;
  virtual double ddKdhh(double h) = 0;
  virtual double C(double h) = 0;
  virtual double dCdh(double h) = 0;

protected:
  double alpha, theta_s, theta_r, k_s;
};

class ConstitutiveRelationsGardner : public ConstitutiveRelations
{
public:
  ConstitutiveRelationsGardner(double alpha, double theta_s, double theta_r, double k_s) : ConstitutiveRelations(alpha, theta_s, theta_r, k_s)
  {}

  // K (Gardner).
  double K(double h)
  {
    if (h < 0) return k_s*exp(alpha*h);
    else return k_s;    
  }

  // dK/dh (Gardner).
  double dKdh(double h)
  {
    if (h < 0) return k_s*alpha*exp(alpha*h);
    else return 0;
  }

  // ddK/dhh (Gardner).
  double ddKdhh(double h)
  {
    if (h < 0) return k_s*alpha*alpha*exp(alpha*h);
    else return 0;
  }

  // C (Gardner).
  double C(double h)
  {
    if (h < 0) return alpha*(theta_s - theta_r)*exp(alpha*h);
    else return alpha*(theta_s - theta_r);    
  }

  // dC/dh (Gardner).
  double dCdh(double h)
  {
    if (h < 0) return alpha*(theta_s - theta_r)*alpha*exp(alpha*h);
    else return 0;    
  }

};

class ConstitutiveRelationsGenuchten : public ConstitutiveRelations
{
public:
  ConstitutiveRelationsGenuchten(double alpha, double m, double n, double theta_s, double theta_r, double k_s, double storativity) : ConstitutiveRelations(alpha, theta_s, theta_r, k_s), m(m), n(n), storativity(storativity)
  {}

  // K (van Genuchten).
  double K(double h)
  {
    if (h < 0) return 
      k_s*std::pow((1 + std::pow((-alpha*h),n)),(-m/2))*std::pow((1 -
      std::pow((-alpha*h),(m*n))*std::pow((1 + std::pow((-alpha*h),n)),(-m))),2) ;
    else return k_s;    
  }

  // dK/dh (van Genuchten).
  double dKdh(double h)
  {
    if (h < 0) return 
      k_s*std::pow((1 + std::pow((-alpha*h),n)),(-m/2))*(1 -
      std::pow((-alpha*h),(m*n))*std::pow((1 +
      std::pow((-alpha*h),n)),(-m)))*(-2*m*n*std::pow((-alpha*h),(m*n))*std::pow((1 +
      std::pow((-alpha*h),n)),(-m))/h +
      2*m*n*std::pow((-alpha*h),n)*std::pow((-alpha*h),(m*n))*std::pow((1 +
      std::pow((-alpha*h),n)),(-m))/(h*(1 + std::pow((-alpha*h),n)))) -
      k_s*m*n*std::pow((-alpha*h),n)*std::pow((1 + std::pow((-alpha*h),n)),(-m/2))*std::pow((1 -
      std::pow((-alpha*h),(m*n))*std::pow((1 + std::pow((-alpha*h),n)),(-m))),2)/(2*h*(1 +
      std::pow((-alpha*h),n))) ;
    else return 0;
  }

  // ddK/dhh (van Genuchten).
  double ddKdhh(double h)
  {
    if (h < 0) return 
      k_s*std::pow((1 + std::pow((-alpha*h),n)),(-m/2))*(1 -
      std::pow((-alpha*h),(m*n))*std::pow((1 +
      std::pow((-alpha*h),n)),(-m)))*(-2*std::pow(m,2)*std::pow(n,2)*std::pow((-alpha*h),(m*n))*std::pow((1
      + std::pow((-alpha*h),n)),(-m))/std::pow(h,2) +
      2*m*n*std::pow((-alpha*h),(m*n))*std::pow((1 + std::pow((-alpha*h),n)),(-m))/std::pow(h,2)
      - 2*m*n*std::pow((-alpha*h),n)*std::pow((-alpha*h),(m*n))*std::pow((1 +
      std::pow((-alpha*h),n)),(-m))/(std::pow(h,2)*(1 + std::pow((-alpha*h),n))) -
      2*m*std::pow(n,2)*std::pow((-alpha*h),(2*n))*std::pow((-alpha*h),(m*n))*std::pow((1 +
      std::pow((-alpha*h),n)),(-m))/(std::pow(h,2)*std::pow((1 + std::pow((-alpha*h),n)),2)) -
      2*std::pow(m,2)*std::pow(n,2)*std::pow((-alpha*h),(2*n))*std::pow((-alpha*h),(m*n))*std::pow((1
      + std::pow((-alpha*h),n)),(-m))/(std::pow(h,2)*std::pow((1 + std::pow((-alpha*h),n)),2)) +
      2*m*std::pow(n,2)*std::pow((-alpha*h),n)*std::pow((-alpha*h),(m*n))*std::pow((1 +
      std::pow((-alpha*h),n)),(-m))/(std::pow(h,2)*(1 + std::pow((-alpha*h),n))) +
      4*std::pow(m,2)*std::pow(n,2)*std::pow((-alpha*h),n)*std::pow((-alpha*h),(m*n))*std::pow((1 +
      std::pow((-alpha*h),n)),(-m))/(std::pow(h,2)*(1 + std::pow((-alpha*h),n)))) +
      k_s*std::pow((1 + std::pow((-alpha*h),n)),(-m/2))*(-m*n*std::pow((-alpha*h),(m*n))*std::pow((1
      + std::pow((-alpha*h),n)),(-m))/h +
      m*n*std::pow((-alpha*h),n)*std::pow((-alpha*h),(m*n))*std::pow((1 +
      std::pow((-alpha*h),n)),(-m))/(h*(1 +
      std::pow((-alpha*h),n))))*(-2*m*n*std::pow((-alpha*h),(m*n))*std::pow((1 +
      std::pow((-alpha*h),n)),(-m))/h +
      2*m*n*std::pow((-alpha*h),n)*std::pow((-alpha*h),(m*n))*std::pow((1 +
      std::pow((-alpha*h),n)),(-m))/(h*(1 + std::pow((-alpha*h),n)))) -
      k_s*m*n*std::pow((-alpha*h),n)*std::pow((1 + std::pow((-alpha*h),n)),(-m/2))*(1 -
      std::pow((-alpha*h),(m*n))*std::pow((1 +
      std::pow((-alpha*h),n)),(-m)))*(-2*m*n*std::pow((-alpha*h),(m*n))*std::pow((1 +
      std::pow((-alpha*h),n)),(-m))/h +
      2*m*n*std::pow((-alpha*h),n)*std::pow((-alpha*h),(m*n))*std::pow((1 +
      std::pow((-alpha*h),n)),(-m))/(h*(1 + std::pow((-alpha*h),n))))/(h*(1 +
      std::pow((-alpha*h),n))) + k_s*m*n*std::pow((-alpha*h),n)*std::pow((1 +
      std::pow((-alpha*h),n)),(-m/2))*std::pow((1 - std::pow((-alpha*h),(m*n))*std::pow((1 +
      std::pow((-alpha*h),n)),(-m))),2)/(2*std::pow(h,2)*(1 + std::pow((-alpha*h),n))) +
      k_s*m*std::pow(n,2)*std::pow((-alpha*h),(2*n))*std::pow((1 +
      std::pow((-alpha*h),n)),(-m/2))*std::pow((1 - std::pow((-alpha*h),(m*n))*std::pow((1 +
      std::pow((-alpha*h),n)),(-m))),2)/(2*std::pow(h,2)*std::pow((1 +
      std::pow((-alpha*h),n)),2)) - k_s*m*std::pow(n,2)*std::pow((-alpha*h),n)*std::pow((1 +
      std::pow((-alpha*h),n)),(-m/2))*std::pow((1 - std::pow((-alpha*h),(m*n))*std::pow((1 +
      std::pow((-alpha*h),n)),(-m))),2)/(2*std::pow(h,2)*(1 + std::pow((-alpha*h),n))) +
      k_s*std::pow(m,2)*std::pow(n,2)*std::pow((-alpha*h),(2*n))*std::pow((1 +
      std::pow((-alpha*h),n)),(-m/2))*std::pow((1 - std::pow((-alpha*h),(m*n))*std::pow((1 +
      std::pow((-alpha*h),n)),(-m))),2)/(4*std::pow(h,2)*std::pow((1 +
      std::pow((-alpha*h),n)),2)) ;

    else return 0;
  }

  // C (van Genuchten).
  double C(double h)
  {
    if (h < 0) return 
      storativity*std::pow((1 + std::pow((-alpha*h),n)),(-m))*(theta_s - theta_r)/theta_s - 
      m*n*std::pow((-alpha*h),n)*std::pow((1 + std::pow((-alpha*h),n)),(-m))*(theta_s - theta_r)/(h*(1 + std::pow((-alpha*h),n)));
    else return storativity;    
  }

  // dC/dh (van Genuchten).
  double dCdh(double h)
  {
    if (h < 0) return
      m*n*std::pow((-alpha*h),n)*std::pow((1 + std::pow((-alpha*h),n)),(-m))*(theta_s - theta_r)/(std::pow(h,2)*(1 + std::pow((-alpha*h),n))) + m*std::pow(n,2)*std::pow((-alpha*h),(2*n))*std::pow((1 + 
      std::pow((-alpha*h),n)),(-m))*(theta_s - theta_r)/(std::pow(h,2)*std::pow((1 + std::pow((-alpha*h),n)),2)) + std::pow(m,2)*std::pow(n,2)*std::pow((-alpha*h),(2*n))*std::pow((1 + std::pow((-alpha*h),n)),(-m))*(theta_s - theta_r)/(std::pow(h,2)*
      std::pow((1 + std::pow((-alpha*h),n)),2)) - m*std::pow(n,2)*std::pow((-alpha*h),n)*std::pow((1 + std::pow((-alpha*h),n)),(-m))*(theta_s - theta_r)/(std::pow(h,2)*(1 + std::pow((-alpha*h),n))) - 
      m*n*storativity*std::pow((-alpha*h),n)*std::pow((1 + std::pow((-alpha*h),n)),(-m))*(theta_s - theta_r)/(theta_s*h*(1 + std::pow((-alpha*h),n))) ;
    else return 0;    
  }

private:
  double m, n, storativity;
};