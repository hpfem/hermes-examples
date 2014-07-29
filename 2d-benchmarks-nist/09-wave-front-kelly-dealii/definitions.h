#include "hermes2d.h"
#include "../NIST-util.h"

using namespace Hermes::Hermes2D;
using namespace WeakFormsH1;
using Hermes::Ord;

/* Right-hand side */

class CustomWeakForm : public WeakForm<double>
{
  class Jacobian : public MatrixFormVol<double>
  {
  public:
    Jacobian() : MatrixFormVol<double>(0, 0) { this->setSymFlag(HERMES_SYM); };

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u,
                        Func<double> *v, GeomVol<double> *e, Func<double>* *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
                    GeomVol<Ord> *e, Func<Ord>* *ext) const;

    MatrixFormVol<double>* clone() const { return new Jacobian(); }
  };
  
  class Residual : public VectorFormVol<double>
  {
    const Hermes::Hermes2DFunction<double>* rhs;
  public:
    Residual(const Hermes::Hermes2DFunction<double>* rhs) : VectorFormVol<double>(0), rhs(rhs) {};

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, 
                         GeomVol<double> *e, Func<double>* *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                    GeomVol<Ord> *e, Func<Ord>* *ext) const;
    VectorFormVol<double>* clone() const { return new Residual(rhs); }
  };
  
  public:
    CustomWeakForm(const Hermes::Hermes2DFunction<double>* rhs)
    {
      add_matrix_form(new Jacobian);
      add_vector_form(new Residual(rhs));
    }
};

class CustomRightHandSide : public Hermes::Hermes2DFunction<double>
{
public:
  CustomRightHandSide(double alpha, double x_loc, double y_loc, double r_zero)
    : Hermes::Hermes2DFunction<double>(), alpha(alpha), x_loc(x_loc), y_loc(y_loc), r_zero(r_zero) 
  { };

  virtual double value(double x, double y) const;
  virtual Ord value (Ord x, Ord y) const { return Ord(8); }
  
  double alpha, x_loc, y_loc, r_zero;
};

/* Exact solution */

class CustomExactSolution : public ExactSolutionScalar<double>
{
public:
  CustomExactSolution(MeshSharedPtr mesh, double alpha, double x_loc, double y_loc, double r_zero)
    : ExactSolutionScalar<double>(mesh), alpha(alpha), x_loc(x_loc), y_loc(y_loc), r_zero(r_zero) 
  { };

  virtual double value(double x, double y) const {
    return atan(alpha * (sqrt(pow(x - x_loc, 2) + pow(y - y_loc, 2)) - r_zero));
  };

  virtual void derivatives (double x, double y, double& dx, double& dy) const;
  virtual Ord ord (double x, double y) const { return Ord(Ord::get_max_order()); }

  MeshFunction<double>* clone() const { return new CustomExactSolution(mesh, alpha, x_loc, y_loc, r_zero); }

  double alpha, x_loc, y_loc, r_zero;
};

/* Bilinear form inducing the energy norm */

class EnergyErrorForm : public Adapt<double>::MatrixFormVolError
{
public:
  EnergyErrorForm(WeakForm<double> *problem_wf) : Adapt<double>::MatrixFormVolError(0, 0, HERMES_UNSET_NORM)
  {
    this->form = problem_wf->get_mfvol()[0];
    this->wf = problem_wf;
  }

  virtual double value(int n, double *wt, Func<double> *u_ext[],
                       Func<double> *u, Func<double> *v, GeomVol<double> *e,
                       Func<double>* *ext) const
  {
    return this->form->value(n, wt, u_ext, u, v, e, ext);
  }

  virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[],
                  Func<Ord> *u, Func<Ord> *v, GeomVol<Ord> *e,
                  Func<Ord>* *ext) const
  {
    return this->form->ord(n, wt, u_ext, u, v, e, ext);
  }

  MatrixFormVol<double>* clone() const { return new EnergyErrorForm(wf); }
private:
  const MatrixFormVol<double>* form;
};

/* Linear form for the residual error estimator */

class ResidualErrorForm : public KellyTypeAdapt<double>::ErrorEstimatorForm
{
public:
  ResidualErrorForm(CustomRightHandSide* rhs) 
    : KellyTypeAdapt<double>::ErrorEstimatorForm(0), rhs(rhs) 
  { };

  double value(int n, double *wt, 
               Func<double> *u_ext[], Func<double> *u, 
               GeomVol<double> *e, Func<double>* *ext) const;

  Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, 
          GeomVol<Ord> *e, Func<Ord>* *ext) const;

private:
  CustomRightHandSide* rhs;

};

class ConvergenceTable
{
  public:  
    void save(const char* filename) const;
     
    void add_column(const std::string& name, const std::string& format) {
      columns.push_back(Column(name, format));
    }
    
    void add_value(unsigned int col, double x) {
      add_value_internal<double>(col, x, "%g");
    }
    void add_value(unsigned int col, int x) {
      add_value_internal<int>(col, x, "%d");
    }
    
    int num_rows() const { 
      return columns[0].data.size(); 
    }
    int num_columns() const { 
      return columns.size(); 
    }
    
  private:
    struct Column
    {
      Column(const std::string& name, const std::string& format) 
        : label(name), format(format) 
      { }
      
      struct Entry
      {
        union
        {
          int ivalue;
          double dvalue;
        };
        enum { INT, DOUBLE } data_type;
        
        Entry(int x) : ivalue(x), data_type(INT) {};
        Entry(double x) : dvalue(x), data_type(DOUBLE) {};
      };
      
      std::string label;
      std::string format;
      std::vector<Entry> data;
    };
    
    std::vector<Column> columns;
    
    template <typename T>
    void add_value_internal(unsigned int col, T x, const std::string& default_fmt) 
    {
      if (columns.size() == 0) add_column("", default_fmt);
      if (col < 0 || col >= columns.size()) 
        throw Hermes::Exceptions::Exception("Invalid column number.");
      
      columns[col].data.push_back(Column::Entry(x));
    }
};

// Convert int to string.
std::string itos(const unsigned int i);