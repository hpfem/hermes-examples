#include "hermes2d.h"


/* Namespaces used */

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;
using namespace Hermes::Hermes2D::RefinementSelectors;

/* Essential boundary conditions */

class EssentialBCNonConst : public EssentialBoundaryCondition<double> 
{
public:
    EssentialBCNonConst(std::string marker);

    ~EssentialBCNonConst() {};

    virtual EssentialBoundaryCondition<double>::EssentialBCValueType get_value_type() const;

    virtual double value(double x, double y, double n_x, double n_y, double t_x, double t_y) const;
};

/* Weak forms */

class WeakFormHelmholtz : public WeakForm<double>
{
public:
    WeakFormHelmholtz(double eps, double mu, double omega, double sigma, double beta, double E0, double h);

private:
    class MatrixFormHelmholtzEquation_real_real : public MatrixFormVol<double>
    {
    public:
        MatrixFormHelmholtzEquation_real_real(int i, int j, double eps, double omega, double mu)
              : MatrixFormVol<double>(i, j), eps(eps), omega(omega), mu(mu) {};

        template<typename Real, typename Scalar>
        Scalar matrix_form(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, 
            Geom<Real> *e, Func<Scalar>* *ext) const;

        virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
            Func<double> *v, Geom<double> *e, Func<double>* *ext) const;

        virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, 
            Func<Ord> *v, Geom<Ord> *e, Func<Ord>* *ext) const;

        MatrixFormVol<double>* clone() const;
    private:
        // Members.
        double eps;
        double omega;
        double mu;
    };

    class MatrixFormHelmholtzEquation_real_imag : public MatrixFormVol<double>
    {
    private:
        // Members.
        double mu;
        double omega;
        double sigma;

    public:
        MatrixFormHelmholtzEquation_real_imag(int i, int j, double mu, double omega, double sigma) 
              : MatrixFormVol<double>(i, j), mu(mu), omega(omega), sigma(sigma) {};

        template<typename Real, typename Scalar>
        Scalar matrix_form(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, 
            Geom<Real> *e, Func<Scalar>* *ext) const;

        virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
            Func<double> *v, Geom<double> *e, Func<double>* *ext) const;

        virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
            Geom<Ord> *e, Func<Ord>* *ext) const;
        MatrixFormVol<double>* clone() const;
    };

    class MatrixFormHelmholtzEquation_imag_real : public MatrixFormVol<double>
    {
    private:
        // Members.
        double mu;
        double omega;
        double sigma;

    public:
        MatrixFormHelmholtzEquation_imag_real(int i, int j, double mu, double omega, double sigma) 
              : MatrixFormVol<double>(i, j), mu(mu), omega(omega), sigma(sigma) {};

        template<typename Real, typename Scalar>
        Scalar matrix_form(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, 
            Geom<Real> *e, Func<Scalar>* *ext) const;

        virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
            Func<double> *v, Geom<double> *e, Func<double>* *ext) const;

        virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
            Geom<Ord> *e, Func<Ord>* *ext) const;
        MatrixFormVol<double>* clone() const;
    };

    class MatrixFormHelmholtzEquation_imag_imag : public MatrixFormVol<double>
    {
    public:
        MatrixFormHelmholtzEquation_imag_imag(int i, int j, double eps, double mu, double omega) 
              : MatrixFormVol<double>(i, j), eps(eps), mu(mu), omega(omega) {};

        template<typename Real, typename Scalar>
        Scalar matrix_form(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, 
            Geom<Real> *e, Func<Scalar>* *ext) const;

        virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
            Func<double> *v, Geom<double> *e, Func<double>* *ext) const;

        virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
            Geom<Ord> *e, Func<Ord>* *ext) const;
        
        MatrixFormVol<double>* clone() const;

        // Members.
        double eps;
        double mu;
        double omega;
    };

    class MatrixFormSurfHelmholtz_real_imag : public MatrixFormSurf<double>
    {
    private:
        double beta;
    public:
        MatrixFormSurfHelmholtz_real_imag(int i, int j, std::string area, double beta)
              : MatrixFormSurf<double>(i, j), beta(beta){ this->set_area(area); };

        template<typename Real, typename Scalar>
        Scalar matrix_form_surf(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, 
            Geom<Real> *e, Func<Scalar>* *ext) const;

        virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
            Func<double> *v, Geom<double> *e, Func<double>* *ext) const;

        virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
            Geom<Ord> *e, Func<Ord>* *ext) const;
        MatrixFormSurf<double>* clone() const;
    };

    class MatrixFormSurfHelmholtz_imag_real : public MatrixFormSurf<double>
    {

    private:
        double beta;
    public:
        MatrixFormSurfHelmholtz_imag_real(int i, int j, std::string area, double beta)
              : MatrixFormSurf<double>(i, j), beta(beta){ this->set_area(area); };

        template<typename Real, typename Scalar>
        Scalar matrix_form_surf(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, 
            Geom<Real> *e, Func<Scalar>* *ext) const;

        virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, 
            Geom<double> *e, Func<double>* *ext) const;

        virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
            Geom<Ord> *e, Func<Ord>* *ext) const;

        MatrixFormSurf<double>* clone() const;
    };

    class VectorFormHelmholtzEquation_real : public VectorFormVol<double>
    {
    public:
    VectorFormHelmholtzEquation_real(int i, double eps, double omega, double mu, double sigma)
      : VectorFormVol<double>(i), eps(eps), omega(omega), mu(mu), sigma(sigma) {};

        template<typename Real, typename Scalar>
        Scalar vector_form(int n, double *wt, Func<Real> *u_ext[], Func<Real> *v, 
            Geom<Real> *e, Func<Scalar>* *ext) const;

        virtual double value(int n, double *wt, Func<double> *u_ext[],
            Func<double> *v, Geom<double> *e, Func<double>* *ext) const;

        virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[],
            Func<Ord> *v, Geom<Ord> *e, Func<Ord>* *ext) const;

        VectorFormVol<double>* clone() const;
    private:
        // Members.
        double eps;
        double omega;
        double mu;
        double sigma;
    };

    class VectorFormHelmholtzEquation_imag : public VectorFormVol<double>
    {
    public:
    VectorFormHelmholtzEquation_imag(int i, double eps, double omega, double mu, double sigma)
      : VectorFormVol<double>(i), eps(eps), omega(omega), mu(mu), sigma(sigma) {};

        template<typename Real, typename Scalar>
        Scalar vector_form(int n, double *wt, Func<Real> *u_ext[], Func<Real> *v, 
            Geom<Real> *e, Func<Scalar>* *ext) const;

        virtual double value(int n, double *wt, Func<double> *u_ext[],
            Func<double> *v, Geom<double> *e, Func<double>* *ext) const;

        virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[],
            Func<Ord> *v, Geom<Ord> *e, Func<Ord>* *ext) const;

        VectorFormVol<double>* clone() const;
    private:
        // Members.
        double eps;
        double omega;
        double mu;
        double sigma;
    };

    class VectorFormSurfHelmholtz_real : public VectorFormSurf<double>
    {
    private:
        double beta;
    public:
        VectorFormSurfHelmholtz_real(int i, std::string area, double beta)
              : VectorFormSurf<double>(i), beta(beta) { this->set_area(area); };

        template<typename Real, typename Scalar>
        Scalar vector_form_surf(int n, double *wt, Func<Real> *u_ext[], Func<Real> *v, 
            Geom<Real> *e, Func<Scalar>* *ext) const;

        virtual double value(int n, double *wt, Func<double> *u_ext[],
            Func<double> *v, Geom<double> *e, Func<double>* *ext) const;

        virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, 
            Geom<Ord> *e, Func<Ord>* *ext) const;
        
        VectorFormSurf<double>* clone() const;
    };

    class VectorFormSurfHelmholtz_imag : public VectorFormSurf<double>
    {
    private:
        double beta;
    public:
        VectorFormSurfHelmholtz_imag(int i, std::string area, double beta)
              : VectorFormSurf<double>(i), beta(beta) { this->set_area(area); };

        template<typename Real, typename Scalar>
        Scalar vector_form_surf(int n, double *wt, Func<Real> *u_ext[], Func<Real> *v, 
            Geom<Real> *e, Func<Scalar>* *ext) const;

        virtual double value(int n, double *wt, Func<double> *u_ext[],
            Func<double> *v, Geom<double> *e, Func<double>* *ext) const;

        virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, 
            Geom<Ord> *e, Func<Ord>* *ext) const;
        
        VectorFormSurf<double>* clone() const;
    };
};
