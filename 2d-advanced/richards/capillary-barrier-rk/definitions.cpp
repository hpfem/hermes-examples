#include "definitions.h"

// The pressure head is raised by H_OFFSET
// so that the initial condition can be taken
// as the zero vector. Note: the resulting
// pressure head will also be greater than the
// true one by this offset.
double H_OFFSET = 1000;

/* Custom weak forms */

CustomWeakFormRichardsRK::CustomWeakFormRichardsRK(ConstitutiveRelations* constitutive) : WeakForm<double>(1)
{
	// Jacobian volumetric part.
	CustomJacobianFormVol* jac_form_vol = new CustomJacobianFormVol(0, 0, constitutive);
	add_matrix_form(jac_form_vol);

	// Residual - volumetric.
	CustomResidualFormVol* res_form_vol = new CustomResidualFormVol(0, constitutive);
	add_vector_form(res_form_vol);
}

double CustomWeakFormRichardsRK::CustomJacobianFormVol::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u,
	Func<double> *v, GeomVol<double> *e, Func<double>* *ext) const
{
	double result = 0;
	Func<double>* h_prev_newton = u_ext[0];
	for (int i = 0; i < n; i++)
	{
		double h_val_i = h_prev_newton->val[i] - H_OFFSET;
		if (std::abs(h_val_i + 1000.) < 1e-7)
		{
			result += -9.8877216176350985e-05;
			continue;
		}

		double C = constitutive->C(h_val_i);
		double K = constitutive->K(h_val_i);
		double dCdh = constitutive->dCdh(h_val_i);
		double dKdh = constitutive->dKdh(h_val_i);
		double ddCdhh = constitutive->ddCdhh(h_val_i);
		double ddKdhh = constitutive->ddKdhh(h_val_i);

		double C2 = C * C;
		double a1_1 = (dKdh * C - dCdh * K) / C2;
		double a1_2 = K / C;

		double a2_1 = ((dKdh * dCdh + K * ddCdhh) * C2
			- 2 * K * C * dCdh * dCdh) / (C2 * C2);
		double a2_2 = 2 * K * dCdh / C2;

		double a3_1 = (ddKdhh * C - dKdh * dCdh) / C2;
		double a3_2 = dKdh / C;

		result += wt[i] * (-a1_1 * u->val[i] * (h_prev_newton->dx[i] * v->dx[i] + h_prev_newton->dy[i] * v->dy[i])
			- a1_2 * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i])
			+ a2_1 * u->val[i] * v->val[i] * (h_prev_newton->dx[i] * h_prev_newton->dx[i]
				+ h_prev_newton->dy[i] * h_prev_newton->dy[i])
			+ a2_2 * v->val[i] * (u->dx[i] * h_prev_newton->dx[i]
				+ u->dy[i] * h_prev_newton->dy[i])
			+ a3_1 * u->val[i] * v->val[i] * h_prev_newton->dy[i]
			+ a3_2 * v->val[i] * u->dy[i]
			);
	}
	return result;
}

Ord CustomWeakFormRichardsRK::CustomJacobianFormVol::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u,
	Func<Ord> *v, GeomVol<Ord> *e, Func<Ord>* *ext) const
{
	return Ord(10);
}

MatrixFormVol<double>* CustomWeakFormRichardsRK::CustomJacobianFormVol::clone() const
{
	CustomJacobianFormVol* toReturn = new CustomJacobianFormVol(*this);
	return toReturn;
}

double CustomWeakFormRichardsRK::CustomResidualFormVol::value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, GeomVol<double> *e,
	Func<double>* *ext) const
{
	double result = 0;
	Func<double>* h_prev_newton = u_ext[0];
	for (int i = 0; i < n; i++)
	{
		double h_val_i = h_prev_newton->val[i] - H_OFFSET;

		if (std::abs(h_val_i + 1000.) < 1e-7)
			continue;

		double C = constitutive->C(h_val_i);
		double K = constitutive->K(h_val_i);
		double dCdh = constitutive->dCdh(h_val_i);
		double dKdh = constitutive->dKdh(h_val_i);

		double r1 = (K / C);
		double r2 = K * dCdh / (C * C);
		double r3 = dKdh / C;

		result += wt[i] * (-r1 * (h_prev_newton->dx[i] * v->dx[i] + h_prev_newton->dy[i] * v->dy[i])
			+ r2 * v->val[i] * (h_prev_newton->dx[i] * h_prev_newton->dx[i]
				+ h_prev_newton->dy[i] * h_prev_newton->dy[i])
			+ r3 * v->val[i] * h_prev_newton->dy[i]
			);
	}
	return result;
}

Ord CustomWeakFormRichardsRK::CustomResidualFormVol::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, GeomVol<Ord> *e, Func<Ord>* *ext) const
{
	return Ord(10);
}

VectorFormVol<double>* CustomWeakFormRichardsRK::CustomResidualFormVol::clone() const
{
	return new CustomResidualFormVol(*this);
}