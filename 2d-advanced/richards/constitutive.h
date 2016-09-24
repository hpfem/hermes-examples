#define HERMES_REPORT_ALL

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
	virtual double ddCdhh(double h) = 0;
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

	// ddC/dhh (Gardner).
	double ddCdhh(double h)
	{
		if (h < 0) return alpha * alpha * (theta_s - theta_r) * alpha * exp(alpha * h);
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
			k_s*std::pow((1 + std::pow((-alpha*h), n)), (-m / 2))*std::pow((1 -
				std::pow((-alpha*h), (m*n))*std::pow((1 + std::pow((-alpha*h), n)), (-m))), 2);
		else return k_s;
	}

	// dK/dh (van Genuchten).
	double dKdh(double h)
	{
		if (h < 0) return
			k_s*std::pow((1 + std::pow((-alpha*h), n)), (-m / 2))*(1 -
				std::pow((-alpha*h), (m*n))*std::pow((1 +
					std::pow((-alpha*h), n)), (-m)))*(-2 * m*n*std::pow((-alpha*h), (m*n))*std::pow((1 +
						std::pow((-alpha*h), n)), (-m)) / h +
						2 * m*n*std::pow((-alpha*h), n)*std::pow((-alpha*h), (m*n))*std::pow((1 +
							std::pow((-alpha*h), n)), (-m)) / (h*(1 + std::pow((-alpha*h), n)))) -
			k_s*m*n*std::pow((-alpha*h), n)*std::pow((1 + std::pow((-alpha*h), n)), (-m / 2))*std::pow((1 -
				std::pow((-alpha*h), (m*n))*std::pow((1 + std::pow((-alpha*h), n)), (-m))), 2) / (2 * h*(1 +
					std::pow((-alpha*h), n)));
		else return 0;
	}

	// ddK/dhh (van Genuchten).
	double ddKdhh(double h)
	{
		if (h < 0) return
			k_s*std::pow((1 + std::pow((-alpha*h), n)), (-m / 2))*(1 -
				std::pow((-alpha*h), (m*n))*std::pow((1 +
					std::pow((-alpha*h), n)), (-m)))*(-2 * std::pow(m, 2)*std::pow(n, 2)*std::pow((-alpha*h), (m*n))*std::pow((1
						+ std::pow((-alpha*h), n)), (-m)) / std::pow(h, 2) +
						2 * m*n*std::pow((-alpha*h), (m*n))*std::pow((1 + std::pow((-alpha*h), n)), (-m)) / std::pow(h, 2)
						- 2 * m*n*std::pow((-alpha*h), n)*std::pow((-alpha*h), (m*n))*std::pow((1 +
							std::pow((-alpha*h), n)), (-m)) / (std::pow(h, 2)*(1 + std::pow((-alpha*h), n))) -
						2 * m*std::pow(n, 2)*std::pow((-alpha*h), (2 * n))*std::pow((-alpha*h), (m*n))*std::pow((1 +
							std::pow((-alpha*h), n)), (-m)) / (std::pow(h, 2)*std::pow((1 + std::pow((-alpha*h), n)), 2)) -
						2 * std::pow(m, 2)*std::pow(n, 2)*std::pow((-alpha*h), (2 * n))*std::pow((-alpha*h), (m*n))*std::pow((1
							+ std::pow((-alpha*h), n)), (-m)) / (std::pow(h, 2)*std::pow((1 + std::pow((-alpha*h), n)), 2)) +
						2 * m*std::pow(n, 2)*std::pow((-alpha*h), n)*std::pow((-alpha*h), (m*n))*std::pow((1 +
							std::pow((-alpha*h), n)), (-m)) / (std::pow(h, 2)*(1 + std::pow((-alpha*h), n))) +
						4 * std::pow(m, 2)*std::pow(n, 2)*std::pow((-alpha*h), n)*std::pow((-alpha*h), (m*n))*std::pow((1 +
							std::pow((-alpha*h), n)), (-m)) / (std::pow(h, 2)*(1 + std::pow((-alpha*h), n)))) +
			k_s*std::pow((1 + std::pow((-alpha*h), n)), (-m / 2))*(-m*n*std::pow((-alpha*h), (m*n))*std::pow((1
				+ std::pow((-alpha*h), n)), (-m)) / h +
				m*n*std::pow((-alpha*h), n)*std::pow((-alpha*h), (m*n))*std::pow((1 +
					std::pow((-alpha*h), n)), (-m)) / (h*(1 +
						std::pow((-alpha*h), n))))*(-2 * m*n*std::pow((-alpha*h), (m*n))*std::pow((1 +
							std::pow((-alpha*h), n)), (-m)) / h +
							2 * m*n*std::pow((-alpha*h), n)*std::pow((-alpha*h), (m*n))*std::pow((1 +
								std::pow((-alpha*h), n)), (-m)) / (h*(1 + std::pow((-alpha*h), n)))) -
			k_s*m*n*std::pow((-alpha*h), n)*std::pow((1 + std::pow((-alpha*h), n)), (-m / 2))*(1 -
				std::pow((-alpha*h), (m*n))*std::pow((1 +
					std::pow((-alpha*h), n)), (-m)))*(-2 * m*n*std::pow((-alpha*h), (m*n))*std::pow((1 +
						std::pow((-alpha*h), n)), (-m)) / h +
						2 * m*n*std::pow((-alpha*h), n)*std::pow((-alpha*h), (m*n))*std::pow((1 +
							std::pow((-alpha*h), n)), (-m)) / (h*(1 + std::pow((-alpha*h), n)))) / (h*(1 +
								std::pow((-alpha*h), n))) + k_s*m*n*std::pow((-alpha*h), n)*std::pow((1 +
									std::pow((-alpha*h), n)), (-m / 2))*std::pow((1 - std::pow((-alpha*h), (m*n))*std::pow((1 +
										std::pow((-alpha*h), n)), (-m))), 2) / (2 * std::pow(h, 2)*(1 + std::pow((-alpha*h), n))) +
			k_s*m*std::pow(n, 2)*std::pow((-alpha*h), (2 * n))*std::pow((1 +
				std::pow((-alpha*h), n)), (-m / 2))*std::pow((1 - std::pow((-alpha*h), (m*n))*std::pow((1 +
					std::pow((-alpha*h), n)), (-m))), 2) / (2 * std::pow(h, 2)*std::pow((1 +
						std::pow((-alpha*h), n)), 2)) - k_s*m*std::pow(n, 2)*std::pow((-alpha*h), n)*std::pow((1 +
							std::pow((-alpha*h), n)), (-m / 2))*std::pow((1 - std::pow((-alpha*h), (m*n))*std::pow((1 +
								std::pow((-alpha*h), n)), (-m))), 2) / (2 * std::pow(h, 2)*(1 + std::pow((-alpha*h), n))) +
			k_s*std::pow(m, 2)*std::pow(n, 2)*std::pow((-alpha*h), (2 * n))*std::pow((1 +
				std::pow((-alpha*h), n)), (-m / 2))*std::pow((1 - std::pow((-alpha*h), (m*n))*std::pow((1 +
					std::pow((-alpha*h), n)), (-m))), 2) / (4 * std::pow(h, 2)*std::pow((1 +
						std::pow((-alpha*h), n)), 2));

		else return 0;
	}

	// C (van Genuchten).
	double C(double h)
	{
		if (h < 0) return
			storativity*std::pow((1 + std::pow((-alpha*h), n)), (-m))*(theta_s - theta_r) / theta_s -
			m*n*std::pow((-alpha*h), n)*std::pow((1 + std::pow((-alpha*h), n)), (-m))*(theta_s - theta_r) / (h*(1 + std::pow((-alpha*h), n)));
		else return storativity;
	}

	// dC/dh (van Genuchten).
	double dCdh(double h)
	{
		if (h < 0) return
			m*n*std::pow((-alpha*h), n)*std::pow((1 + std::pow((-alpha*h), n)), (-m))*(theta_s - theta_r) / (std::pow(h, 2)*(1 + std::pow((-alpha*h), n))) + m*std::pow(n, 2)*std::pow((-alpha*h), (2 * n))*std::pow((1 +
				std::pow((-alpha*h), n)), (-m))*(theta_s - theta_r) / (std::pow(h, 2)*std::pow((1 + std::pow((-alpha*h), n)), 2)) + std::pow(m, 2)*std::pow(n, 2)*std::pow((-alpha*h), (2 * n))*std::pow((1 + std::pow((-alpha*h), n)), (-m))*(theta_s - theta_r) / (std::pow(h, 2)*
					std::pow((1 + std::pow((-alpha*h), n)), 2)) - m*std::pow(n, 2)*std::pow((-alpha*h), n)*std::pow((1 + std::pow((-alpha*h), n)), (-m))*(theta_s - theta_r) / (std::pow(h, 2)*(1 + std::pow((-alpha*h), n))) -
			m*n*storativity*std::pow((-alpha*h), n)*std::pow((1 + std::pow((-alpha*h), n)), (-m))*(theta_s - theta_r) / (theta_s*h*(1 + std::pow((-alpha*h), n)));
		else return 0;
	}

	// ddC/dhh (Gardner).
	/// \todo This is not correct, the correct relation should be here.
	// it is e.g. used in basic-rk-newton.
	double ddCdhh(double h)
	{
		return 0.0;
	}
protected:
	double m, n, storativity;
};

class ConstitutiveRelationsGenuchtenWithLayer : public ConstitutiveRelationsGenuchten
{
public:
	ConstitutiveRelationsGenuchtenWithLayer(int constitutive_table_method, int num_inside_pts, double low_limit, double table_precision,
		double table_limit, double* k_s_vals, double* alpha_vals, double* n_vals, double* m_vals, double* theta_r_vals, double* theta_s_vals, double* storativity_vals) : ConstitutiveRelationsGenuchten(0, 0, 0, 0, 0, 0, 0), constitutive_table_method(constitutive_table_method),
		num_inside_pts(num_inside_pts), polynomials_ready(false), low_limit(low_limit), table_precision(table_precision), table_limit(table_limit), k_s_vals(k_s_vals), alpha_vals(alpha_vals), n_vals(n_vals), m_vals(m_vals), theta_r_vals(theta_r_vals), theta_s_vals(theta_s_vals), storativity_vals(storativity_vals), constitutive_tables_ready(false),
		polynomials_allocated(false), k_table(NULL), dKdh_table(NULL), ddKdhh_table(NULL), c_table(NULL), dCdh_table(NULL), polynomials(NULL), pol_search_help(NULL), k_pols(NULL), c_pols(NULL)
	{}
	ConstitutiveRelationsGenuchtenWithLayer(double alpha, double m, double n, double theta_s, double theta_r, double k_s, double storativity, int constitutive_table_method, int num_inside_pts, double low_limit, double table_precision,
		double table_limit, double* k_s_vals, double* alpha_vals, double* n_vals, double* m_vals, double* theta_r_vals, double* theta_s_vals, double* storativity_vals)
		: ConstitutiveRelationsGenuchten(alpha, m, n, theta_s, theta_r, k_s, storativity), constitutive_table_method(constitutive_table_method),
		num_inside_pts(num_inside_pts), polynomials_ready(false), low_limit(low_limit), table_precision(table_precision), table_limit(table_limit), k_s_vals(k_s_vals), alpha_vals(alpha_vals), n_vals(n_vals), m_vals(m_vals), theta_r_vals(theta_r_vals), theta_s_vals(theta_s_vals), storativity_vals(storativity_vals), constitutive_tables_ready(false),
		polynomials_allocated(false), k_table(NULL), dKdh_table(NULL), ddKdhh_table(NULL), c_table(NULL), dCdh_table(NULL), polynomials(NULL), pol_search_help(NULL), k_pols(NULL), c_pols(NULL)
	{}

	// This function uses Horner scheme to efficiently evaluate values of 
	// polynomials in approximating K(h) function close to full saturation.
	double horner(double *pol, double x, int n) {
		double px = 0.0;
		for (int i = 0; i < n; i++) {
			px = px*x + pol[n - 1 - i];
		}
		return px;
	}

	// K (van Genuchten).
	double K(double h) {
		return K(h, 0);
	}
	double K(double h, int layer)
	{
		double value;
		int location;
		if (h > low_limit && h < 0 && polynomials_ready && constitutive_table_method == 1) {
			return horner(polynomials[layer][0], h, (6 + num_inside_pts));
		}
		else
		{
			if (constitutive_table_method == 0 || !constitutive_tables_ready || h < table_limit) {
				alpha = alpha_vals[layer];
				n = n_vals[layer];
				m = m_vals[layer];
				k_s = k_s_vals[layer];
				theta_r = theta_r_vals[layer];
				theta_s = theta_s_vals[layer];
				storativity = storativity_vals[layer];
				if (h < 0) return
					k_s*(std::pow(1 - std::pow(-(alpha*h), m*n) /
						std::pow(1 + std::pow(-(alpha*h), n), m), 2) /
						std::pow(1 + std::pow(-(alpha*h), n), m / 2.));
				else return k_s;
			}
			else if (h < 0 && constitutive_table_method == 1) {
				location = -int(h / table_precision);
				value = (k_table[layer][location + 1] - k_table[layer][location])*(-h / table_precision - location) + k_table[layer][location];
				return value;
			}
			else if (h < 0 && constitutive_table_method == 2) {
				location = pol_search_help[int(-h)];
				return horner(k_pols[location][layer][0], h, 6);
			}
			else return k_s_vals[layer];
		}
	}

	// dK/dh (van Genuchten).
	double dKdh(double h) {
		return dKdh(h, 0);
	}
	double dKdh(double h, int layer)
	{
		double value;
		int location;

		if (h > low_limit && h < 0 && polynomials_ready && constitutive_table_method == 1) {
			return horner(polynomials[layer][1], h, (5 + num_inside_pts));
		}
		else
		{
			if (constitutive_table_method == 0 || !constitutive_tables_ready || h < table_limit) {
				alpha = alpha_vals[layer];
				n = n_vals[layer];
				m = m_vals[layer];
				k_s = k_s_vals[layer];
				theta_r = theta_r_vals[layer];
				theta_s = theta_s_vals[layer];
				storativity = storativity_vals[layer];
				if (h < 0) return
					k_s*((alpha*std::pow(-(alpha*h), -1 + n)*
						std::pow(1 + std::pow(-(alpha*h), n), -1 - m / 2.)*
						std::pow(1 - std::pow(-(alpha*h), m*n) /
							std::pow(1 + std::pow(-(alpha*h), n), m), 2)*m*n) / 2. +
							(2 * (1 - std::pow(-(alpha*h), m*n) /
								std::pow(1 + std::pow(-(alpha*h), n), m))*
								(-(alpha*std::pow(-(alpha*h), -1 + n + m*n)*
									std::pow(1 + std::pow(-(alpha*h), n), -1 - m)*m*n) +
									(alpha*std::pow(-(alpha*h), -1 + m*n)*m*n) /
									std::pow(1 + std::pow(-(alpha*h), n), m))) /
						std::pow(1 + std::pow(-(alpha*h), n), m / 2.));
				else return 0;
			}
			else if (h < 0 && constitutive_table_method == 1) {
				location = -int(h / table_precision);
				value = (dKdh_table[layer][location + 1] - dKdh_table[layer][location])*(-h / table_precision - location) + dKdh_table[layer][location];
				return value;
			}
			else if (h < 0 && constitutive_table_method == 2) {
				location = pol_search_help[int(-h)];
				return horner(k_pols[location][layer][1], h, 5);
			}
			else return 0;
		}
	}

	// ddK/dhh (van Genuchten).
	double ddKdhh(double h) {
		return ddKdhh(h, 0);
	}
	double ddKdhh(double h, int layer)
	{

		int location;
		double value;


		if (h > low_limit && h < 0 && polynomials_ready && constitutive_table_method == 1) {
			return horner(polynomials[layer][2], h, (4 + num_inside_pts));
		}
		else
		{
			if (constitutive_table_method == 0 || !constitutive_tables_ready || h < table_limit) {
				alpha = alpha_vals[layer];
				n = n_vals[layer];
				m = m_vals[layer];
				k_s = k_s_vals[layer];
				theta_r = theta_r_vals[layer];
				theta_s = theta_s_vals[layer];
				storativity = storativity_vals[layer];

				if (h < 0) return
					k_s*(-(std::pow(alpha, 2)*std::pow(-(alpha*h), -2 + n)*
						std::pow(1 + std::pow(-(alpha*h), n), -1 - m / 2.)*
						std::pow(1 - std::pow(-(alpha*h), m*n) /
							std::pow(1 + std::pow(-(alpha*h), n), m), 2)*m*(-1 + n)*n) / 2.\
						- (std::pow(alpha, 2)*std::pow(-(alpha*h), -2 + 2 * n)*
							std::pow(1 + std::pow(-(alpha*h), n), -2 - m / 2.)*
							std::pow(1 - std::pow(-(alpha*h), m*n) /
								std::pow(1 + std::pow(-(alpha*h), n), m), 2)*(-1 - m / 2.)*m*
							std::pow(n, 2)) / 2. + 2 * alpha*std::pow(-(alpha*h), -1 + n)*
						std::pow(1 + std::pow(-(alpha*h), n), -1 - m / 2.)*
						(1 - std::pow(-(alpha*h), m*n) /
							std::pow(1 + std::pow(-(alpha*h), n), m))*m*n*
							(-(alpha*std::pow(-(alpha*h), -1 + n + m*n)*
								std::pow(1 + std::pow(-(alpha*h), n), -1 - m)*m*n) +
								(alpha*std::pow(-(alpha*h), -1 + m*n)*m*n) /
								std::pow(1 + std::pow(-(alpha*h), n), m)) +
								(2 * std::pow(-(alpha*std::pow(-(alpha*h), -1 + n + m*n)*
									std::pow(1 + std::pow(-(alpha*h), n), -1 - m)*m*n) +
									(alpha*std::pow(-(alpha*h), -1 + m*n)*m*n) /
									std::pow(1 + std::pow(-(alpha*h), n), m), 2)) /
						std::pow(1 + std::pow(-(alpha*h), n), m / 2.) +
						(2 * (1 - std::pow(-(alpha*h), m*n) /
							std::pow(1 + std::pow(-(alpha*h), n), m))*
							(std::pow(alpha, 2)*std::pow(-(alpha*h), -2 + 2 * n + m*n)*
								std::pow(1 + std::pow(-(alpha*h), n), -2 - m)*(-1 - m)*m*
								std::pow(n, 2) + std::pow(alpha, 2)*
								std::pow(-(alpha*h), -2 + n + m*n)*
								std::pow(1 + std::pow(-(alpha*h), n), -1 - m)*std::pow(m, 2)*
								std::pow(n, 2) - (std::pow(alpha, 2)*
									std::pow(-(alpha*h), -2 + m*n)*m*n*(-1 + m*n)) /
								std::pow(1 + std::pow(-(alpha*h), n), m) +
								std::pow(alpha, 2)*std::pow(-(alpha*h), -2 + n + m*n)*
								std::pow(1 + std::pow(-(alpha*h), n), -1 - m)*m*n*
								(-1 + n + m*n))) / std::pow(1 + std::pow(-(alpha*h), n), m / 2.));

				else return 0;
			}
			else if (h < 0 && constitutive_table_method == 1) {
				location = -int(h / table_precision);
				value = (ddKdhh_table[layer][location + 1] -
					ddKdhh_table[layer][location])*(-h / table_precision - location) + ddKdhh_table[layer][location];
				return value;
			}
			else if (h < 0 && constitutive_table_method == 2) {
				location = pol_search_help[int(-h)];
				return horner(k_pols[location][layer][2], h, 4);
			}
			else return 0;
		}
	}

	// C (van Genuchten).
	double C(double h) {
		return C(h, 0);
	}
	double C(double h, int layer)
	{
		int location;
		double value;

		if (constitutive_table_method == 0 || !constitutive_tables_ready || h < table_limit) {
			alpha = alpha_vals[layer];
			n = n_vals[layer];
			m = m_vals[layer];
			k_s = k_s_vals[layer];
			theta_r = theta_r_vals[layer];
			theta_s = theta_s_vals[layer];
			storativity = storativity_vals[layer];
			if (h < 0) return
				alpha*std::pow(-(alpha*h), -1 + n)*
				std::pow(1 + std::pow(-(alpha*h), n), -1 - m)*m*n*(theta_s - theta_r) +
				(storativity*((theta_s - theta_r) / std::pow(1 + std::pow(-(alpha*h), n), m) +
					theta_r)) / theta_s;
			else return storativity;
		}
		else if (h < 0 && constitutive_table_method == 1) {
			location = -int(h / table_precision);
			value = (c_table[layer][location + 1] - c_table[layer][location])*(-h / table_precision - location) + c_table[layer][location];
			return value;
		}
		else if (h < 0 && constitutive_table_method == 2) {
			location = pol_search_help[int(-h)];
			return horner(c_pols[location][layer][0], h, 4);
		}
		else return storativity_vals[layer];
	}

	// dC/dh (van Genuchten).
	double dCdh(double h) {
		return dCdh(h, 0);
	}
	double dCdh(double h, int layer)
	{
		int location;
		double value;

		if (constitutive_table_method == 0 || !constitutive_tables_ready || h < table_limit) {
			alpha = alpha_vals[layer];
			n = n_vals[layer];
			m = m_vals[layer];
			k_s = k_s_vals[layer];
			theta_r = theta_r_vals[layer];
			theta_s = theta_s_vals[layer];
			storativity = storativity_vals[layer];
			if (h < 0) return
				-(std::pow(alpha, 2)*std::pow(-(alpha*h), -2 + n)*
					std::pow(1 + std::pow(-(alpha*h), n), -1 - m)*m*(-1 + n)*n*
					(theta_s - theta_r)) - std::pow(alpha, 2)*std::pow(-(alpha*h), -2 + 2 * n)*
				std::pow(1 + std::pow(-(alpha*h), n), -2 - m)*(-1 - m)*m*std::pow(n, 2)*
				(theta_s - theta_r) + (alpha*std::pow(-(alpha*h), -1 + n)*
					std::pow(1 + std::pow(-(alpha*h), n), -1 - m)*m*n*storativity*
					(theta_s - theta_r)) / theta_s;
			else return 0;
		}
		else if (h < 0 && constitutive_table_method == 1) {
			location = -int(h / table_precision);
			value = (dCdh_table[layer][location + 1] -
				dCdh_table[layer][location])*(-h / table_precision - location) + dCdh_table[layer][location];
			return value;
		}
		else if (h < 0 && constitutive_table_method == 2) {
			location = pol_search_help[int(-h)];
			return horner(c_pols[location][layer][1], h, 3);
		}
		else return 0;

	}

	int constitutive_table_method, num_inside_pts;
	bool polynomials_ready;
	double low_limit;
	double table_precision;
	double table_limit;

	double* k_s_vals;
	double* alpha_vals;
	double* n_vals;
	double* m_vals;
	double* theta_r_vals;
	double* theta_s_vals;
	double* storativity_vals;
	bool constitutive_tables_ready;
	bool polynomials_allocated;

	double** k_table;
	double** dKdh_table;
	double** ddKdhh_table;
	double** c_table;
	double** dCdh_table;
	double*** polynomials;
	int* pol_search_help;
	double**** k_pols;
	double**** c_pols;
};


bool init_polynomials(int n, double low_limit, double *points, int n_inside_points, int layer, ConstitutiveRelationsGenuchtenWithLayer* constitutive, int material_count, int num_of_intervals, double* intervals_4_approx);

bool get_constitutive_tables(int method, ConstitutiveRelationsGenuchtenWithLayer* constitutive, int material_count);