#ifndef EULER_UTIL_H
#define EULER_UTIL_H

#include "hermes2d.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;

// Class calculating various quantities
class QuantityCalculator
{
public:
  // Calculates energy from other quantities.
  static double calc_energy(double rho, double rho_v_x, double rho_v_y, double pressure, double kappa);
 
  // Calculates pressure from other quantities.
  static double calc_pressure(double rho, double rho_v_x, double rho_v_y, double energy, double kappa);
 
  // Calculates speed of sound.
  static double calc_sound_speed(double rho, double rho_v_x, double rho_v_y, double energy, double kappa);
};

class CFLCalculation
{
public:
  CFLCalculation(double CFL_number, double kappa);

  // If the time step is necessary to decrease / possible to increase, the value time_step will be rewritten.
  void calculate(Hermes::vector<Solution<double>*> solutions, Mesh* mesh, double & time_step) const;
  void calculate_semi_implicit(Hermes::vector<Solution<double>*> solutions, Mesh* mesh, double & time_step) const;

  void set_number(double new_CFL_number);
  
protected:
  double CFL_number;
  double kappa;
};

class ADEStabilityCalculation
{
public:
  ADEStabilityCalculation(double AdvectionRelativeConstant, double DiffusionRelativeConstant, double epsilon);

  // The method in fact returns half teh length of the shortest edge.
  double approximate_inscribed_circle_radius(Element * e);

  // If the time step is necessary to decrease / possible to increase, the value time_step will be rewritten.
  void calculate(Hermes::vector<Solution<double>*> solutions, Mesh* mesh, double & time_step);
  
protected:
  double AdvectionRelativeConstant;
  double DiffusionRelativeConstant;
  double epsilon;
};

class DiscontinuityDetector
{
public:
  /// Constructor.
  DiscontinuityDetector(Hermes::vector<const Space<double> *> spaces, 
                        Hermes::vector<Solution<double> *> solutions);

  /// Destructor.
   ~DiscontinuityDetector();

  /// Return a reference to the inner structures.
  virtual std::set<int>& get_discontinuous_element_ids() = 0;

protected:
  /// Members.
  Hermes::vector<const Space<double> *> spaces;
  Hermes::vector<Solution<double> *> solutions;
  std::set<int> discontinuous_element_ids;
  Mesh* mesh;
};

class KrivodonovaDiscontinuityDetector : public DiscontinuityDetector
{
public:
  /// Constructor.
  KrivodonovaDiscontinuityDetector(Hermes::vector<const Space<double> *> spaces, 
                        Hermes::vector<Solution<double> *> solutions);

  /// Destructor.
   ~KrivodonovaDiscontinuityDetector();

  /// Return a reference to the inner structures.
  std::set<int>& get_discontinuous_element_ids();
  std::set<int>& get_discontinuous_element_ids(double threshold);

protected:
  /// Calculates relative (w.r.t. the boundary edge_i of the Element e).
  double calculate_relative_flow_direction(Element* e, int edge_i);

  /// Calculates jumps of all solution components across the edge edge_i of the Element e.
  void calculate_jumps(Element* e, int edge_i, double result[4]);

  /// Calculates h.
  double calculate_h(Element* e, int polynomial_order);

  /// Calculates the norm of the solution on the central element.
  void calculate_norms(Element* e, int edge_i, double result[4]);
};

class KuzminDiscontinuityDetector : public DiscontinuityDetector
{
public:
  /// Constructor.
  KuzminDiscontinuityDetector(Hermes::vector<const Space<double> *> spaces, 
                        Hermes::vector<Solution<double> *> solutions, bool limit_all_orders_independently = false);

  /// Destructor.
   ~KuzminDiscontinuityDetector();

  /// Return a reference to the inner structures.
  std::set<int>& get_discontinuous_element_ids();

  /// Return a reference to the inner structures.
  std::set<int>& get_second_order_discontinuous_element_ids();

  /// Returns info about the method.
  bool get_limit_all_orders_independently();
protected:
  /// Center.
  void find_centroid_values(Hermes::Hermes2D::Element* e, double u_c[4]);
  void find_centroid_derivatives(Hermes::Hermes2D::Element* e, double u_dx_c[4], double u_dy_c[4]);
  void find_second_centroid_derivatives(Hermes::Hermes2D::Element* e, double u_dxx_c[4], double u_dxy_c[4], double u_dyy_c[4]);

  /// Vertices.
  void find_vertex_values(Hermes::Hermes2D::Element* e, double vertex_values[4][4]);
  void find_vertex_derivatives(Hermes::Hermes2D::Element* e, double vertex_derivatives[4][4][2]);

  /// Logic - 1st order.
  void find_u_i_min_max_first_order(Hermes::Hermes2D::Element* e, double u_i_min[4][4], double u_i_max[4][4]);
  void find_alpha_i_first_order(double u_i_min[4][4], double u_i_max[4][4], double u_c[4], double u_i[4][4], double alpha_i[4]);
  void find_alpha_i_first_order_real(Hermes::Hermes2D::Element* e, double u_i[4][4], double u_c[4], double u_dx_c[4], double u_dy_c[4], double alpha_i_real[4]);

  /// Logic - 2nd order.
  void find_u_i_min_max_second_order(Hermes::Hermes2D::Element* e, double u_d_i_min[4][4][2], double u_d_i_max[4][4][2]);
  void find_alpha_i_second_order(double u_d_i_min[4][4][2], double u_d_i_max[4][4][2], double u_dx_c[4], double u_dy_c[4], double u_d_i[4][4][2], double alpha_i[4]);
  void find_alpha_i_second_order_real(Hermes::Hermes2D::Element* e, double u_i[4][4][2], double u_dx_c[4], double u_dy_c[4], double u_dxx_c[4], double u_dxy_c[4], double u_dyy_c[4], double alpha_i_real[4]);

private:
  /// For limiting of second order terms.
  std::set<int> second_order_discontinuous_element_ids;
  bool limit_all_orders_independently;
};

class FluxLimiter
{
public:
  /// Enumeration of types.
  /// Used to pick the proper DiscontinuityDetector.
  enum LimitingType
  {
    Krivodonova,
    Kuzmin
  };
  /// Constructor.
  FluxLimiter(LimitingType type, double* solution_vector, Hermes::vector<const Space<double> *> spaces, bool Kuzmin_limit_all_orders_independently = false);

  /// Destructor.
   ~FluxLimiter();

  /// Do the limiting.
  /// With the possibility to also limit the spaces from which the spaces in the constructors are refined.
  virtual void limit_according_to_detector(Hermes::vector<Space<double> *> coarse_spaces_to_limit = Hermes::vector<Space<double> *>());
  
  /// For Kuzmin's detector.
  virtual void limit_second_orders_according_to_detector(Hermes::vector<Space<double> *> coarse_spaces_to_limit = Hermes::vector<Space<double> *>());
  
  void get_limited_solutions(Hermes::vector<Solution<double>*> solutions_to_limit);
protected:
  /// Members.
  double* solution_vector;
  Hermes::vector<const Space<double> *> spaces;
  DiscontinuityDetector* detector;
  Hermes::vector<Solution<double>*> limited_solutions;
};

// Filters.
class MachNumberFilter : public Hermes::Hermes2D::SimpleFilter<double>
{
public: 
  MachNumberFilter(Hermes::vector<MeshFunction<double>*> solutions, double kappa) : SimpleFilter<double>(solutions), kappa(kappa) {};
  ~MachNumberFilter() {};
protected:
  virtual void filter_fn(int n, Hermes::vector<double*> values, double* result);

  double kappa;
};

class PressureFilter : public Hermes::Hermes2D::SimpleFilter<double>
{
public: 
  PressureFilter(Hermes::vector<MeshFunction<double>*> solutions, double kappa) : SimpleFilter<double>(solutions), kappa(kappa) {};
  ~PressureFilter() {};
protected:
  virtual void filter_fn(int n, Hermes::vector<double*> values, double* result);

  double kappa;
};

class EntropyFilter : public Hermes::Hermes2D::SimpleFilter<double>
{
public: 
  EntropyFilter(Hermes::vector<MeshFunction<double>*> solutions, double kappa, double rho_ext, double p_ext) : SimpleFilter<double>(solutions), kappa(kappa), rho_ext(rho_ext), p_ext(p_ext) {};
  ~EntropyFilter() {};
protected:
  virtual void filter_fn(int n, Hermes::vector<double*> values, double* result);

  double kappa, rho_ext, p_ext;
};

class EulerFluxes
{
public:
  EulerFluxes(double kappa) : kappa(kappa) {}

  template<typename Scalar>
  Scalar A_1_0_0(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy) {

    return Scalar(0.0);
  }

  template<typename Scalar>
  Scalar A_1_0_1(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy) {
    return Scalar(1.0);
  }

  template<typename Scalar>
  Scalar A_1_0_2(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy) {
    return Scalar(0.0);
  }

  template<typename Scalar>
  Scalar A_1_0_3(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy) {
    return Scalar(0.0);
  }

  template<typename Scalar>
  Scalar A_2_0_0(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy) {
    return Scalar(0.0);
  }

  template<typename Scalar>
  Scalar A_2_0_1(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy) {
    return Scalar(0.0);
  }

  template<typename Scalar>
  Scalar A_2_0_2(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy) {
    return Scalar(1.0);
  }

  template<typename Scalar>
  Scalar A_2_0_3(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy) {
    return Scalar(0.0);
  }

  template<typename Scalar>
  Scalar A_1_1_0(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy) {
    return Scalar(- ((rho_v_x * rho_v_x) / (rho * rho)) + 0.5 * (kappa - 1.0) * 
            ((rho_v_x * rho_v_x + rho_v_y * rho_v_y) /   (rho * rho)));
  }

  template<typename Scalar>
  Scalar A_1_1_1(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy){
    return Scalar((3. - kappa) * (rho_v_x / rho));
  }

  template<typename Scalar>
  Scalar A_1_1_2(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy){
    return Scalar((1.0 - kappa) * (rho_v_y / rho));
  }

  template<typename Scalar>
  Scalar A_1_1_3(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy){
    return Scalar(kappa - 1.);
  }

  template<typename Scalar>
  Scalar A_2_1_0(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy){
    return Scalar(- rho_v_x * rho_v_y / (rho * rho));
  }

  template<typename Scalar>
  Scalar A_2_1_1(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy){
    return Scalar(rho_v_y / rho);
  }

  template<typename Scalar>
  Scalar A_2_1_2(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy){
    return Scalar(rho_v_x / rho);
  }

  template<typename Scalar>
  Scalar A_2_1_3(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy){
    return Scalar(0);
  }

  template<typename Scalar>
  Scalar A_1_2_0(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy){
    return Scalar(- rho_v_x * rho_v_y / (rho * rho));
  }

  template<typename Scalar>
  Scalar A_1_2_1(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy){
    return Scalar(rho_v_y / rho);
  }

  template<typename Scalar>
  Scalar A_1_2_2(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy){
    return Scalar(rho_v_x / rho);
  }

  template<typename Scalar>
  Scalar A_1_2_3(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy){
    return Scalar(0);
  }

  template<typename Scalar>
  Scalar A_2_2_0(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy){
    return Scalar(- ((rho_v_y * rho_v_y) / (rho * rho)) + 0.5 * (kappa - 1.0) 
            * ((rho_v_x * rho_v_x + rho_v_y * rho_v_y) /   (rho * rho)));
  }

  template<typename Scalar>
  Scalar A_2_2_1(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy){
    return Scalar((1.0 - kappa) * (rho_v_x / rho));
  }

  template<typename Scalar>
  Scalar A_2_2_2(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy){
    return Scalar((3.0 - kappa) * (rho_v_y / rho));
  }

  template<typename Scalar>
  Scalar A_2_2_3(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy){
    return Scalar(kappa - 1.);
  }

  template<typename Scalar>
  Scalar A_1_3_0(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy){
    return Scalar((rho_v_x / rho) * (((kappa - 1.0) * ((rho_v_x * rho_v_x + rho_v_y * rho_v_y) / (rho * rho)))
      - (kappa * energy / rho)));
  }

  template<typename Scalar>
  Scalar A_1_3_1(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy){
    return Scalar((kappa * energy / rho) - (kappa - 1.0) * rho_v_x * rho_v_x / (rho * rho)
      - 0.5 * (kappa - 1.0) * (rho_v_x * rho_v_x + rho_v_y * rho_v_y) / (rho * rho));
  }

  template<typename Scalar>
  Scalar A_1_3_2(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy){
    return Scalar((1.0 - kappa) * (rho_v_x * rho_v_y) / (rho * rho));
  }

  template<typename Scalar>
  Scalar A_1_3_3(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy){
    return Scalar(kappa * (rho_v_x / rho));
  }

  template<typename Scalar>
  Scalar A_2_3_0(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy){
    return Scalar(- (rho_v_y * energy) / (rho * rho) - (rho_v_y / (rho * rho)) * (kappa - 1.0) 
            * (energy - ((rho_v_x * rho_v_x + rho_v_y * rho_v_y) / (2 * rho))) + (rho_v_y / rho) 
            * (kappa - 1.0) * ((rho_v_x * rho_v_x + rho_v_y * rho_v_y) / (2 * rho * rho)));
  }

  template<typename Scalar>
  Scalar A_2_3_1(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy){
    return Scalar((1.0 - kappa) * (rho_v_x * rho_v_y) / (rho * rho));
  }

  template<typename Scalar>
  Scalar A_2_3_2(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy){
    return Scalar((energy / rho) + (1 / rho) * (kappa - 1.0) * ( energy - ((rho_v_x * rho_v_x 
            + rho_v_y * rho_v_y) / (2 * rho))) + (1.0 - kappa) * ((rho_v_y * rho_v_y) / (rho * rho)));
  }

  template<typename Scalar>
  Scalar A_2_3_3(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy){
    return Scalar(kappa * rho_v_y / rho);
  }
  protected:
    double kappa;
};

#endif
