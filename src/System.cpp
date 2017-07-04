#include "System.h"


System::System(Jet *newjet,
               Vector3d &newpoint_in,
               Vector3d &newray_direction,
               double newnu) {
  jet = newjet;
  point_in = newpoint_in;
  ray_direction = newray_direction;
  nu = newnu;
}


void Tau::operator()(const double &x, double &dxdt, const double t) {
  Vector3d point = point_in + t * ray_direction;
  dxdt = jet->getKI(point, ray_direction, nu);
}


void I::operator()(const double &x, double &dxdt, const double t) {
  Vector3d point = point_in + t * ray_direction;
  dxdt = jet->getEtaI(point, ray_direction, nu) -
      jet->getKI(point, ray_direction, nu) * x;
	// FIXME: Do i need this?
	double x0 = x;
	if (x0 < 0) {
		x0 = 0.0;
	}
}


void write_cout(const double &x, const double t) {
  std::cout << "write_cout " << t << '\t' << x << std::endl;
}


bool check_opt_depth(double tau_max, const double &x) {
	return x >= tau_max;
}


bool check_opt_depth_xt(double tau_max, std::pair<double, double> xt) {
	auto x = xt.first;
	return x >= tau_max;
}