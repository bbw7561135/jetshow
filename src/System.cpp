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
  Vector3d point = point_in + t * ray_direction / pc;
//  std::cout << "calculating k_I at point " << point << std::endl;
  dxdt = jet->getKI(point, ray_direction, nu);
//  std::cout << "Step at t = " << t << std::endl;
}



void I::operator()(const double &x, double &dxdt, const double t) {
  Vector3d point = point_in + t * ray_direction / pc;
	std::cout << "t inside I =" << t << std::endl;
  //std::cout << "calculating transfer at point " << point << std::endl;
	double eta = jet->getEtaI(point, ray_direction, nu);
	std::cout << "eta inside I = " << eta << std::endl;
  dxdt = jet->getEtaI(point, ray_direction, nu) +
      jet->getKI(point, ray_direction, nu) * x;
}


void write_cout(const double &x, const double t) {
  std::cout << "write_cout " << t << '\t' << x << std::endl;
}

bool check_opt_depth(double tau_max, const double &x) {
  return x >= tau_max;
}
