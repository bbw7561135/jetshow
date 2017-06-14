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
  //std::cout << "calculating k_I at point " << point << std::endl;
  dxdt = jet->getKI(point, ray_direction, nu);
}

void write_cout(const double &x, const double t) {
  std::cout << t << '\t' << x << std::endl;
}
