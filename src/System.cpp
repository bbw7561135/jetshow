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
//  std::cout << "calculating k_I at point " << point << std::endl;
  dxdt = jet->getKI(point, ray_direction, nu);
//  std::cout << "Step at t = " << t << std::endl;
//	double kI = jet->getKI(point, ray_direction, nu);
//  std::cout << "At t = " << t << " k_I = " << kI << std::endl;
}


void I::operator()(const double &x, double &dxdt, const double t) {
  Vector3d point = point_in + t * ray_direction;
//	std::cout << "t inside I =" << t << std::endl;
  //std::cout << "calculating transfer at point " << point << std::endl;
//	double eta = jet->getEtaI(point, ray_direction, nu);
//	double k = jet->getKI(point, ray_direction, nu);
//	double sf = eta/k;
//	std::cout << "k inside I = " << k << std::endl;
//	std::cout << "eta inside I = " << eta << std::endl;
//	std::cout << "SF = " << sf << std::endl;
  dxdt = jet->getEtaI(point, ray_direction, nu) -
      jet->getKI(point, ray_direction, nu) * x;
	// FIXME: Do i need this?
	if (dxdt < 0.0) {
		dxdt = 0.0;
	}
}


void write_cout(const double &x, const double t) {
  std::cout << "write_cout " << t << '\t' << x << std::endl;
}

bool check_opt_depth(double tau_max, const double &x) {
//	bool doStop = x >= tau_max;
//	if (doStop) {
//		std::cout << "Stopping because of tau becoming > tau_max" << std::endl;
//		std::cout << "Now tau is " << x << std::endl;
//	}
//  return doStop;
	return x >= tau_max;
}
