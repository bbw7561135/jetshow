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


//FullStokes::FullStokes(Jet *newjet, Vector3d &newpoint_in,
//                       Vector3d &newray_direction, double newnu) {
//	jet = newjet;
//	point_in = newpoint_in;
//	ray_direction = newray_direction;
//	nu = newnu;
//}
//
//void FullStokes::operator()(const state_type &x, state_type &dxdt,
//                            const double t) {
//	Vector3d point = point_in + t * ray_direction;
//	dxdt[0] = jet->getEtaI(point, ray_direction, nu) -
//			jet->getKI(point, ray_direction, nu) * x[0] -
//			jet->getKQ(point, ray_direction, nu) * x[1] -
//			jet->getKU(point, ray_direction, nu) * x[2] -
//			jet->getKV(point, ray_direction, nu) * x[3];
//	dxdt[1] = jet->getEtaQ(point, ray_direction, nu) -
//	          jet->getKI(point, ray_direction, nu) * x[1] -
//	          jet->getKQ(point, ray_direction, nu) * x[0] -
//	          jet->getKF(point, ray_direction, nu) * x[2] -
//	          jet->gethQ(point, ray_direction, nu) * x[3];
//	dxdt[2] = jet->getEtaU(point, ray_direction, nu) -
//	          jet->getKI(point, ray_direction, nu) * x[2] -
//	          jet->getKU(point, ray_direction, nu) * x[0] +
//	          jet->getKF(point, ray_direction, nu) * x[1] -
//	          jet->getKC(point, ray_direction, nu) * x[3];
//	dxdt[3] = jet->getEtaV(point, ray_direction, nu) -
//	          jet->getKI(point, ray_direction, nu) * x[3] -
//	          jet->getKV(point, ray_direction, nu) * x[0] +
//	          jet->gethQ(point, ray_direction, nu) * x[1] +
//	          jet->getKC(point, ray_direction, nu) * x[2];
//
//	state_type x0 = x;
//	if (x0[0] < 0) {
//		x0[0] = 0.0;
//	}
//}


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

