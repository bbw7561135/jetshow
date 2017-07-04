#include <math.h>
#include <string>
#include <iostream>
#include <fstream>
#include "utils.h"



//double eta_0(double n, double b) {
//    return  pi*nu_p(n)*nu_p(n)*nu_b(b)*m_e/c;
//}

//double k_0(double nu, double n, double b) {
//    return pi*nu_p(n)*nu_p(n)*nu_b(b)/(c*nu*nu);
//}

//double eta_i(double nu, double n, double b, double sin_theta,
//             double s) {
//    return eta_0(n, b)*sin_theta*pow(nu_b(b)*sin_theta/nu, (s-1.)/2.)*pow(3., s/2.)/(2.*(s+1.))*tgamma(s/4.+19./12.)*tgamma(s/4.-1./12.);
//}
double nu_p(double n) {return sqrt(n*q_e*q_e / (pi*m_e));}

double nu_b(Vector3d &b, Vector3d &n_los) {
    return q_e*(n_los.cross(b)).norm()/(2.*pi*m_e*c);
}

double nu_b_value(Vector3d &b) {
	return q_e*b.norm()/(2.*pi*m_e*c);
}

double k_0(Vector3d &b, Vector3d &n_los, double nu, double n) {
    return pi*nu_p(n)*nu_p(n)*nu_b(b, n_los)/(c*nu*nu);
}

double k_0_value(Vector3d &b, double nu, double n) {
	return pi*nu_p(n)*nu_p(n)*nu_b_value(b)/(c*nu*nu);
}

double k_i(Vector3d &b, Vector3d &n_los, double nu, double n, double s) {
    double factor = (pow(3., (s+1.)/2.)/4.)*tgamma(s/4.+11./6.)*tgamma(s/4.+1./6.);
    return k_0(b, n_los, nu, n) * pow(nu_b(b, n_los)/nu, s/2.) * factor;
}

double k_q(Vector3d &b, Vector3d &n_los, double nu, double n, double s) {
	return (s+2.)/(s+10./3)*k_i(b, n_los, nu, n, s);
}

double k_u(Vector3d &b, Vector3d &n_los, double nu, double n, double s) {
	return 0;
}

double k_v(Vector3d &b, Vector3d &n_los, double nu, double n, double s) {
	double cos_theta = b.dot(n_los)/b.norm();
	double factor = pow(3., s/2.)*(s+3.)*(s+2.)/(4.*(s+1.))*tgamma(s/4.+11./12.)*tgamma(s/4.+7./12.);
	return -k_0_value(b, nu, n)*cos_theta*pow(nu_b(b, n_los)/nu, (s+1.)/2.) * factor;
}

double k_F_c(Vector3d &b, Vector3d &n_los, double nu, double n, double s) {
	double cos_theta = b.dot(n_los)/b.norm();
	return 2.*pi*nu_p(n)*nu_p(n)*nu_b_value(b)*cos_theta/(c*nu*nu);
}

double k_C_c(Vector3d &b, Vector3d &n_los, double nu, double n, double s) {
	return -pi*nu_p(n)*nu_p(n)*nu_b(b, n_los)*nu_b(b, n_los)/(c*nu*nu*nu);
}

double
k_F_r(Vector3d &b, Vector3d &n_los, double nu, double n, double gamma_min,
      double s) {
	return (s+2.)*log(gamma_min)/((s+1.)*pow(gamma_min, s+1.))*k_F_c(b, n_los, nu,
	                                                                 n, s);
}

double
k_C_r(Vector3d &b, Vector3d &n_los, double nu, double n, double gamma_min,
      double s) {
	return (2./(s-2.))*(pow(gamma_min, 2.-s) -
			pow(nu_b(b, n_los)/nu, (s-2.)/2.))*k_C_c(b, n_los, nu, n, s);
}


double h_Q(Vector3d &b, Vector3d &n_los, double nu, double n, double s) {
	return 0;
}



double eta_0(Vector3d &b, Vector3d &n_los, double n) {
    return pi*nu_p(n)*nu_p(n)*nu_b(b, n_los)*m_e/c;
}

double eta_0_value(Vector3d &b, double n) {
	return pi*nu_p(n)*nu_p(n)*nu_b_value(b)*m_e/c;
}

double eta_i(Vector3d &b, Vector3d &n_los, double nu, double n, double s) {
    double factor = pow(3., s/2.)/(2.*(s+1))*tgamma(s/4.+19./12.)*tgamma(s/4.-1./12.);
    return eta_0(b, n_los, n) * pow(nu_b(b, n_los)/nu, (s-1.)/2.) * factor;
}

double eta_q(Vector3d &b, Vector3d &n_los, double nu, double n, double s) {
	return (s+1.0)/(s+7./3.)*eta_i(b, n_los, nu, n, s);
}

double eta_u(Vector3d &b, Vector3d &n_los, double nu, double n, double s) {
	return 0;
}

double eta_v(Vector3d &b, Vector3d &n_los, double nu, double n, double s) {
	double cos_theta = b.dot(n_los)/b.norm();
	double factor = pow(3., (s-1.)/2.)*(s+2.)/(2.*s)*tgamma(s/4.+2./3.)*tgamma(s/4.+1./3.);
	return -eta_0_value(b, n)*cos_theta*pow(nu_b(b, n_los)/nu, s/2.) * factor;
}


double getG(Vector3d &v) {
    Vector3d beta = v/c;
    return sqrt(1./(1.- beta.squaredNorm()));
}

double getD(Vector3d &n_los, Vector3d &v) {
    Vector3d beta = v/c;
    return 1./(getG(v)*(1.-beta.dot(n_los)));
}

Vector3d get_n_los_prime(Vector3d &n_los, Vector3d &v) {
    double df = getD(n_los, v);
    double gamma = getG(v);
    Vector3d beta = v/c;
    return df*n_los - (df+1.)*(gamma/(gamma+1.))*beta;
}

Vector3d getBprime(Vector3d &b, Vector3d &v) {
    double gamma = getG(v);
    Vector3d result = b/gamma + gamma/((1.+gamma)*(c*c))*v*v.dot(b);
    return result;
}


double comoving_transfer_distance(double z, double H0, double omega_M,
																	double omega_V) {
	Ctd ctd(z);
	double ctd_state = 0.0;
	using namespace boost::numeric::odeint;
	integrate<double>(ctd, ctd_state, 0., z, 0.00001);
	double result = (100./H0) * 3. * pow(10., 9.) * ctd_state;
	return result;
}

double pc_to_mas(double z) {
	double d_a = comoving_transfer_distance(z)/(1.+z);
	double angle_rads = 1./d_a;
	return rad_to_mas*angle_rads;
}

double mas_to_pc(double z) {
	double d_a = comoving_transfer_distance(z)/(1.+z);
	return mas_to_rad*d_a;
}


std::ostream &
write_2dvector(std::ostream &os, std::vector<std::vector<double>> &v,
							 double scale) {
	for (int i = 0; i < v.size(); ++i)
	{
		for (int j = 0; j < v[i].size(); ++j)
		{
			double value = v[i][j]/scale;
			os << value <<" ";
		}
		os<<"\n";
	}
	return os;
}

double pixel_solid_angle(double pixel_size_mas, double z) {
	double pixel_size_pc = pixel_size_mas * mas_to_pc(z);
	double d_a = comoving_transfer_distance(z)/(1.+z);
	return pixel_size_pc*pixel_size_pc/(d_a*d_a);

}


// Generates ``n`` random directions (normalized vectors).
std::vector<Vector3d> generate_random_directions(int n, unsigned int seed) {
	std::vector<Vector3d> points;
  boost::mt19937 gen;
	gen.seed(seed);
  boost::uniform_on_sphere<double> dist(3);

	for (int i = 0; i < n; ++i) {
		std::vector<double> res = dist(gen);
		Vector3d v = Vector3d(res.data());
		points.push_back(v);
	}
	return points;
}


// Generates ``n`` random points inside spherical volume with restrictions on
// radius and polar angle. Density profile ``exponent``.
// FIXME: Implement ``theta_min`` > 0 case.
std::vector<Vector3d> generate_random_points_sphere(int n, double r_max,
                                                    double exponent,
                                                    unsigned int seed,
                                                    double r_min,
                                                    double theta_max,
                                                    double theta_min)
{
	std::vector<Vector3d> points;
  boost::mt19937 gen;
	gen.seed(seed);
	boost::random::uniform_real_distribution<double> boost_distrib(0.0, 1.0);

	// Default unrestricted polar angle
	double costheta_lim = 1.0;
	// If polar angle is restricted
	if (theta_max > 0) {
		costheta_lim = cos(theta_max);
	}

	double r;
	double phi;
	double u;
	double theta;

	for (int j = 0; j < n; ++j) {
		r = r_min + (r_max - r_min) * boost_distrib(gen);
		phi = 2.0 * pi * boost_distrib(gen);
		u = (1.0 - costheta_lim) * (2.0 * boost_distrib(gen) - 1.0);
		if (u < 0.0) {
			u = u - costheta_lim;
		} else {
			u = u + costheta_lim;
		}
		theta = acos(u);

		// dV = d(r3)d(cos\theta)d(\phi) => for uniform density we must take a cube
		// root of ``r``.
		// FIXME: What with r^(-3) density profiles?
		r = pow(r, 1./(3.0-exponent));
		Vector3d v = Vector3d{r*sin(theta)*cos(phi),
		                      r*sin(theta)*sin(phi),
		                      r*cos(theta)};
		points.push_back(v);
	}

	return points;
}

std::vector<Vector3d>
generate_random_points_general(int n, double r_max, Geometry *geo,
                               double exponent, unsigned int seed) {
	std::vector<Vector3d> points;
	boost::mt19937 gen;
	gen.seed(seed);
	boost::random::uniform_real_distribution<double> boost_distrib(0.0, 1.0);


	double r;
	double phi;
	double u;
	double costheta;
	double theta;
	int j = 0;
	Vector3d v;

	while (j < n) {
		r = r_max*boost_distrib(gen);
		phi = 2.0 * pi * boost_distrib(gen);
		costheta = 2.0 * boost_distrib(gen) - 1.0;
		theta = acos(costheta);

		// dV = d(r3)d(cos\theta)d(\phi) => for uniform density we must take a cube
		// root of ``r``.
		// FIXME: What with r^(-3) density profiles?
		r = pow(r, 1./(3.0-exponent));
		v = Vector3d{r*sin(theta)*cos(phi),
		             r*sin(theta)*sin(phi),
		             r*cos(theta)};
		if (geo->is_within(v)) {
			points.push_back(v);
			++j;
		}
	}

	return points;
}

int steps_schedule(double tau, int n_min, int n_max, double tau_min,
                      double tau_max) {
	if (tau <= tau_min) {
		return n_min;
	} else if (tau >= tau_max) {
		return n_max;
	} else {
	double b = -log(tau_min);
	double k = n_max/pow((log(tau_max) + b), 2.0);
	double n = k * pow((log(tau) + b), 2.0) + n_min;
	return ceil(n);
	}
}


Ctd::Ctd(double z, double H0, double omega_M, double omega_V) : z(z), H0(H0),
																																omega_M(omega_M),
																																omega_V(omega_V) {}

void Ctd::operator()(const double &x, double &dxdt, const double t) {
	dxdt = pow(omega_M*(1.+t*t*t)+omega_V, -0.5);
};