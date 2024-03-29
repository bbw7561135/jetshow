#ifndef JETSHOW_UTILS_H
#define JETSHOW_UTILS_H

#include <math.h>
#include <vector>
#include <Eigen/Eigen>
#include <boost/math/constants/constants.hpp>
#include <boost/numeric/odeint.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_on_sphere.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include "Geometry.h"

using Eigen::Vector3d;


const double mas_to_rad = 4.84813681109536e-09;
const double rad_to_mas = 206264806.24709633;
// Parsec [cm]
const double pc = 3.0856775814671913e+18;
// Mass of electron [g]
const double m_e = 9.10938356e-28;
// Mass of proton [g]
const double m_p = 1.672621898e-24;
// Charge of electron [C]
 const double q_e = 4.80320467299766e-10;
//const double q_e = 1.6*1E-19;
// Charge of proton [C]
 const double q_p = q_e;
//const double q_p = 1.6*1E-19;
// Speed of light [cm / s]
const double c = 29979245800.0;
// Jy in cgc
const double to_jy = 1E+23;
// pi
const double pi = boost::math::constants::pi<double>();


double nu_p(double n);

double nu_b(Vector3d &b, Vector3d &n_los);

double sin_theta(Vector3d &b, Vector3d &n_los);

double k_0(Vector3d &b, Vector3d &n_los, double nu, double n);

// For random B-field
double k_0(double b, Vector3d &n_los, double nu, double n);

// Absorption coefficient for given vector of magnetic field ``b``, unit LOS
// vector ``n_los`` and others measured in emission frame
double k_i(Vector3d &b, Vector3d &n_los, double nu, double n, double s=2.5);

// For random B-field
double k_i(double b, Vector3d &n_los, double nu, double n, double s=2.5);

// For random B-field - alternative formulation
double k_i_(double b, Vector3d &n_los, double nu, double n, double s=2.5);

double k_q(Vector3d &b, Vector3d &n_los, double nu, double n, double s=2.5);

// For random B-field
double k_q(double b, Vector3d &n_los, double nu, double n, double s=2.5);

double k_u(Vector3d &b, Vector3d &n_los, double nu, double n, double s=2.5);

// For random B-field
double k_u(double b, Vector3d &n_los, double nu, double n, double s=2.5);

double k_v(Vector3d &b, Vector3d &n_los, double nu, double n, double s=2.5);

// For random B-field
double k_v(double b, Vector3d &n_los, double nu, double n, double s=2.5);

double k_F_c(Vector3d &b, Vector3d &n_los, double nu, double n, double s=2.5);

// For random B-field
double k_F_c(double b, Vector3d &n_los, double nu, double n, double s=2.5);

double k_C_c(Vector3d &b, Vector3d &n_los, double nu, double n, double s=2.5);

// For random B-field
double k_C_c(double b, Vector3d &n_los, double nu, double n, double s=2.5);

// TODO: For e+ ``k_F_r`` will have different sign, but ``k_C_r`` - the same
// sign
double k_F_r(Vector3d &b, Vector3d &n_los, double nu, double n,
             double gamma_min=100., double s=2.5);

// For random B-field
double k_F_r(double b, Vector3d &n_los, double nu, double n,
             double gamma_min=100., double s=2.5);

double k_C_r(Vector3d &b, Vector3d &n_los, double nu, double n,
             double gamma_min=100., double s=2.5);

// For random B-field
double k_C_r(double b, Vector3d &n_los, double nu, double n,
             double gamma_min=100., double s=2.5);

double h_Q(Vector3d &b, Vector3d &n_los, double nu, double n, double s=2.5);

// For random B-field
double h_Q(double b, Vector3d &n_los, double nu, double n, double s=2.5);


double eta_0(Vector3d &b, Vector3d &n_los, double nu, double n);

// For random B-field
double eta_0(double b, Vector3d &n_los, double nu, double n);

// Emission coefficient for given vector of magnetic field ``b``, unit LOS
// vector ``n_los`` and others measured in emission frame
double eta_i(Vector3d &b, Vector3d &n_los, double nu, double n, double s=2.5);

// For random B-field
double eta_i(double b, Vector3d &n_los, double nu, double n, double s=2.5);

double eta_q(Vector3d &b, Vector3d &n_los, double nu, double n, double s=2.5);

// For random B-field
double eta_q(double b, Vector3d &n_los, double nu, double n, double s=2.5);

double eta_u(Vector3d &b, Vector3d &n_los, double nu, double n, double s=2.5);

// For random B-field
double eta_u(double b, Vector3d &n_los, double nu, double n, double s=2.5);

double eta_v(Vector3d &b, Vector3d &n_los, double nu, double n, double s=2.5);

// For random B-field
double eta_v(double b, Vector3d &n_los, double nu, double n, double s=2.5);


// Lorentz factor for the velocity ``v``
double getG(Vector3d &v);

// Doppler factor for LOS unit vector ``n_los`` in lab frame and bulk velocity
// ``v`` relative to the lab frame.
double getD(Vector3d &n_los, Vector3d &v);

// LOS unit vector in the comoving frame. This is relativistic aberration.
Vector3d get_n_los_prime(Vector3d &n_los, Vector3d &v);

// Vector of magnetic field in comoving (emission) frame in terms of the
// magnetic field ``b`` and velocity ``v`` in the observer frame.
Vector3d getBprime(Vector3d &b, Vector3d &v);

// This dummy function is needed to keep code of transport coefficients the same
// for both Vector and Scalar B-fields. It just returns ``b`` so it assumes that
// B-field is specified in plasma frame
double getBprime(double &b, Vector3d &v);


class Ctd {
public:
	//From astropy.cosmology.WMAP9
	Ctd(double z, double H0=69.32, double omega_M=0.2865, double omega_V=0.7134130719051658,
			double gamma_nu=8.69280948342326e-05);
	void operator() (const double &x, double &dxdt, const double t);
	double z;
	double H0;
	double omega_M;
	double omega_V;
	double gamma_nu;
};

// Given redshift ``z``, Hubble constant ``H_0`` [km/s/Mpc] and density
// parameters ``omega_M`` and ``omega_V``, returns comoving transverse distance
// (see arXiv:astro-ph/9905116v4 formula 14). Angular diameter distance is
// factor (1 + z) lower and luminosity distance is the same factor higher.
double comoving_transfer_distance(double z, double H0=69.32, double omega_M=0.2865, double omega_V=0.7134130719051658, double gamma_nu=8.69280948342326e-05);

double comoving_transfer_distance2(double z, double H0=69.32, double omega_M=0.2865, double omega_V=0.7134130719051658, double gamma_nu=8.69280948342326e-05);


double da_old(double z, double H0=50., double q0=0.05);


// Return scale factor that converts from parsecs to milliarcseconds
double pc_to_mas(double z);

// Return scale factor that converts from milliarcseconds to parsecs.
double mas_to_pc(double z);


// Return solid angle of one pixel. We need this to convert specific intensities
// I_{/nu} to flux densities S_{\nu}[Jy] = I_{\nu} * this
double pixel_solid_angle(double pixel_size_mas, double z);


std::ostream& write_2dvector(std::ostream& os,
														 std::vector<std::vector<double>>& v,
														 double scale=1.0);

std::ostream& write_vector(std::ostream& os,
                           std::vector<double>& v,
                           double scale=1.0);


std::vector<Vector3d> generate_random_directions(int n, unsigned int seed=0);


std::vector<Vector3d> generate_random_points_sphere(int n, double r_max,
                                                    double exponent,
                                                    unsigned int seed,
                                                    double r_min=0.0,
                                                    double theta_max=0.0,
                                                    double theta_min=0.0);

std::vector<Vector3d> generate_random_points_general(int n, double r_max,
                                                     Geometry* geo,
                                                     double exponent,
                                                     unsigned int seed);

int steps_schedule(double tau, int n_min, int n_max, double tau_min=0.1,
                   double tau_max=100.0);

//// return an evenly spaced 1-d grid of doubles.
//std::vector<double> linspace(double first, double last, int len) {
//	std::vector<double> result(len);
//	double step = (last-first) / (len - 1);
//	for (int i=0; i<len; i++) { result[i] = first + i*step; }
//	return result;
//}

void read_from_txt(std::string fn, std::vector< std::vector<double> >& properties);


// Linear interpolation following MATLAB linspace
std::vector<double> MyLinearSpacedArray(double a, double b, std::size_t N);

#endif //JETSHOW_UTILS_H
