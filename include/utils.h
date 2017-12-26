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


const double mas_to_rad = 4.8481368*1E-9;
const double rad_to_mas = 1./mas_to_rad;
// Parsec [cm]
const double pc = 3.0857*1E18;
// Mass of electron [g]
const double m_e = 9.109382*1E-28;
// Mass of proton [g]
const double m_p = 1.672621*1E-24;
// Charge of electron [C]
 const double q_e = 4.8*1E-10;
//const double q_e = 1.6*1E-19;
// Charge of proton [C]
 const double q_p = 4.8*1E-10;
//const double q_p = 1.6*1E-19;
// Speed of light [cm / s]
const double c = 3.*1E+10;
// Jy in cgc
const double to_jy = 1E+23;
// pi
const double pi = boost::math::constants::pi<double>();


double nu_p(double n);

double nu_b(Vector3d &b, Vector3d &n_los);

double sin_theta(Vector3d &b, Vector3d &n_los);

double k_0(Vector3d &b, Vector3d &n_los, double nu, double n);

// Absorption coefficient for given vector of magnetic field ``b``, unit LOS
// vector ``n_los`` and others measured in emission frame
double k_i(Vector3d &b, Vector3d &n_los, double nu, double n, double s=2.5);

double k_q(Vector3d &b, Vector3d &n_los, double nu, double n, double s=2.5);

double k_u(Vector3d &b, Vector3d &n_los, double nu, double n, double s=2.5);

double k_v(Vector3d &b, Vector3d &n_los, double nu, double n, double s=2.5);

double k_F_c(Vector3d &b, Vector3d &n_los, double nu, double n, double s=2.5);

double k_C_c(Vector3d &b, Vector3d &n_los, double nu, double n, double s=2.5);

// TODO: For e+ ``k_F_r`` will have different sign, but ``k_C_r`` - the same
// sign
double k_F_r(Vector3d &b, Vector3d &n_los, double nu, double n,
             double gamma_min=100., double s=2.5);

double k_C_r(Vector3d &b, Vector3d &n_los, double nu, double n,
             double gamma_min=100., double s=2.5);

double h_Q(Vector3d &b, Vector3d &n_los, double nu, double n, double s=2.5);



double eta_0(Vector3d &b, Vector3d &n_los, double nu, double n);

// Emission coefficient for given vector of magnetic field ``b``, unit LOS
// vector ``n_los`` and others measured in emission frame
double eta_i(Vector3d &b, Vector3d &n_los, double nu, double n, double s=2.5);

double eta_q(Vector3d &b, Vector3d &n_los, double nu, double n, double s=2.5);

double eta_u(Vector3d &b, Vector3d &n_los, double nu, double n, double s=2.5);

double eta_v(Vector3d &b, Vector3d &n_los, double nu, double n, double s=2.5);


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


class Ctd {
public:
		Ctd(double z, double H0=73., double omega_M=0.3, double omega_V=0.7);
		void operator() (const double &x, double &dxdt, const double t);
		double z;
		double H0;
		double omega_M;
		double omega_V;
};

// Given redshift ``z``, Hubble constant ``H_0`` [km/s/Mpc] and density
// parameters ``omega_M`` and ``omega_V``, returns comoving transverse distance
// (see arXiv:astro-ph/9905116v4 formula 14). Angular diameter distance is
// factor (1 + z) lower and luminosity distance is the same factor higher.
double comoving_transfer_distance(double z, double H0=73., double omega_M=0.3,
																  double omega_V=0.7);


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


#endif //JETSHOW_UTILS_H
