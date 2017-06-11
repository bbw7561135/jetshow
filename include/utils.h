#ifndef JETSHOW_UTILS_H
#define JETSHOW_UTILS_H

#include <math.h>
#include <Eigen/Eigen>
#include <boost/math/constants/constants.hpp>

using Eigen::Vector3d;


const double mas_to_rad = 4.8481368*pow(10., -9.);
const double rad_to_mas = 1./mas_to_rad;
// Parsec [cm]
const double pc = 3.0857*pow(10., 18.);
// Mass of electron [g]
const double m_e = 9.109382*pow(10., -28.);
// Mass of proton [g]
const double m_p = 1.672621*pow(10., -24.);
// Charge of electron [C]
// const double q_e = 1.602176*pow(10.,-19.);
const double q_e = 4.8*pow(10.,-10.);
// Charge of proton [C]
const double q_p = 4.8*pow(10.,-10.);
// Speed of light [cm / s]
const double c = 3.*pow(10.,10.);
// Jy in cgc
const double to_jy = pow(10.,23.);
// pi
const double pi = boost::math::constants::pi<double>();


double nu_p(double n);

double nu_b(Vector3d &b, Vector3d &n_los);;

double k_0(Vector3d &b, Vector3d &n_los, double nu, double n);

// Absorption coefficient for given vector of magnetic field ``b``, unit LOS
// vector ``n_los`` and others measured in emission frame
double k_i(Vector3d &b, Vector3d &n_los, double nu, double n, double s=2.5);

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

#endif //JETSHOW_UTILS_H
