#include <math.h>
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
    return q_e*(n_los.cross(b)).norm()/(2.*pi*m_e*c*b.norm());
}

double k_0(Vector3d &b, Vector3d &n_los, double nu, double n) {
    return pi*nu_p(n)*nu_p(n)*nu_b(b, n_los)/(c*nu*nu);
}

double k_i(Vector3d &b, Vector3d &n_los, double nu, double n, double s) {
    double factor = pow(3., (s+1.)/2.)/(4.)*tgamma(s/4.+11./6.)*tgamma(s/4.+1./6.);
    return k_0(b, n_los, nu, n) * pow(nu_b(b, n_los)/nu, s/2.) * factor;
}

double eta_0(Vector3d &b, Vector3d &n_los, double n) {
    return pi*nu_p(n)*nu_p(n)*nu_b(b, n_los)*m_e/c;
}

double eta_i(Vector3d &b, Vector3d &n_los, double nu, double n, double s) {
    double factor = pow(3., s/2.)/(2.*(s+1))*tgamma(s/4.+19./12.)*tgamma(s/4.-0.5);
    return eta_0(b, n_los, n) * pow(nu_b(b, n_los)/nu, (s-1.)/2.) * factor;
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

