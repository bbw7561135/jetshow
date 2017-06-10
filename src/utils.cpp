#include <math.h>
#include "utils.h"


double nu_p(double n) {
    return sqrt(n*q_e*q_e / (pi*m_e));
}

double nu_b(double b) {
    return q_e*b/(2.*pi*m_e*c);
}

double eta_0(double n, double b) {
    return  pi*nu_p(n)*nu_p(n)*nu_b(b)*m_e/c;
}

double k_0(double nu, double n, double b) {
    return pi*nu_p(n)*nu_p(n)*nu_b(b)/(c*nu*nu);
}

double eta_i(double nu, double n, double b, double sin_theta,
             double s) {
    return eta_0(n, b)*sin_theta*pow(nu_b(b)*sin_theta/nu, (s-1.)/2.)*pow(3., s/2.)/(2.*(s+1.))*tgamma(s/4.+19./12.)*tgamma(s/4.-1./12.);
}