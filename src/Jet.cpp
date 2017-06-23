#include "Jet.h"


Jet::Jet(Geometry *newgeo, VField *newvfield, BField *newbField,
         NField *newnField) {
    geometry_ = newgeo;
    vfield_ = newvfield;
    bfield_ = newbField;
    nfield_ = newnField;
}

// This is k_i in lab frame that could be integrated along LOS.
double Jet::getKI(Vector3d &point, Vector3d &n_los, double nu) {
    // First, comoving frame ``k_i_prime`` (in the rest frame of the emission
    // element) is connected to this ``k_i`` as ``k_i = k_i_prime / D``.
    // Second, in ``k_i_prime`` we need all quantities in comoving frame
    // (primed) in terms of lab frame:
    // b_prime = f(b, v)
    // n_los_prime = f(n_los, v)
    // nu_prime = f(nu, n_los, v) = nu/getD
    // n_prime = f(n, v) = n/Gamma
    Vector3d b = getB(point);
    Vector3d v = getV(point);
    auto D = getD(n_los, v);
    auto gamma = getG(v);
    double n = getN(point);
    auto b_prime = getBprime(b, v);
    auto n_los_prime = get_n_los_prime(n_los, v);
    auto nu_prime = nu/D;
    auto n_prime = n/gamma;
    auto k_i_prime = k_i(b_prime, n_los_prime, nu_prime, n_prime);

		//
//		std::cout << "b = " << b.norm() << " v = " << v.norm() << " G = " << gamma <<
//							" D = " << D << " n = " << n << " n_pr = " << n_prime <<
//							" b_pr= " << b_prime.norm() << std::endl;

    return k_i_prime/D;
}

// This is k_i in lab frame that could be integrated along LOS.
double Jet::getEtaI(Vector3d &point, Vector3d &n_los, double nu) {
    // First, comoving frame ``eta_i_prime`` (in the rest frame of the emission
    // element) is connected to this ``eta_i`` as ``eta_i = D^2 * eta_i_prime``.
    // Second, in ``eta_i_prime`` we need all quantities in comoving frame
    // (primed) in terms of lab frame:
    // b_prime = f(b, v)
    // n_los_prime = f(n_los, v)
    // nu_prime = f(nu, n_los, v) = nu/getD
    // n_prime = f(n, v) = n/Gamma
    Vector3d b = getB(point);
    Vector3d v = getV(point);
    auto D = getD(n_los, v);
    auto gamma = getG(v);
    double n = getN(point);
    auto b_prime = getBprime(b, v);
    auto n_los_prime = get_n_los_prime(n_los, v);
    auto nu_prime = nu/D;
    auto n_prime = n/gamma;
    auto eta_i_prime = eta_i(b_prime, n_los_prime, nu_prime, n_prime);
//		std::cout << "b = " << b.norm() << " v = " << v.norm() << " G = " << gamma <<
//							" D = " << D << " n = " << n << " n_pr = " << n_prime <<
//							" b_pr= " << b_prime.norm() << std::endl;
//		std::cout << "n_los in jet.getEtaI " << n_los << std::endl;
//		std::cout << "n_los_prime in jet.getEtaI " << n_los_prime << std::endl;
    return eta_i_prime*D*D;
}

const Vector3d Jet::getB(const Vector3d &point) {
    return bfield_->bf(point);
}

const Vector3d Jet::getV(const Vector3d &point) {
    return vfield_->v(point);
}

const double Jet::getN(const Vector3d &point) {
    return nfield_->n(point);
}


std::list<Intersection> Jet::hit(Ray &ray) {
  return geometry_->hit(ray);
}
