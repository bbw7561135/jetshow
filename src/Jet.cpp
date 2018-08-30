#include "Jet.h"
#include <math.h>


//AnalyticalScalarBFieldJet::AnalyticalScalarBFieldJet(Geometry *newgeo, VField *newvfield, RandomScalarBField *newbField,
//         NField *newnField) {
//    geometry_ = newgeo;
//    vfield_ = newvfield;
//    bfield_ = newbField;
//    nfield_ = newnField;
//}

Jet::Jet(BaseGeometry *newgeo, VField *newvfield, VectorBField *newbField,
         NField *newnField) {
    geometry_ = newgeo;
    vfield_ = newvfield;
    bfield_ = newbField;
    nfield_ = newnField;
}

//Jet::Jet(Geometry *newgeo, VField *newvfield, VectorBField *newbField,
//         NField *newnField) {
//    geometry_ = newgeo;
//    vfield_ = newvfield;
//    bfield_ = newbField;
//    nfield_ = newnField;
//}


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

    Vector3d v = getV(point);
    auto gamma = getG(v);
    auto b_prime = bfield_->bf_plasma_frame(point, v);
    auto D = getD(n_los, v);
    double n_prime = nfield_->nf_plasma_frame(point, gamma);
    auto n_los_prime = get_n_los_prime(n_los, v);
    auto nu_prime = nu/D;
    auto k_i_prime = k_i(b_prime, n_los_prime, nu_prime, n_prime);

    if (std::isnan(k_i_prime)) {
        std::cout << "Nan k_i_prime!" << std::endl;
        std::cout << "b_prime = " << b_prime << std::endl;
        std::cout << "n_prime = " << n_prime << std::endl;
        std::cout << "nu_prime = " << nu_prime << std::endl;
        std::cout << "n_los_prime = " << n_los_prime << std::endl;
        std::cout << "D = " << D << std::endl;
        std::cout << "v = " << v << std::endl;
        std::cout << "gamma = " << gamma << std::endl;
    }
    return k_i_prime/D;
}


//// This is k_i in lab frame that could be integrated along LOS.
//double AnalyticalScalarBFieldJet::getKI(Vector3d &point, Vector3d &n_los, double nu) {
//    // First, comoving frame ``k_i_prime`` (in the rest frame of the emission
//    // element) is connected to this ``k_i`` as ``k_i = k_i_prime / D``.
//    // Second, in ``k_i_prime`` we need all quantities in comoving frame
//    // (primed) in terms of lab frame:
//    // b_prime = f(b, v)
//    // n_los_prime = f(n_los, v)
//    // nu_prime = f(nu, n_los, v) = nu/getD
//    // n_prime = f(n, v) = n/Gamma
//
//    Vector3d v = getV(point);
//    auto D = getD(n_los, v);
//    auto gamma = getG(v);
//    double n = getN(point);
//    // This means that we are using B-field specification in the plasma frame
//    auto b_prime = bfield_->bf(point);
//    auto n_los_prime = get_n_los_prime(n_los, v);
//    auto nu_prime = nu/D;
////  	auto n_prime = n/gamma;
//    // This means that we are now using n specification in plasma frame
//    auto n_prime = n;
//    auto k_i_prime = k_i(b_prime, n_los_prime, nu_prime, n_prime);
//
//    if (std::isnan(k_i_prime)) {
//        std::cout << "Nan k_i_prime!" << std::endl;
//        std::cout << "b_prime = " << b_prime << std::endl;
//        std::cout << "n_prime = " << n_prime << std::endl;
//        std::cout << "nu_prime = " << nu_prime << std::endl;
//        std::cout << "n_los_prime = " << n_los_prime << std::endl;
//        std::cout << "D = " << D << std::endl;
//        std::cout << "v = " << v << std::endl;
//        std::cout << "gamma = " << gamma << std::endl;
//    }
//    return k_i_prime/D;
//}
//
//
//// This is eta_i in lab frame that could be integrated along LOS.
//double AnalyticalScalarBFieldJet::getEtaI(Vector3d &point, Vector3d &n_los, double nu) {
//    // First, comoving frame ``eta_i_prime`` (in the rest frame of the emission
//    // element) is connected to this ``eta_i`` as ``eta_i = D^2 * eta_i_prime``.
//    // Second, in ``eta_i_prime`` we need all quantities in comoving frame
//    // (primed) in terms of lab frame:
//    // b_prime = f(b, v)
//    // n_los_prime = f(n_los, v)
//    // nu_prime = f(nu, n_los, v) = nu/getD
//    // n_prime = f(n, v) = n/Gamma
//
//    Vector3d v = getV(point);
//    auto D = getD(n_los, v);
//    auto gamma = getG(v);
//    double n = getN(point);
//    // This means that we are using B-field specification in the plasma frame
//    auto b_prime = bfield_->bf(point);
//    auto n_los_prime = get_n_los_prime(n_los, v);
//    auto nu_prime = nu/D;
////    auto n_prime = n/gamma;
//    // This means that we are now using n specification in plasma frame
//    auto n_prime = n;
//    auto eta_i_prime = eta_i(b_prime, n_los_prime, nu_prime, n_prime);
//    return eta_i_prime*D*D;
//}


//double Jet::getKQ(Vector3d &point, Vector3d &n_los, double nu) {
//	// This conflicts with scalar valued random B-field
////	Vector3d b = getB(point);
//	Vector3d v = getV(point);
//	auto D = getD(n_los, v);
//	auto gamma = getG(v);
//	double n = getN(point);
////	auto b_prime = getBprime(b, v);
//	// This means that we are using B-fields specification in the plasma frame
//	auto b_prime = bfield_->bf(point);
//	auto n_los_prime = get_n_los_prime(n_los, v);
//	auto nu_prime = nu/D;
//	auto n_prime = n/gamma;
//	auto k_q_prime = k_q(b_prime, n_los_prime, nu_prime, n_prime);
//	return k_q_prime/D;
//}
//
//double Jet::getKU(Vector3d &point, Vector3d &n_los, double nu) {
//	// This conflicts with scalar valued random B-field
////	Vector3d b = getB(point);
//	Vector3d v = getV(point);
//	auto D = getD(n_los, v);
//	auto gamma = getG(v);
//	double n = getN(point);
////	auto b_prime = getBprime(b, v);
//	// This means that we are using B-fields specification in the plasma frame
//	auto b_prime = bfield_->bf(point);
//	auto n_los_prime = get_n_los_prime(n_los, v);
//	auto nu_prime = nu/D;
//	auto n_prime = n/gamma;
//	auto k_u_prime = k_u(b_prime, n_los_prime, nu_prime, n_prime);
//	return k_u_prime/D;
//}
//
//double Jet::getKV(Vector3d &point, Vector3d &n_los, double nu) {
//	// This conflicts with scalar valued random B-field
////	Vector3d b = getB(point);
//	Vector3d v = getV(point);
//	auto D = getD(n_los, v);
//	auto gamma = getG(v);
//	double n = getN(point);
////	auto b_prime = getBprime(b, v);
//	// This means that we are using B-fields specification in the plasma frame
//	auto b_prime = bfield_->bf(point);
//	auto n_los_prime = get_n_los_prime(n_los, v);
//	auto nu_prime = nu/D;
//	auto n_prime = n/gamma;
//	auto k_v_prime = k_v(b_prime, n_los_prime, nu_prime, n_prime);
//	return k_v_prime/D;
//}
//
//double Jet::getKF(Vector3d &point, Vector3d &n_los, double nu) {
//	// This conflicts with scalar valued random B-field
////	Vector3d b = getB(point);
//	Vector3d v = getV(point);
//	auto D = getD(n_los, v);
//	auto gamma = getG(v);
//	double n = getN(point);
////	auto b_prime = getBprime(b, v);
//	// This means that we are using B-fields specification in the plasma frame
//	auto b_prime = bfield_->bf(point);
//	auto n_los_prime = get_n_los_prime(n_los, v);
//	auto nu_prime = nu/D;
//	auto n_prime = n/gamma;
//	auto k_F_prime = k_F_r(b_prime, n_los_prime, nu_prime, n_prime);
//	return k_F_prime/D;
//}
//
//double Jet::getKC(Vector3d &point, Vector3d &n_los, double nu) {
//	// This conflicts with scalar valued random B-field
////	Vector3d b = getB(point);
//	Vector3d v = getV(point);
//	auto D = getD(n_los, v);
//	auto gamma = getG(v);
//	double n = getN(point);
////	auto b_prime = getBprime(b, v);
//	// This means that we are using B-fields specification in the plasma frame
//	auto b_prime = bfield_->bf(point);
//	auto n_los_prime = get_n_los_prime(n_los, v);
//	auto nu_prime = nu/D;
//	auto n_prime = n/gamma;
//	auto k_C_prime = k_C_r(b_prime, n_los_prime, nu_prime, n_prime);
//	return k_C_prime/D;
//}
//
//double Jet::gethQ(Vector3d &point, Vector3d &n_los, double nu) {
//	// This conflicts with scalar valued random B-field
////	Vector3d b = getB(point);
//	Vector3d v = getV(point);
//	auto D = getD(n_los, v);
//	auto gamma = getG(v);
//	double n = getN(point);
////	auto b_prime = getBprime(b, v);
//	// This means that we are using B-fields specification in the plasma frame
//	auto b_prime = bfield_->bf(point);
//	auto n_los_prime = get_n_los_prime(n_los, v);
//	auto nu_prime = nu/D;
//	auto n_prime = n/gamma;
//	auto h_Q_prime = h_Q(b_prime, n_los_prime, nu_prime, n_prime);
//	return h_Q_prime/D;
//}


// This is eta_i in lab frame that could be integrated along LOS.
double Jet::getEtaI(Vector3d &point, Vector3d &n_los, double nu) {
    // First, comoving frame ``eta_i_prime`` (in the rest frame of the emission
    // element) is connected to this ``eta_i`` as ``eta_i = D^2 * eta_i_prime``.
    // Second, in ``eta_i_prime`` we need all quantities in comoving frame
    // (primed) in terms of lab frame:
    // b_prime = f(b, v)
    // n_los_prime = f(n_los, v)
    // nu_prime = f(nu, n_los, v) = nu/getD
    // n_prime = f(n, v) = n/Gamma
    Vector3d v = getV(point);
    auto gamma = getG(v);
    auto b_prime = bfield_->bf_plasma_frame(point, v);
    auto D = getD(n_los, v);
    double n_prime = nfield_->nf_plasma_frame(point, gamma);
    auto n_los_prime = get_n_los_prime(n_los, v);
    auto nu_prime = nu/D;

    auto eta_i_prime = eta_i(b_prime, n_los_prime, nu_prime, n_prime);
    return eta_i_prime*D*D;
}

//double Jet::getEtaQ(Vector3d &point, Vector3d &n_los, double nu) {
//	// This conflicts with scalar valued random B-field
////	Vector3d b = getB(point);
//	Vector3d v = getV(point);
//	auto D = getD(n_los, v);
//	auto gamma = getG(v);
//	double n = getN(point);
////	auto b_prime = getBprime(b, v);
//	// This means that we are using B-fields specification in the plasma frame
//	auto b_prime = bfield_->bf(point);
//	auto n_los_prime = get_n_los_prime(n_los, v);
//	auto nu_prime = nu/D;
//	auto n_prime = n/gamma;
//	auto eta_q_prime = eta_q(b_prime, n_los_prime, nu_prime, n_prime);
//	return eta_q_prime*D*D;
//}
//
//double Jet::getEtaU(Vector3d &point, Vector3d &n_los, double nu) {
//	// This conflicts with scalar valued random B-field
////	Vector3d b = getB(point);
//	Vector3d v = getV(point);
//	auto D = getD(n_los, v);
//	auto gamma = getG(v);
//	double n = getN(point);
////	auto b_prime = getBprime(b, v);
//	// This means that we are using B-fields specification in the plasma frame
//	auto b_prime = bfield_->bf(point);
//	auto n_los_prime = get_n_los_prime(n_los, v);
//	auto nu_prime = nu/D;
//	auto n_prime = n/gamma;
//	auto eta_u_prime = eta_u(b_prime, n_los_prime, nu_prime, n_prime);
//	return eta_u_prime*D*D;
//}
//
//double Jet::getEtaV(Vector3d &point, Vector3d &n_los, double nu) {
//	// This conflicts with scalar valued random B-field
////	Vector3d b = getB(point);
//	Vector3d v = getV(point);
//	auto D = getD(n_los, v);
//	auto gamma = getG(v);
//	double n = getN(point);
////	auto b_prime = getBprime(b, v);
//	// This means that we are using B-fields specification in the plasma frame
//	auto b_prime = bfield_->bf(point);
//	auto n_los_prime = get_n_los_prime(n_los, v);
//	auto nu_prime = nu/D;
//	auto n_prime = n/gamma;
//	auto eta_v_prime = eta_v(b_prime, n_los_prime, nu_prime, n_prime);
//	return eta_v_prime*D*D;
//}


const Vector3d Jet::getB(const Vector3d &point) {
    return bfield_->bf(point);
}



const Vector3d Jet::getV(const Vector3d &point) {
    return vfield_->v(point);
}

const double Jet::getN(const Vector3d &point) {
    return nfield_->nf(point);
}


std::list<Intersection> Jet::hit(Ray &ray) {
  return geometry_->hit(ray);
}



