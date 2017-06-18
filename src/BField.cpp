#include <iostream>
#include <math.h>
#include <Eigen/Eigen>
#include <boost/math/special_functions/bessel.hpp>
#include <BField.h>
#include <utils.h>


using Eigen::Vector3d;


ConstCylinderBField::ConstCylinderBField(double b_0, double n_b) : b_0_(b_0),
                                                                   n_b_(n_b) {};

Vector3d ConstCylinderBField::bf(const Vector3d &point) const {
    return Vector3d(0.0, 0.0, b_0_*pow(point[2]/pc, -n_b_));
}



RadialConicalBField::RadialConicalBField(double b_0, double n_b) : b_0_(b_0),
                                                                   n_b_(n_b) {};

Vector3d RadialConicalBField::bf(const Vector3d &point) const {
    double r = point.norm();
    return Vector3d(b_0_*pow(r/pc, -n_b_)*point[0]/r,
                    b_0_*pow(r/pc, -n_b_)*point[1]/r,
                    b_0_*pow(r/pc, -n_b_)*point[2]/r);
}


HelicalCylinderBField::HelicalCylinderBField(double b_0, double pitch_angle) :
        b_0_(b_0), pitch_angle_(pitch_angle) {};

Vector3d HelicalCylinderBField::bf(const Vector3d &point) const {
    double r = sqrt(point[0]*point[0]+ point[1]*point[1]);
    return Vector3d(b_0_*tan(pitch_angle_*point[1]/r),
                    -b_0_*tan(pitch_angle_*point[0]/r),
                    b_0_);
}


SpiralConicalBField::SpiralConicalBField(double b_0, double pitch_angle) :
        b_0_(b_0), pitch_angle_(pitch_angle) {};

Vector3d SpiralConicalBField::bf(const Vector3d &point) const {
    double z = point[2];
    double x = point[0];
    double y = point[1];
    double b_z = b_0_/(z*z);
    return Vector3d(b_z*(x/z + y*tan(pitch_angle_)),
                    b_z*(y/z - x*tan(pitch_angle_)),
                    b_z);
}


ForceFreeCylindricalBField::ForceFreeCylindricalBField(double b_0, double mu) :
        b_0_(b_0), mu_(mu) {};

Vector3d ForceFreeCylindricalBField::bf(const Vector3d &point) const {
    double x = point[0];
    double y = point[1];
    double r = point.norm();
    double atan_term = atan(y/x);
    double bessel_0 = boost::math::cyl_bessel_i(0, mu_);
    double bessel_1 = boost::math::cyl_bessel_i(1, mu_);
    return Vector3d(-b_0_*bessel_1*sin(atan_term),
                    b_0_*bessel_1*sin(atan_term),
                    b_0_*bessel_0);
}



// This is old code for helical field. Note that it is likely to be wrong.
//// fi-value of magnetic field in some point ``p``
//double BField::fiValue(Vector3d &p) {
//
//    return fiValue0 * z0 / p[2];
//
//}
//
//// z-value of magnetic field in some point ``p``
//double BField::zValue(Vector3d &p) {
//
//    return zValue0 * (z0 / p[2]) * (z0 / p[2]);
//
//}
//
//// Vector of magnetic field at some point ``p``
//Vector3d BField::bf(Vector3d &p) {
//
//    double tmp = sqrt(p[0] * p[0] + p[1] * p[1]);
//    return Eigen::Vector3d(-zValue(p) * p[1] / tmp, fiValue(p) * p[0] / tmp,
//    zValue(p));
//}
//
//// Constructor
//BField::BField(double newz0, double newfiValue0, double newzValue0) {
//    z0 = newz0;
//    fiValue0 = newfiValue0;
//    zValue0 = newzValue0;
//}
