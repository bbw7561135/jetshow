#include <math.h>
#include <Eigen/Eigen>
#include <utils.h>
#include "VField.h"

using Eigen::Vector3d;


ConstFlatVField::ConstFlatVField(double gamma) : gamma_(gamma) {};

Vector3d ConstFlatVField::v(const Vector3d &point) const {
    return Vector3d(0, 0, c*sqrt(1. - 1./(gamma_*gamma_)));
}


ShearedFlatVField::ShearedFlatVField(double gamma_axis, double gamma_border,
                                     double r) : gamma_axis_(gamma_axis),
                                                 gamma_border_(gamma_border),
                                                 r_(r) {}

Vector3d ShearedFlatVField::v(const Vector3d &point) const {
	double x = point[0];
	double y = point[1];
	double r = sqrt(x*x + y*y);
	double gamma = gamma_axis_-(gamma_axis_-gamma_border_)*r/r_;
	return Vector3d(0, 0, c*sqrt(1. - 1./(gamma*gamma)));
}


SheathFlatVField::SheathFlatVField(double gamma_spine, double gamma_sheath,
                                   double r_sheath) :
		gamma_spine_(gamma_spine), gamma_sheath_(gamma_sheath),
		r_sheath_(r_sheath) {}

Vector3d SheathFlatVField::v(const Vector3d &point) const {
	double x = point[0];
	double y = point[1];
	double r = sqrt(x*x + y*y);
	double gamma;
	if (r < r_sheath_) {
		gamma = gamma_spine_;
	} else {
		gamma = gamma_sheath_;
	}
	return Vector3d(0, 0, c*sqrt(1. - 1./(gamma*gamma)));
}


ConstCentralVField::ConstCentralVField(double gamma) : gamma_(gamma) {}

Vector3d ConstCentralVField::v(const Vector3d &point) const {
    double x = point[0];
    double y = point[1];
    double z = point[2];
    double r = point.norm();
    double v_r = c*sqrt(1. - 1./(gamma_*gamma_));
    return Vector3d(v_r*x/r, v_r*y/r, v_r*z/r);
};


ShearedCentralVField::ShearedCentralVField(double gamma_axis,
                                           double gamma_border, double theta) :
		gamma_axis_(gamma_axis), gamma_border_(gamma_border), theta_(theta) {}

Vector3d ShearedCentralVField::v(const Vector3d &point) const {
	double x = point[0];
	double y = point[1];
	double z = point[2];
	double r = point.norm();
	double theta = acos(z/r);
	double gamma = gamma_axis_-(gamma_axis_-gamma_border_)*theta/theta_;
	double v_r = c*sqrt(1. - 1./(gamma*gamma));
	return Vector3d(v_r*x/r, v_r*y/r, v_r*z/r);
}


SheathCentralVField::SheathCentralVField(double gamma_spine,
                                         double gamma_sheath,
                                         double theta_sheath) :
		gamma_spine_(gamma_spine), gamma_sheath_(gamma_sheath),
		theta_sheath_(theta_sheath) {}

Vector3d SheathCentralVField::v(const Vector3d &point) const {
	double x = point[0];
	double y = point[1];
	double z = point[2];
	double r = point.norm();
	double theta = acos(z/r);
	double gamma;
	if (theta < theta_sheath_) {
		gamma = gamma_spine_;
	} else {
		gamma = gamma_sheath_;
	}
	double v_r = c*sqrt(1. - 1./(gamma*gamma));
	return Vector3d(v_r*x/r, v_r*y/r, v_r*z/r);
}


ConstParabolicVField::ConstParabolicVField(double gamma, double R0) : gamma_(gamma), R0_(R0) {}

Vector3d ConstParabolicVField::v(const Vector3d &point) const {
    double x = point[0];
    double y = point[1];
    double r = std::hypot(x, y);
    double v_r = c*sqrt(1. - 1./(gamma_*gamma_));
    double alpha = atan(2*r/sqrt(R0_));
    return Vector3d(v_r*cos(alpha)*x/r, v_r*cos(alpha)*y/r, v_r*sin(alpha));
};


AccParabolicVField::AccParabolicVField(double gamma0, double R0) : gamma0_(gamma0), R0_(R0) {}

Vector3d AccParabolicVField::v(const Vector3d &point) const {
    double x = point[0];
    double y = point[1];
    double z = point[2];
    double r = std::hypot(x, y);
    double gamma = gamma0_*sqrt(z);
    double v_r = c*sqrt(1. - 1./(gamma*gamma));
    double alpha = atan(2*r/sqrt(R0_));
    return Vector3d(v_r*cos(alpha)*x/r, v_r*cos(alpha)*y/r, v_r*sin(alpha));
};


ShearedAccParabolicVField::ShearedAccParabolicVField(double gamma_axis0, double gamma_border0, double R0,
        double R0_border) :
        gamma_axis0_(gamma_axis0), gamma_border0_(gamma_border0), R0_(R0), R0_border_(R0_border) {}

Vector3d ShearedAccParabolicVField::v(const Vector3d &point) const {
    double x = point[0];
    double y = point[1];
    double z = point[2];
    double r = std::hypot(x, y);
    double gamma0 = gamma_axis0_-(gamma_axis0_-gamma_border0_)*r/R0_border_;
    double gamma = gamma0*sqrt(z);
    double v_r = c*sqrt(1. - 1./(gamma*gamma));
    double alpha = atan(2*r/sqrt(R0_));
    return Vector3d(v_r*cos(alpha)*x/r, v_r*cos(alpha)*y/r, v_r*sin(alpha));
}
//ConstParabolicConstConeVField::ConstParabolicConstConeVField(double gamma, double R0, double z0) :
//z0_(z0), conev(gamma), parav(gamma, R0) {}
//
//Vector3d ConstParabolicConstConeVField::v(const Vector3d &point) const {
//    if (point[2] > z0_) {
//        return conev.v(point);
//    }
//    else {
//        return parav.v(point);
//    }
//};


AccParabolicConstConeVField::AccParabolicConstConeVField(double gamma0, double R0, double z0) :
        z0_(z0), conev(gamma0), parav(gamma0, R0) {}

Vector3d AccParabolicConstConeVField::v(const Vector3d &point) const {
    if (point[2] > z0_) {
        return conev.v(point);
    }
    else {
        return parav.v(point);
    }
};


ShearedAccParabolicConstConeVField::ShearedAccParabolicConstConeVField(double gamma_axis0, double gamma_border0,
        double R0, double R0_border, double z0) :
        z0_(z0),
        conev(gamma_axis0, gamma_border0, acos(z0/R0_border)),
        parav(gamma_axis0, gamma_border0, R0, R0_border) {}

Vector3d ShearedAccParabolicConstConeVField::v(const Vector3d &point) const {
    if (point[2] > z0_) {
        return conev.v(point);
    }
    else {
        return parav.v(point);
    }
};