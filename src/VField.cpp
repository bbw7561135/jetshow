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
