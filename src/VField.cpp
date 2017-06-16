//
// Created by ilya on 6/6/17.
//

#include <math.h>
#include <Eigen/Eigen>
#include "VField.h"

using Eigen::Vector3d;


ConstFlatVField::ConstFlatVField(double gamma) : gamma_(gamma) {};

Vector3d ConstFlatVField::v(const Vector3d &point) const {
    return Vector3d(0, 0, sqrt(1. - 1./(gamma_*gamma_)));
}


ConstCentralVField::ConstCentralVField(double gamma) : gamma_(gamma) {}

Vector3d ConstCentralVField::v(const Vector3d &point) const {
    double x = point[0];
    double y = point[1];
    double z = point[2];
    double r = point.norm();
    double v_r = sqrt(1. - 1./(gamma_*gamma_));
    return Vector3d(v_r*x/r, v_r*y/r, v_r*z/r);
};

