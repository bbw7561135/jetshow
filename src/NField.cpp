//
// Created by ilya on 6/6/17.
//

#include "NField.h"


BKNField::BKNField (double z_0, double n_0) : z_0_(z_0), n_0_(n_0) {};


double BKNField::n(const Vector3d &point) const {
    double z = point[2];
    return n_0_*(z_0_/z)*(z_0_/z);
}
