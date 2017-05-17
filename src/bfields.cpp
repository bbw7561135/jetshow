#include <iostream>
#include <math.h>
#include <Eigen/Eigen>
#include <bfields.h>


// fi-value of magnetic field in some point ``p``
double BField::fiValue(Vector3d &p) {

    return fiValue0 * z0 / p[2];

}

// z-value of magnetic field in some point ``p``
double BField::zValue(Vector3d &p) {

    return zValue0 * (z0 / p[2]) * (z0 / p[2]);

}

// Vector of magnetic field at some point ``p``
Vector3d BField::bf(Vector3d &p) {

    double tmp = sqrt(p[0] * p[0] + p[1] * p[1]);
    return Eigen::Vector3d(-zValue(p) * p[1] / tmp, fiValue(p) * p[0] / tmp,
    zValue(p));
}

// Constructor
BField::BField(double newz0, double newfiValue0, double newzValue0) {
    z0 = newz0;
    fiValue0 = newfiValue0;
    zValue0 = newzValue0;
}
