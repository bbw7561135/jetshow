#include <utils.h>
#include "NField.h"


BKNField::BKNField (double n_0, double n_n) : n_0_(n_0), n_n_(n_n) {};


double BKNField::n(const Vector3d &point) const {
		double r = point.norm();
    return n_0_*pow(r/pc, -n_n_);
}


CompositeBKNField::CompositeBKNField(double n_0, double n_n_inner, double n_n_outer, double z0) :
    inner_field_(n_0, n_n_inner), outer_field_(n_0, n_n_outer), z0_(z0) {};


double CompositeBKNField::n(const Vector3d &point) const {
    double z = point[2];
    if (z > z0_) {
        return outer_field_.n(point);
    }
    else {
        return inner_field_.n(point);
    }
}


ConstNField::ConstNField(double n) : n_(n) {};

double ConstNField::n(const Vector3d &point) const {
	return n_;
}
