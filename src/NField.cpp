#include <utils.h>
#include "NField.h"


BKNField::BKNField (double n_0, double n_n) : n_0_(n_0), n_n_(n_n) {};


double BKNField::n(const Vector3d &point) const {
		double r = point.norm();
    return n_0_*pow(r/pc, -n_n_);
}

ConstNField::ConstNField(double n) : n_(n) {};

double ConstNField::n(const Vector3d &point) const {
	return n_;
}
