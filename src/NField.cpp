#include <utils.h>
#include "NField.h"


NField::NField (bool in_plasma_frame, double s) :
    in_plasma_frame_(in_plasma_frame),
    s_(s) {}


double NField::nf_plasma_frame(const Vector3d &point, double &gamma) const {
    double n = nf(point);
    if (in_plasma_frame_) {
        return n;
    } else {
        return n/gamma;
    }
}


double NField::get_s() const {
    return s_;
}


ConstNField::ConstNField(double n, bool in_plasma_frame, double s) :
    NField(in_plasma_frame, s),
    n_(n) {}


double ConstNField::nf(const Vector3d &point) const {
    return n_;
}


BKNField::BKNField (double n_0, double n_n, bool in_plasma_frame, double s) :
    NField(in_plasma_frame, s),
    n_0_(n_0),
    n_n_(n_n) {}


double BKNField::nf(const Vector3d &point) const {
    double r = point.norm();
    return n_0_*pow(r/pc, -n_n_);
}


BKNFieldZ::BKNFieldZ (double n_0, double n_n, bool in_plasma_frame, double s) :
        NField(in_plasma_frame, s),
        n_0_(n_0),
        n_n_(n_n) {}


double BKNFieldZ::nf(const Vector3d &point) const {
    double z = abs(point[2]);
    return n_0_*pow(z/pc, -n_n_);
}


CompositeBKNField::CompositeBKNField(double n_0, double n_n_inner, double n_n_outer, double z0, bool in_plasma_frame,
    double s) :
    NField(in_plasma_frame, s),
    inner_field_(n_0, n_n_inner, in_plasma_frame, s),
    outer_field_(n_0, n_n_outer, in_plasma_frame, s),
    z0_(z0) {}


double CompositeBKNField::nf(const Vector3d &point) const {
    double z = abs(point[2]);
    if (z > z0_) {
        return outer_field_.nf(point);
    }
    else {
        return inner_field_.nf(point);
    }
}


SimulationNField::SimulationNField(Delaunay_triangulation *tr, bool in_plasma_frame, double s) :
    NField(in_plasma_frame, s),
    interp_(tr) {}

double SimulationNField::nf(const Vector3d &point) const {
    return interp_.interpolated_value(point);
};