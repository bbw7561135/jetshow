#ifndef JETSHOW_NFIELD_H
#define JETSHOW_NFIELD_H

#include <Eigen/Eigen>
#include "SimulationInterpolater.h"


using Eigen::Vector3d;

class NField {
public:
    virtual double nf(const Vector3d &point) const = 0;
    double nf_plasma_frame(const Vector3d &point, double &gamma) const;

protected:
    explicit NField(bool in_plasma_frame);
    bool in_plasma_frame_;
};

class ConstNField: public NField {
public:
    ConstNField(double n, bool in_plasma_frame);
    double nf(const Vector3d &point) const override;
private:
    double n_;
};

class BKNField: public NField {
public:
    BKNField(double n_0, double n_n, bool in_plasma_frame);
    double nf(const Vector3d &point) const override;
private:
    double n_0_;
    double n_n_;

};


class CompositeBKNField: public NField {
public:
    CompositeBKNField(double n_0, double n_n_inner, double n_n_outer, double z0, bool in_plasma_frame);
    double nf(const Vector3d &point) const override;

private:
    double z0_;
    BKNField inner_field_;
    BKNField outer_field_;
};


class SimulationNField: public NField {
public:
    SimulationNField(Delaunay_triangulation *tr, bool in_plasma_frame);
    double nf(const Vector3d &point) const override;

private:
    SimulationInterpolater interp_;
};

#endif //JETSHOW_NFIELD_H
