#ifndef JETSHOW_NFIELD_H
#define JETSHOW_NFIELD_H

#include <Eigen/Eigen>

using Eigen::Vector3d;

class NField {
public:
    virtual double n(const Vector3d& point) const = 0;
};

class ConstNField: public NField {
public:
		explicit ConstNField(double n);
		double n(const Vector3d &point) const override;

private:
		double n_;
};

class BKNField: public NField {

    public:
    BKNField(double n_0, double n_n);

    double n(const Vector3d &point) const override;

private:
    double n_0_;
    double n_n_;

};


class CompositeBKNField: public NField {
public:
    CompositeBKNField(double n_0, double n_n_inner, double n_n_outer, double z0);
    double n(const Vector3d &point) const override;

private:
    double z0_;
    BKNField inner_field_;
    BKNField outer_field_;
};

#endif //JETSHOW_NFIELD_H
