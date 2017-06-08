#ifndef JETSHOW_CILINDER_H
#define JETSHOW_CILINDER_H


#include <Eigen/Eigen>
#include "Geometry.h"

using Eigen::Vector3d;


class Cilinder : public Geometry {
public:

    Cilinder(const Vector3d &origin_, const Vector3d &direction_, double r_);

    Intersection hit(Ray &ray) const override;

private:
    Vector3d origin_;
    Vector3d direction_;
    double r_;
public:
    const Vector3d& origin() const;
    const Vector3d& direction() const;
    double r() const;

};


#endif //JETSHOW_CILINDER_H
