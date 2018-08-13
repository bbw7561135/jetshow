//
// Created by ilya on 13.08.18.
//

#ifndef JETSHOW_PARABALOID_H
#define JETSHOW_PARABALOID_H


#include <Eigen/Eigen>
#include "Geometry.h"

using Eigen::Vector3d;


class Parabaloid : public Geometry {
public:

    Parabaloid(const Vector3d &origin_, const Vector3d &direction_,
             const double &r0_);

    std::list<Intersection> hit(Ray &ray) const override;
    const Vector3d& origin() const override;
    const Vector3d& direction() const;
    double r0() const;
    bool is_within(Vector3d& point) const override;
    const double big_scale() const override ;

private:
    Vector3d origin_;
    Vector3d direction_;
    double r0_;
};


#endif //JETSHOW_PARABALOID_H
