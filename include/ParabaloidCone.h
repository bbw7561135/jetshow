//
// Created by ilya on 13.08.18.
//

#ifndef JETSHOW_PARABALOIDCONE_H
#define JETSHOW_PARABALOIDCONE_H

#include <Eigen/Eigen>
#include "Geometry.h"
#include "Cone.h"
#include "Parabaloid.h"

using Eigen::Vector3d;


class ParabaloidCone : public Geometry {
public:

    ParabaloidCone(const Vector3d &origin_, const Vector3d &direction_,
                   const double &r0_, const double &z0_, const double &scale_);

    std::list<Intersection> hit(Ray &ray) const override;
    const Vector3d& origin() const override;
    const Vector3d& direction() const;
    double z0() const;
    bool is_within(Vector3d& point) const override;
    const double big_scale() const override ;


private:
    double z0_;
    double big_scale_;
    Cone cone_;
    Parabaloid parabaloid_;
};


#endif //JETSHOW_PARABALOIDCONE_H
