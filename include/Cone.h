//
// Created by ilya on 5/31/17.
//

#ifndef JETSHOW_CONE_H
#define JETSHOW_CONE_H


#include "Geometry.h"


class Cone: public Geometry {
public:
    Cone(const Vector3d &origin, const Vector3d &direction, const double &angle,
         const double &scale);
    const Vector3d& origin() const ;
    const Vector3d& direction() const ;
    const double angle() const;
    std::list<Intersection> hit(Ray &ray) const override;
    bool is_within(Vector3d& point) const override;
    double const big_scale() const override ;

private:
    Vector3d origin_;
    Vector3d direction_;
    double angle_;
    double cos_angle_;
    double big_scale_;
};

void Print(const std::vector<Vector3d>& v);


#endif //JETSHOW_CONE_H
