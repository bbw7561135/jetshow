//
// Created by ilya on 13.08.18.
//

#include <math.h>
#include <utils.h>
#include "Parabaloid.h"
#include "Geometry.h"


Parabaloid::Parabaloid(const Vector3d &origin, const Vector3d &direction, const double &r0) {
    origin_ = origin;
    direction_ = direction;
    direction_.normalize();
    r0_ = r0;
}

const Vector3d& Parabaloid::origin() const {
    return origin_;
}

const Vector3d& Parabaloid::direction() const {
    return direction_;
}

double Parabaloid::r0() const {
    return r0_;
}

bool Parabaloid::is_within(Vector3d& point) const {
    return std::hypot(point[0], point[1]) < r0_*sqrt(point[2]);
}

std::list<Intersection> Parabaloid::hit(Ray &ray) const {
    std::list<double> ts = intersection(ray.origin(), ray.direction(), 1., 1., 0., 0., 0., 0., 0., 0., -r0_*r0_);

    // No intersection - most common case
    if (ts.empty()) {
        return std::list<Intersection>{};

    }
    // Two intersections
    else if (ts.size() == 2) {
        double t1 = ts.front();
        double t2 = ts.back();
        Vector3d point_in = ray.point(pc*std::min(t1, t2));
        Vector3d point_out = ray.point(pc*std::max(t1, t2));
        return std::list<Intersection>{Intersection(ray, point_in, point_out)};
    }
    // One intersection - ignore it for now
    else {
        return std::list<Intersection>{};
    }
}