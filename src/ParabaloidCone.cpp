//
// Created by ilya on 13.08.18.
//

#include <math.h>
#include <utils.h>
#include "Geometry.h"
#include "ParabaloidCone.h"
#include "Cone.h"


ParabaloidCone::ParabaloidCone(const Vector3d &origin, const Vector3d &direction, const double &r0, const double &z0,
        const double &big_scale) :
        cone_(Vector3d(0, 0, -z0*(1/(0.5*r0*pow(z0, -0.5)) - 1)), direction, atan(0.5*r0*pow(z0, -0.5)), big_scale),
        parabaloid_(origin, direction, r0, big_scale)  {
    z0_ = z0;
}

const Vector3d& ParabaloidCone::origin() const {
    return parabaloid_.origin();
}

const Vector3d& ParabaloidCone::direction() const {
    return cone_.direction();
}

double ParabaloidCone::z0() const {
    return z0_;
}

const double ParabaloidCone::big_scale() const {
    return cone_.big_scale();
}

bool ParabaloidCone::is_within(Vector3d& point) const {
    if (point[2] < z0_) {
        return parabaloid_.is_within(point);
    }
    else {
        return cone_.is_within(point);
    }
}

std::list<Intersection> ParabaloidCone::hit(Ray &ray) const {
    auto result_cone = cone_.hit(ray);
    auto result_para = parabaloid_.hit(ray);
    Vector3d point_in;
    Vector3d point_out;

    if (result_cone.empty() and result_para.empty()) {
        return std::list<Intersection>{};
    }

    // If there's intersection with Cone but not with Para than it is case
    // when both ``point_in.z`` and ``point_out.`` > or < ``z0``
    if (!result_cone.empty() and result_para.empty()) {
        auto points = result_cone.front().get_path();
        Vector3d point_in = points.first;
        Vector3d point_out = points.second;

        if (point_in[2] > z0_) {
            return result_cone;
        }
        else {
            return std::list<Intersection>{};
        }
    }

    // If there's intersection with Para than always intersect Cone
    if (!result_para.empty()) {
        auto points_cone = result_cone.front().get_path();
        Vector3d point_in_cone = points_cone.first;
        Vector3d point_out_cone = points_cone.second;
        auto points_para = result_para.front().get_path();
        Vector3d point_in_para = points_para.first;
        Vector3d point_out_para = points_para.second;

        if (point_in_para[2] > z0_) {
            point_in = point_in_cone;
        }
        else {
            point_in = point_in_para;
        }

        if (point_out_para[2] < z0_) {
            point_out = point_out_para;
        }
        else {
            point_out = point_out_cone;
        }
        return std::list<Intersection>{Intersection(ray, point_in, point_out)};
    }
}