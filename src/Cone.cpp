//
// Created by ilya on 5/31/17.
//

#include <math.h>
#include <iostream>
#include "Cone.h"
#include <vector>

using std::min;
using std::max;
using Eigen::Matrix3d;


Cone::Cone(const Vector3d &neworigin, const Vector3d &newdirection,
           const double &newangle, const double &newscale) {
    origin_ = neworigin;
    direction_ = newdirection;
    direction_.normalize();
    angle_ = newangle;
    cos_angle_ = cos(angle_);

}


const Vector3d& Cone::origin() const {
    return origin_;
}


const Vector3d& Cone::direction() const {
    return direction_;
}


const double Cone::angle() const {
    return angle_;
}


const double Cone::big_scale() const {
    return big_scale_;
}


bool Cone::is_within(Vector3d& point) const {
    Vector3d diff = point-origin_;
    diff.normalize();
    double cosp = direction_.dot(diff);
    return cosp > cos_angle_ || cosp < -cos_angle_;
}


std::list<Intersection> Cone::hit(Ray &ray) const {
    const Vector3d& ray_origin = ray.origin();
    const Vector3d& ray_direction = ray.direction();
    const Vector3d& cone_origin = origin();
    const Vector3d& cone_direction = direction();
    Matrix3d eye_matrix;
    eye_matrix << 1, 0, 0,
                  0, 1, 0,
                  0, 0, 1;
    double cone_angle = angle();
    // DP
    const Vector3d delta = ray_origin - cone_origin;
    // M
    Matrix3d M = cone_direction * cone_direction.transpose() - cos(cone_angle)*cos(cone_angle)*eye_matrix;
    double c2 = ray_direction.transpose() * M * ray_direction;
    double c1 = ray_direction.transpose() * M * delta;
    double c0 = delta.transpose() * M * delta;
    std::cout << "c2 = " << c2 << " c1 = " << c1 << " c0 = " << c0 << std::endl;

    if (c2 == 0.) {
        std::cout << "Along border" << std::endl;

        return std::list<Intersection>{Intersection(ray, *this)};
    }

    double d = c1*c1 - c0*c2;
    if (d < 0) {
        std::cout << "No intercestions, d = " << d << std::endl;
        return std::list<Intersection>{};
    }
    else if (d == 0) {
        // One intersection if ray goes through apex of cone or it parallels to
        // any of it's generating lines.
        double t = -c1/c2;
        Vector3d point = ray.point(t);
        if (point == origin_) {
            std::cout << "One intersection at apex" << std::endl;
            return std::list<Intersection>{Intersection(ray, point, point)};
        } else {
            std::cout << "One intersection and infinite path before/past" << std::endl;
            return std::list<Intersection>{Intersection(ray, point, *this)};
        }
    } else {
        // Infinite intersections only if ray parallel to ray direction.
        double eps = (c1 > 0 ? 1 : -1);
        double t1 = (-c1 - eps*sqrt(c1*c1-c2*c0))/(c2);
        double t2 = c0/(-c1 - eps*sqrt(c1*c1-c2*c0));
        std::cout << "t1: " << t1 << std::endl;
        std::cout << "t2: " << t2 << std::endl;
        Vector3d point_in = ray.point(min(t1, t2));
        Vector3d point_out = ray.point(max(t1, t2));
        double cos_ray_cone = ray_direction.dot(cone_direction);
        // TODO: Case of ray parallel to cone direction is very rare - can
        // safely ignore this case in cse of poor performance.
        if (abs(cos_ray_cone) != 1.) {
            std::cout << "Two intersection - finite case" << std::endl;
            return std::list<Intersection>{Intersection(ray, point_in,
                                                        point_out)};
        } else {
            std::cout << "Two intersection - two half-infinite cases" << std::endl;
            return std::list<Intersection>{Intersection(ray, point_out, *this),
                                           Intersection(ray, point_in, *this)};
        }
    }

}

void Print(const std::vector<Vector3d>& v) {
    for (unsigned i = 0; i < v.size(); i++) {
        std::cout << v[i] << std::endl;
    }
}
