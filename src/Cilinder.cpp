//
// Created by ilya on 6/8/17.
//

#include <iostream>
#include "Cilinder.h"

using std::min;
using std::max;


Cilinder::Cilinder(const Vector3d &origin_, const Vector3d &direction_,
                   double r_) : origin_(origin_), direction_(direction_),
                                r_(r_) {}

const Vector3d& Cilinder::origin() const {
    return origin_;
}

const Vector3d& Cilinder::direction() const {
    return direction_;
}

double Cilinder::r() const {
    return r_;
}


Intersection Cilinder::hit(Ray &ray) const {
    const Vector3d ray_origin = ray.origin();
    const Vector3d ray_direction = ray.direction();
    const Vector3d &cilinder_origin = origin();
    const Vector3d &cilinder_direction = direction();

    const Vector3d delta = ray_origin - cilinder_origin;

    const Vector3d v = ray_direction -
                 ray_direction.dot(cilinder_direction) * cilinder_direction;
    double c2 = v.squaredNorm();
    const Vector3d w = delta - delta.dot(cilinder_direction) * cilinder_direction;
    double c1 = v.dot(w);
    double c0 = w.squaredNorm() - r();

    std::cout << "c2 = " << c2 << " c1 = " << c1 << " c0 = " << c0 << std::endl;

    double d = c1*c1 - c0*c2;
    if (d < 0) {
        std::cout << "No intercestions, d = " << d << std::endl;
        vector<Vector3d> points;
        return Intersection(ray_direction, points);
    }
    else if (d > 0) {
        double eps = (c1 > 0 ? 1 : -1);
        double sqrt_d = sqrt(d);
        double t1 = (-c1 - eps*sqrt_d)/(c2);
        double t2 = c0/(-c1 - eps*sqrt_d);
        std::cout << "t1: " << t1 << std::endl;
        std::cout << "t2: " << t2 << std::endl;
        Vector3d point_in = ray.point(min(t1, t2));
        Vector3d point_out = ray.point(max(t1, t2));
        vector<Vector3d> points = {point_in, point_out};
        return Intersection(ray_direction, points);
    } else {
        if (c2 == 1 && c1 ==0 && c0 == 0) {
            std::cout << "One intersection" << std::endl;
            double t = -c1/c2;
            Vector3d point = ray.point(t);
            vector<Vector3d> points = {point};
            return Intersection(ray_direction, points);
        }
        else if (c2 ==0 && c1 == 0 && c0 == 0) {
            std::cout << "Along border" << std::endl;
            vector<Vector3d> points;
            return Intersection(ray_direction, points);
        }
        else if (c2 == 0 && c1 == 0 && c0 > 0) {
            std::cout << "No interception. Externally along border" << std::endl;
            vector<Vector3d> points;
            return Intersection(ray_direction, points);
        }
        else if (c2 == 0 && c1 == 0 && c0 <= 0) {
            std::cout << "No interception. Internally along border" << std::endl;
            vector<Vector3d> points;
            return Intersection(ray_direction, points);
        }
    }
}
