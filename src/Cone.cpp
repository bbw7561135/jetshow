//
// Created by ilya on 5/31/17.
//

#include <math.h>
#include <iostream>
#include "Cone.h"
#include <vector>

using std::min;
using std::max;


Cone::Cone(Vector3d &neworigin, Vector3d &newdirection, double &newangle) {
    origin_ = neworigin;
    direction_ = newdirection;
    direction_.normalize();
    angle_ = newangle;

}

Vector3d Cone::origin() {
    return origin_;
}

Vector3d Cone::direction() {
    return direction_;
}

double Cone::angle() {
    return angle_;
}


Intersection Cone::hit(Ray &ray) {
//  (0,0,1)-(0,0,0)
//    std::cout << "Ray origin:" << std::endl;
//    std::cout << ray.origin() << std::endl;
//    std::cout << "Cone origin:" << std::endl;
//    std::cout << origin() << std::endl;
    const Vector3d dp = ray.origin() - origin();
//    std::cout << "dp:" << std::endl;
//    std::cout << dp << std::endl;
    Vector3d expr1 = ray.direction() - ray.direction().dot(direction()) * direction();
    Vector3d expr2 = dp - dp.dot(direction()) * direction();
    double a = pow(cos(angle()),2)*expr1.dot(expr2) - pow(sin(angle()),2)*pow(ray.direction().dot(direction()),2);
    double b = 2.*pow(cos(angle()),2)*expr1.dot(expr2) - 2.*pow(sin(angle()),2)*ray.direction().dot(direction())*dp.dot(direction());
    // TODO: Check line below - ``expr2.dot(expr2)``
    double c = pow(cos(angle()),2)*expr2.dot(expr2) - pow(sin(angle()),2)*pow(dp.dot(direction()),2);
    double d = b*b - 4.*a*c;
    std::cout << "a: " << a << std::endl;
    std::cout << "b: " << b << std::endl;
    std::cout << "c: " << c << std::endl;
    std::cout << "d: " << d << std::endl;

    return getIntersection(a, b, c, d, ray);
}

Intersection Cone::getIntersection(double a, double b, double c, double d, Ray &ray) {
    Vector3d ray_direction = ray.direction();
    // No intersection case
    if (d < 0.) {
        std::cout << "No intercestions" << std::endl;
        vector<Vector3d> points;
        return Intersection(ray_direction, points);

    }

        // Ray goes along cone border
    else if (a == 0.) {
        std::cout << "Along border" << std::endl;
        vector<Vector3d> points;
        return Intersection(ray_direction, points);
    } else {
        double t1 = (-b + sqrt(d)) / (2. * a);
        double t2 = (-b - sqrt(d)) / (2. * a);
        std::cout << "t1: " << t1 << std::endl;
        std::cout << "t2: " << t2 << std::endl;
        Vector3d point_in = ray.point(min(t1, t2));
        Vector3d point_out = ray.point(max(t1, t2));
        vector<Vector3d> points = {point_in, point_out};
        return Intersection(ray_direction, points);
    }
}

void Print(const std::vector<Vector3d>& v) {
    for (unsigned i = 0; i < v.size(); i++) {
        std::cout << v[i] << std::endl;
    }
}