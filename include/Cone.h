//
// Created by ilya on 5/31/17.
//

#ifndef JETSHOW_CONE_H
#define JETSHOW_CONE_H


#include "Geometry.h"


class Cone: public Geometry {
public:
    Cone(Vector3d &origin, Vector3d &direction, double &angle);
    Vector3d origin();
    Vector3d direction();
    double angle();
    Intersection hit(Ray &ray);

private:
    Vector3d origin_;
    Vector3d direction_;
    double angle_;
    Intersection getIntersection(double a, double b, double c, double d, Ray &ray);


};

void Print(const std::vector<Vector3d>& v);


void Print(const Vector3d& v);

#endif //JETSHOW_CONE_H
