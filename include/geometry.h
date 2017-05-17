//
// Created by ilya on 12/15/16.
//

#ifndef JETSHOW_GEOMETRY_H
#define JETSHOW_GEOMETRY_H

#include <Eigen/Eigen>

using Eigen::Vector3d;
using std::pair;

class Ray
{
public:
    Vector3d point(double t);

private:
    Vector3d origin;
    Vector3d direction;
};


class Geometry
{
public:
    virtual pair<double, double>  hit(Ray &ray) = 0;

};


class Cone: public Geometry
{
public:

    Cone(Vector3d &origin, Vector3d &direction, double &angle);

private:
    Vector3d origin;
    Vector3d direction;
    double angle;
};


/*class Intersection
{
public:
    Intersection(Geometry &geometry, Ray &ray);

private:
    Ray ray;


};*/
#endif //JETSHOW_GEOMETRY_H
