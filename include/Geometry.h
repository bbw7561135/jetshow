//
// Created by ilya on 5/31/17.
//

#ifndef JETSHOW_GEOMETRY_H
#define JETSHOW_GEOMETRY_H

#include <Eigen/Eigen>
#include "Ray.h"
#include "Intersection.h"

using Eigen::Vector3d;

class Ray;
class Intersection;


class Geometry {
    public:
        virtual const Vector3d& origin() const = 0;
        virtual const double big_scale() const = 0;
        virtual std::list<Intersection> hit(Ray &ray) const = 0;
        virtual bool is_within(Vector3d& point) const = 0;
        std::pair<Vector3d, Vector3d> full_infinite_path(Ray &ray) const;
        std::pair<Vector3d, Vector3d> half_infinite_path(Ray &ray,
                                                         const Vector3d &point) const;

};


#endif //JETSHOW_GEOMETRY_H
