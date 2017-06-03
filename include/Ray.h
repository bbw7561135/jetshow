//
// Created by ilya on 5/31/17.
//

#ifndef JETSHOW_RAY_H
#define JETSHOW_RAY_H

#include <Eigen/Eigen>
#include "Geometry.h"


using Eigen::Vector3d;


class Ray {
    friend class Cone;
public:
        Ray(Vector3d &origin, Vector3d &direction);
        Vector3d point(double t);
        Vector3d origin();
        Vector3d direction();

protected:
        Vector3d origin_;
        Vector3d direction_;
};


#endif //JETSHOW_RAY_H
