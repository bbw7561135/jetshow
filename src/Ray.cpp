//
// Created by ilya on 5/31/17.
//

#include <iostream>
#include "Ray.h"


Ray::Ray(Vector3d &neworigin, Vector3d &newdirection) {
    origin_ = neworigin;
    direction_ = newdirection;
    direction_.normalize();
//    std::cout << "Initialized with origin:" << std::endl;
//    std::cout << origin_ << std::endl;
}


Vector3d Ray::point(double t) {
    return origin_ + direction_*t;
}

Vector3d Ray::direction() {
    return direction_;
}

Vector3d Ray::origin() {
    return origin_;
}