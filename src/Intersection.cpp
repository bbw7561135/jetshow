//
// Created by ilya on 5/31/17.
//

#include "Intersection.h"

Intersection::Intersection(const Vector3d& newdirection, vector<Vector3d> &newpoints) {
    direction_ = newdirection;
    points_ = newpoints;
}

const Vector3d& Intersection::direction() const {
    return direction_;
}

vector<Vector3d> Intersection::points() {
    return points_;
}

bool Intersection::has_intersection() const {
    return !points_.empty();
}