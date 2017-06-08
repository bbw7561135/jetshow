//
// Created by ilya on 5/31/17.
//

#ifndef JETSHOW_INTERCEPTION_H
#define JETSHOW_INTERCEPTION_H

#include <Eigen/Eigen>

using Eigen::Vector3d;
using std::vector;


class Intersection {
public:
    Intersection(const Vector3d &direction, vector<Vector3d> &points);
    const Vector3d& direction() const;
    vector<Vector3d> points();
    bool has_intersection() const;

private:
    Vector3d direction_;
    vector<Vector3d> points_;

};


#endif //JETSHOW_INTERCEPTION_H
