//
// Created by ilya on 6/6/17.
//

#ifndef JETSHOW_VFIELD_H
#define JETSHOW_VFIELD_H

#include <Eigen/Eigen>

using Eigen::Vector3d;


class VField {
public:
    virtual Vector3d v(const Vector3d& point) const = 0;
};

class ConstFlatVField: public VField {
public:
    ConstFlatVField(double gamma);
    Vector3d v(const Vector3d& point) const override;

private:
    double gamma_;
};


class ConstCenralVField: public VField {
public:
    ConstCenralVField(double gamma);
    Vector3d v(const Vector3d& point) const override;

private:
    double gamma_;
};

#endif //JETSHOW_VFIELD_H
