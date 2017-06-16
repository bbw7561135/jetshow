//
// Created by ilya on 6/6/17.
//

#ifndef JETSHOW_VFIELD_H
#define JETSHOW_VFIELD_H

#include <Eigen/Eigen>

using Eigen::Vector3d;

// TODO: What is more efficient?
// 1. Returning type Vector3d. Thus returning a value. Copying can be avoided
// via RVO?
// 2. Returning a reference. I need lvalue expression inside member function to
// reference it. Thus creating and possibly copying. Is it optimized?
// 3. Using move semantic?
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


class ConstCentralVField: public VField {
public:
    ConstCentralVField(double gamma);
    Vector3d v(const Vector3d& point) const override;

private:
    double gamma_;
};

#endif //JETSHOW_VFIELD_H
