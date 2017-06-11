#ifndef JETSHOW_NFIELD_H
#define JETSHOW_NFIELD_H

#include <Eigen/Eigen>

using Eigen::Vector3d;

class NField {
public:
    virtual double n(const Vector3d& point) const = 0;
};

class BKNField: public NField {

    public:
    BKNField(double z_0, double n_0);

    double n(const Vector3d &point) const override;

private:
    double z_0_;
    double n_0_;

};


#endif //JETSHOW_NFIELD_H
