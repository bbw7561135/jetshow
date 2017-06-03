#ifndef JETSHOW_BFIELDS_H
#define JETSHOW_BFIELDS_H

#include <Eigen/Eigen>

using Eigen::Vector3d;

class BField
{
public:
    BField(double z0, double fiValue0, double zValue0);

    Vector3d bf(Vector3d &p);
    double fiValue(Vector3d &p);
    double zValue(Vector3d &p);
    double getZ0() { return z0;}
    double getFiValue0() { return fiValue0;}
    double getZValue0() { return zValue0;}


private:
    double z0;
    double fiValue0;
    double zValue0;
};

#endif //JETSHOW_BFIELDS_H
