#ifndef JETSHOW_BFIELDS_H
#define JETSHOW_BFIELDS_H

#include <Eigen/Eigen>

using Eigen::Vector3d;


class BField {
public:
    virtual Vector3d bf(const Vector3d &point) const = 0 ;
};


class ConstCylinderBField : public BField {
public:
    ConstCylinderBField(double b_0, double n_b) ;
    Vector3d bf(const Vector3d &point) const override ;
private:
    double b_0_;
    double n_b_;

};


class RadialConicalBField : public BField {
public:
    RadialConicalBField(double b_0, double n_b) ;
    Vector3d bf(const Vector3d &point) const override ;
private:
    double b_0_;
    double n_b_;
};


class HelicalCylinderBField : public BField {
public:
    HelicalCylinderField(double b_0, double pitch_angle) ;
    Vector3d bf(const Vector3d &point) const override ;
private:
    double b_0_;
    double pitch_angle_;
};


class SpiralConicalBField : public BField {
public:
    SpiralConicalBField(double b_0, double pitch_angle) ;
    Vector3d bf(const Vector3d &point) const override ;
private:
    double b_0_;
    double pitch_angle_;
};


class ForceFreeCylindricalBField : public BField {
public:
    ForceFreeCylindricalBField(double b_0, double mu) ;
    Vector3d bf(const Vector3d &point) const override ;
private:
    double b_0_;
    double mu_;
};

// This is old code for helical field. Note that it is likely to be wrong.
//class BField
//{
//public:
//    BField(double z0, double fiValue0, double zValue0);
//
//    Vector3d bf(Vector3d &p);
//    double fiValue(Vector3d &p);
//    double zValue(Vector3d &p);
//    double getZ0() { return z0;}
//    double getFiValue0() { return fiValue0;}
//    double getZValue0() { return zValue0;}
//
//
//private:
//    double z0;
//    double fiValue0;
//    double zValue0;
//};

#endif //JETSHOW_BFIELDS_H
