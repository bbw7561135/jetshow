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
		explicit ConstFlatVField(double gamma);
    Vector3d v(const Vector3d& point) const override;

private:
    double gamma_;
};


class ShearedFlatVField: public VField {
public:
		ShearedFlatVField(double gamma_axis, double gamma_border, double r);
		Vector3d v(const Vector3d& point) const override ;

private:
		double gamma_axis_;
		double gamma_border_;
		double r_;
};


class SheathFlatVField: public VField {
public:
		SheathFlatVField(double gamma_spine, double gamma_sheath, double r_sheath);
		Vector3d v(const Vector3d& point) const override ;

private:
		double gamma_spine_;
		double gamma_sheath_;
		double r_sheath_;
};


class ConstCentralVField: public VField {
public:
		explicit ConstCentralVField(double gamma);
    Vector3d v(const Vector3d& point) const override;

private:
    double gamma_;
};


class ShearedCentralVField: public VField {
public:
		ShearedCentralVField(double gamma_axis, double gamma_border, double theta);
		Vector3d v(const Vector3d& point) const override ;

private:
		double gamma_axis_;
		double gamma_border_;
		double theta_;
};


class SheathCentralVField: public VField {
public:
		SheathCentralVField(double gamma_spine, double gamma_sheath,
		                    double theta_sheath);
		Vector3d v(const Vector3d& point) const override ;

private:
		double gamma_spine_;
		double gamma_sheath_;
		double theta_sheath_;
};


class ConstParabolicVField: public VField {
public:
    ConstParabolicVField(double gamma, double R0);
    Vector3d v(const Vector3d& point) const override;

private:
    double gamma_;
    double R0_;
};


// gamma(z) = gamma0*sqrt(z)
class AccParabolicVField: public VField {
public:
    AccParabolicVField(double gamma0, double R0);
    Vector3d v(const Vector3d& point) const override;

private:
    double gamma0_;
    double R0_;
};


class ShearedAccParabolicVField: public VField {
public:
    ShearedAccParabolicVField(double gamma_axis0, double gamma_border0, double R0, double R0_border);
    Vector3d v(const Vector3d& point) const override ;

private:
    double gamma_axis0_;
    double gamma_border0_;
    double R0_;
    double R0_border_;
};


//class ConstParabolicConstConeVField: public VField {
//public:
//    ConstParabolicConstConeVField(double gamma, double R0, double z0);
//    Vector3d v(const Vector3d& point) const override;
//
//private:
//    double z0_;
//    ConstCentralVField conev;
//    ConstParabolicVField parav;
//};


class AccParabolicConstConeVField: public VField {
public:
    AccParabolicConstConeVField(double gamma0, double R0, double z0);
    Vector3d v(const Vector3d& point) const override;

private:
    double z0_;
    ConstCentralVField conev;
    AccParabolicVField parav;
};


class ShearedAccParabolicConstConeVField: public VField {
public:
    ShearedAccParabolicConstConeVField(double gamma_axis0, double gamma_border0, double R0, double R0_border,
            double z0);
    Vector3d v(const Vector3d& point) const override;

private:
    double z0_;
    ShearedCentralVField conev;
    ShearedAccParabolicVField parav;
};


#endif //JETSHOW_VFIELD_H
