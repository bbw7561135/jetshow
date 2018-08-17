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
		explicit ConstCentralVField(double gamma, Vector3d origin={0, 0, 0});
    Vector3d v(const Vector3d& point) const override;

private:
    double gamma_;
    Vector3d origin_;
};


class ShearedCentralVField: public VField {
public:
		ShearedCentralVField(double gamma_axis, double gamma_border, double theta,
                             Vector3d origin={0, 0, 0});
		Vector3d v(const Vector3d& point) const override ;

private:
		double gamma_axis_;
		double gamma_border_;
		double theta_;
        Vector3d origin_;

};


class SheathCentralVField: public VField {
public:
		SheathCentralVField(double gamma_spine, double gamma_sheath,
		                    double theta_sheath, Vector3d origin={0, 0, 0});
		Vector3d v(const Vector3d& point) const override ;

private:
		double gamma_spine_;
		double gamma_sheath_;
		double theta_sheath_;
		Vector3d origin_;
};

// Rz0 - radius of parabaloid at z=1pc
class ConstParabolicVField: public VField {
public:
    ConstParabolicVField(double gamma, double Rz0);
    Vector3d v(const Vector3d& point) const override;

private:
    double gamma_;
    double Rz0_;
};


// gamma(z) = 1+a*sqrt(z)
// gamma0 - speed at z=R0
// Rz0 - radius of parabaloid at z=1pc
class AccParabolicVField: public VField {
public:
    AccParabolicVField(double gamma0, double R0, double Rz0);
    Vector3d v(const Vector3d& point) const override;

private:
    double gamma0_;
    double R0_;
    double Rz0_;
};

// gamma(z) = 1+a*sqrt(z)
// gamma_axis0 - speed at z=R0 at radius=0
// gamma_border0 - speed at z=R0 at the border of paraboloid
// Rz0 - radius of parabaloid at z=1pc
class ShearedAccParabolicVField: public VField {
public:
    ShearedAccParabolicVField(double gamma_axis0, double gamma_border0, double R0, double Rz0);
    Vector3d v(const Vector3d& point) const override ;

private:
    double gamma_axis0_;
    double gamma_border0_;
    double R0_;
    double Rz0_;
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


// gamma(z) = 1+a*sqrt(z)
// gamma0 - maximum speed (at z=z0)
// Rz0 - radius of parabaloid at z=1pc
// z0 - z-coordinate where paraboloid goes into cone
class AccParabolicConstConeVField: public VField {
public:
    AccParabolicConstConeVField(double gamma0, double Rz0, double z0);
    Vector3d v(const Vector3d& point) const override;

private:
    double z0_;
    ConstCentralVField conev;
    AccParabolicVField parav;
};


// gamma(z) = 1+a*sqrt(z)
// gamma_axis0 - speed at z=z0 at radius=0
// gamma_border0 - speed at z=z0 at the border of paraboloid
// Rz0 - radius of parabaloid at z=1pc
// z0 - z-coordinate where paraboloid goes into cone
class ShearedAccParabolicConstConeVField: public VField {
public:
    ShearedAccParabolicConstConeVField(double gamma_axis0, double gamma_border0, double Rz0, double z0);
    Vector3d v(const Vector3d& point) const override;

private:
    double z0_;
    ShearedCentralVField conev;
    ShearedAccParabolicVField parav;
};


#endif //JETSHOW_VFIELD_H
