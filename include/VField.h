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
    ConstCentralVField(double gamma);
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

#endif //JETSHOW_VFIELD_H
