#ifndef JETSHOW_BFIELDS_H
#define JETSHOW_BFIELDS_H

#include <Eigen/Eigen>
#include "Cells.h"
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_on_sphere.hpp>
#include <boost/random/variate_generator.hpp>


using Eigen::Vector3d;
typedef boost::random::mt19937 gen_type;


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
    HelicalCylinderBField(double b_0, double pitch_angle) ;
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


class RandomBField : public BField {
public:
		RandomBField(BField* bfield, double rnd_fraction);
		virtual Vector3d direction(const Vector3d &point) const = 0 ;
		Vector3d bf(const Vector3d &point) const override ;

protected:
		BField* bfield_;
		double rnd_fraction_;
};


class RandomCellsBField : public RandomBField {
public:
		RandomCellsBField(Cells* cells, BField* bfield, double rnd_fraction);
		Vector3d direction(const Vector3d &point) const override ;

private:
		Cells* cells_;
};


class RandomPointBField : public RandomBField {
public:
		RandomPointBField(BField* bfield, double rnd_fraction, unsigned int seed);
		Vector3d direction(const Vector3d &point) const override ;

private:
		mutable std::vector<boost::variate_generator<gen_type, boost::uniform_on_sphere<double>>> randoms_on_sphere;
};

#endif //JETSHOW_BFIELDS_H
