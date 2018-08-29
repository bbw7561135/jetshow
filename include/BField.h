#ifndef JETSHOW_BFIELDS_H
#define JETSHOW_BFIELDS_H

#include <Eigen/Eigen>
#include "Cells.h"
#include "SimulationInterpolater.h"
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_on_sphere.hpp>
#include <boost/random/variate_generator.hpp>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Interpolation_traits_2.h>
#include <CGAL/natural_neighbor_coordinates_2.h>
#include <CGAL/interpolation_functions.h>
#include <CGAL/Barycentric_coordinates_2/Triangle_coordinates_2.h>

typedef CGAL::Simple_cartesian<double>                                   K_;
typedef K_::Point_2                                                Point_;
typedef CGAL::Triangulation_vertex_base_with_info_2<double, K_>      Vb;
typedef CGAL::Triangulation_data_structure_2<Vb>                  Tds;
typedef CGAL::Delaunay_triangulation_2<K_, Tds>                    Delaunay_triangulation;
typedef K_::FT                                               Coord_type;
typedef std::vector<Coord_type >                            Scalar_vector;
typedef CGAL::Barycentric_coordinates::Triangle_coordinates_2<K_> Triangle_coordinates;


using Eigen::Vector3d;
typedef boost::random::mt19937 gen_type;


// B-field with vector values, e.g. ordered.
class VectorBField {
public:
    virtual Vector3d bf(const Vector3d &point) const = 0 ;
};


// B-field that has no preferred direction, e.g. random.
class ScalarBField {
public:
    virtual double bf(const Vector3d &point) const = 0 ;
};


class RandomScalarBField : public ScalarBField {
public:
		RandomScalarBField(double b_0, double m_b);
		double bf(const Vector3d &point) const override;

private:
		double b_0_;
		double m_b_;
};


// B-Field like ``RandomScalar`` that depends on z-coordinate only
class RandomScalarBFieldZ : public ScalarBField {
public:
		RandomScalarBFieldZ(double b_0, double m_b);
		double bf(const Vector3d &point) const override;

private:
		double b_0_;
		double m_b_;
};


// Random B-field with exponent changing at some distance (e.g. due to jet collimation)
class CompositeRandomScalarBFieldZ : public ScalarBField {
public:
    CompositeRandomScalarBFieldZ(double b_0, double m_b_inner, double m_b_outer, double z0);
    double bf(const Vector3d &point) const override;

private:
    double z0_;
    RandomScalarBFieldZ inner_field_;
    RandomScalarBFieldZ outer_field_;
};


class ConstCylinderBField : public VectorBField {
public:
    ConstCylinderBField(double b_0, double n_b) ;
    Vector3d bf(const Vector3d &point) const override ;
private:
    double b_0_;
    double n_b_;

};


// B-Field like ``ConstCylinder`` that depends on z-coordinate only
class ConstCylinderBFieldZ : public VectorBField {
public:
		ConstCylinderBFieldZ (double b_0, double n_b) ;
		Vector3d bf(const Vector3d &point) const override ;
private:
		double b_0_;
		double n_b_;

};


class RadialConicalBField : public VectorBField {
public:
    RadialConicalBField(double b_0, double n_b) ;
    Vector3d bf(const Vector3d &point) const override ;
private:
    double b_0_;
    double n_b_;
};


class ToroidalBField : public VectorBField {
public:
		ToroidalBField(double b_0, double n_b) ;
		Vector3d bf(const Vector3d &point) const override ;
private:
		double b_0_;
		double n_b_;
};


class HelicalCylinderBField : public VectorBField {
public:
    HelicalCylinderBField(double b_0, double pitch_angle) ;
    Vector3d bf(const Vector3d &point) const override ;
private:
    double b_0_;
    double pitch_angle_;
};


class SpiralConicalBField : public VectorBField {
public:
    SpiralConicalBField(double b_0, double pitch_angle) ;
    Vector3d bf(const Vector3d &point) const override ;
private:
    double b_0_;
    double pitch_angle_;
};


class ForceFreeCylindricalBField : public VectorBField {
public:
    ForceFreeCylindricalBField(double b_0, double mu) ;
    Vector3d bf(const Vector3d &point) const override ;
private:
    double b_0_;
    double mu_;
};


// B-field not completely random, e.g. made of finite cells.
class RandomVectorBField : public VectorBField {
public:
		RandomVectorBField(VectorBField* bfield, double rnd_fraction);
		virtual Vector3d direction(const Vector3d &point) const = 0 ;
		Vector3d bf(const Vector3d &point) const override ;

protected:
		VectorBField* bfield_;
		double rnd_fraction_;
};


class CellsRandomVectorBField : public RandomVectorBField {
public:
		CellsRandomVectorBField(Cells* cells, VectorBField* bfield, double rnd_fraction);
		Vector3d direction(const Vector3d &point) const override ;

private:
		Cells* cells_;
};


class PointsRandomVectorBField : public RandomVectorBField {
public:
		PointsRandomVectorBField(VectorBField* bfield, double rnd_fraction, unsigned int seed);
		Vector3d direction(const Vector3d &point) const override ;

private:
		mutable std::vector<boost::variate_generator<gen_type, boost::uniform_on_sphere<double>>> randoms_on_sphere;
};


// Special class for Elena's calculations of axisymmetric flows.
class SimulationBField : public VectorBField {
public:
    SimulationBField(Delaunay_triangulation *tr_p, Delaunay_triangulation *tr_fi);
    Vector3d bf(const Vector3d &point) const override ;

private:
    SimulationInterpolater interp_p_;
    SimulationInterpolater interp_fi_;
};
#endif //JETSHOW_BFIELDS_H
