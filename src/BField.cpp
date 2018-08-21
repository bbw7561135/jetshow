#include <iostream>
#include <math.h>
#include <Eigen/Eigen>
#include <boost/math/special_functions/bessel.hpp>
#include <BField.h>
#include <utils.h>
#include <Cells.h>
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
typedef boost::variate_generator<gen_type&, boost::uniform_on_sphere<double>> gg;


RandomScalarBField::RandomScalarBField(double b_0, double m_b) :
		b_0_(b_0), m_b_(m_b) {}

double RandomScalarBField::bf(const Vector3d &point) const {
	double r = point.norm();
	return b_0_ * pow(r/pc, -m_b_);
};


RandomScalarBFieldZ::RandomScalarBFieldZ(double b_0, double m_b) :
		b_0_(b_0), m_b_(m_b) {}

double RandomScalarBFieldZ::bf(const Vector3d &point) const {
	double r = abs(point[2]);
	return b_0_ * pow(r/pc, -m_b_);
};


CompositeRandomScalarBFieldZ::CompositeRandomScalarBFieldZ(double b_0, double m_b_inner, double m_b_outer, double z0) :
        z0_(z0), inner_field_(b_0, m_b_inner), outer_field_(b_0, m_b_outer) {}

double CompositeRandomScalarBFieldZ::bf(const Vector3d &point) const {
    double r = abs(point[2]);
    if (r > z0_) {
        return outer_field_.bf(point);
    }
    else {
        return inner_field_.bf(point);
    }
};


ConstCylinderBField::ConstCylinderBField(double b_0, double n_b) : b_0_(b_0),
                                                                   n_b_(n_b) {};

Vector3d ConstCylinderBField::bf(const Vector3d &point) const {
	double r = point.norm();
	return Vector3d(0.0, 0.0, b_0_*pow(r/pc, -n_b_));
}


ConstCylinderBFieldZ::ConstCylinderBFieldZ(double b_0, double n_b) : b_0_(b_0),
                                                                   n_b_(n_b) {};

Vector3d ConstCylinderBFieldZ::bf(const Vector3d &point) const {
	double r = abs(point[2]);
	return Vector3d(0.0, 0.0, b_0_*pow(r/pc, -n_b_));
}


RadialConicalBField::RadialConicalBField(double b_0, double n_b) : b_0_(b_0),
                                                                   n_b_(n_b) {};

Vector3d RadialConicalBField::bf(const Vector3d &point) const {
    double r = point.norm();
    return Vector3d(b_0_*pow(r/pc, -n_b_)*point[0]/r,
                    b_0_*pow(r/pc, -n_b_)*point[1]/r,
                    b_0_*pow(r/pc, -n_b_)*point[2]/r);
}



ToroidalBField::ToroidalBField(double b_0, double n_b) : b_0_(b_0),
                                                         n_b_(n_b){};

Vector3d ToroidalBField::bf(const Vector3d &point) const {
	double x = point[0];
	double y = point[1];
	double z = point[2];
	double fi = atan(y/x);
	double b = b_0_*pow(z/pc, -n_b_);
	return Vector3d(-sin(fi)*b, cos(fi)*b, 0);
}


HelicalCylinderBField::HelicalCylinderBField(double b_0, double pitch_angle) :
        b_0_(b_0), pitch_angle_(pitch_angle) {};

Vector3d HelicalCylinderBField::bf(const Vector3d &point) const {
    double r = sqrt(point[0]*point[0]+ point[1]*point[1]);
    return Vector3d(b_0_*tan(pitch_angle_*point[1]/r),
                    -b_0_*tan(pitch_angle_*point[0]/r),
                    b_0_);
}


SpiralConicalBField::SpiralConicalBField(double b_0, double pitch_angle) :
        b_0_(b_0), pitch_angle_(pitch_angle) {};

Vector3d SpiralConicalBField::bf(const Vector3d &point) const {
    double z = point[2];
    double x = point[0];
    double y = point[1];
    double b_z = b_0_/(z*z/(pc*pc));
    return Vector3d(b_z*(x/z + y*tan(pitch_angle_)/pc),
                    b_z*(y/z - x*tan(pitch_angle_)/pc),
                    b_z);
}


ForceFreeCylindricalBField::ForceFreeCylindricalBField(double b_0, double mu) :
        b_0_(b_0), mu_(mu) {};

Vector3d ForceFreeCylindricalBField::bf(const Vector3d &point) const {
    double x = point[0];
    double y = point[1];
    double atan_term = atan(y/x);
    double bessel_0 = boost::math::cyl_bessel_i(0, mu_);
    double bessel_1 = boost::math::cyl_bessel_i(1, mu_);
    return Vector3d(-b_0_*bessel_1*sin(atan_term),
                    b_0_*bessel_1*sin(atan_term),
                    b_0_*bessel_0);
}


RandomBField::RandomBField(BField *bfield, double rnd_fraction) :
		rnd_fraction_(rnd_fraction)
{
	bfield_ = bfield;
};

Vector3d RandomBField::bf(const Vector3d &point) const {
	Vector3d b = bfield_->bf(point);
	double bv = b.norm();
	Vector3d n = direction(point);
//	std::cout << "Rnd direction = " << n << std::endl;
	b = b + rnd_fraction_*bv*n;
	return b;
}


RandomCellsBField::RandomCellsBField(Cells* cells, BField* bfield,
                                     double rnd_fraction) :
		RandomBField(bfield, rnd_fraction) {
	cells_ = cells;
}

Vector3d RandomCellsBField::direction(const Vector3d &point) const {
	return cells_->getID(point);
}


RandomPointBField::RandomPointBField(BField* bfield, double rnd_fraction,
                                     unsigned int seed) :
		RandomBField(bfield, rnd_fraction),
		randoms_on_sphere() {
	for (int j = 0; j < omp_get_max_threads(); ++j) {
		gen_type rand_gen;
		rand_gen.seed(j+seed);
		boost::uniform_on_sphere<double> unif_sphere(3);
		boost::variate_generator<gen_type, boost::uniform_on_sphere<double>> random_on_sphere(rand_gen, unif_sphere);
		randoms_on_sphere.push_back(random_on_sphere);
	}

};


Vector3d RandomPointBField::direction(const Vector3d &point) const {
	std::vector<double> res = randoms_on_sphere[omp_get_thread_num()]();
	return std::move(Vector3d(res.data()));
}


SimulationBField::SimulationBField(Delaunay_triangulation *tr_p,
        Delaunay_triangulation *tr_fi) : interp_p_(tr_p), interp_fi_(tr_fi) {}

Vector3d SimulationBField::bf(const Vector3d &point) const {
    double x = point[0];
    double y = point[1];
    double fi = atan(y/x);

    double interpolated_value_p = interp_p_.interpolated_value(point);
    double interpolated_value_fi = interp_fi_.interpolated_value(point);
    return Vector3d{-sin(fi)*interpolated_value_fi, cos(fi)*interpolated_value_fi, interpolated_value_p};
}