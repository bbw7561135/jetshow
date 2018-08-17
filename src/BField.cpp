#include <iostream>
#include <math.h>
#include <Eigen/Eigen>
#include <boost/math/special_functions/bessel.hpp>
#include <BField.h>
#include <utils.h>
#include <Cells.h>
#include <boost/random/variate_generator.hpp>


#include <CGAL/Cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Interpolation_traits_2.h>
#include <CGAL/natural_neighbor_coordinates_2.h>
#include <CGAL/interpolation_functions.h>
#include <CGAL/Barycentric_coordinates_2/Triangle_coordinates_2.h>

typedef CGAL::Cartesian<double>                                   K_;
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
	return std::move(Vector3d(std::move(res.data())));
}


SimulationBField::SimulationBField(Delaunay_triangulation *tr_p, Delaunay_triangulation *tr_fi) {
    tr_p_ = tr_p;
    tr_fi_ = tr_fi;
}

Vector3d SimulationBField::bf(const Vector3d &point) const {
    // Conver 3D point to (r, r_p) coordinates
    double x = point[0];
    double y = point[1];
    double z = point[2];
    double fi = atan(y/x);
    double r_p = hypot(x, y);
    Point_ pt(z, r_p);

    Delaunay_triangulation::Face_handle fh_p = tr_p_->locate(pt);
    Delaunay_triangulation::Face_handle fh_fi = tr_fi_->locate(pt);

    std::vector<Point_ > vertexes_p;
    std::vector<Point_ > vertexes_fi;
    std::vector<double> info_p;
    std::vector<double> info_fi;


    for (int i=0; i<3; i++) {
        vertexes_p.push_back(fh_p->vertex(i)->point());
        info_p.push_back(fh_p->vertex(i)->info());

        std::cout << "Triangle:\t" << tr_p_->triangle(fh_p) << std::endl;
        std::cout << "Vertex 0:\t" << tr_p_->triangle(fh_p)[i] << std::endl;
        std::cout << "B_p data:\t" << fh_p->vertex(i)->info() << std::endl;

        vertexes_fi.push_back(fh_p->vertex(i)->point());
        info_fi.push_back(fh_fi->vertex(i)->info());
        std::cout << "Triangle:\t" << tr_fi_->triangle(fh_fi) << std::endl;
        std::cout << "Vertex 0:\t" << tr_fi_->triangle(fh_fi)[i] << std::endl;
        std::cout << "B_fi data:\t" << fh_fi->vertex(i)->info() << std::endl;
    }

    // Create an std::vector to store coordinates.
    Scalar_vector coordinates_p;
    Scalar_vector coordinates_fi;
    // Instantiate the class Triangle_coordinates_2 for the triangle defined above.
    Triangle_coordinates triangle_coordinates_p(vertexes_p[0], vertexes_p[1], vertexes_p[2]);
    Triangle_coordinates triangle_coordinates_fi(vertexes_fi[0], vertexes_fi[1], vertexes_fi[2]);


    triangle_coordinates_p(pt, std::inserter(coordinates_p, coordinates_p.end()));
    triangle_coordinates_fi(pt, std::inserter(coordinates_fi, coordinates_fi.end()));

    double interpolated_value_p = 0;
    double interpolated_value_fi = 0;
    for(int j = 0; j < 3; ++j) {
        std::cout << "coordinate for B_p" << j + 1 << " = " << coordinates_p[j] << "; ";
        std::cout << "coordinate for B_phi" << j + 1 << " = " << coordinates_fi[j] << "; ";
        interpolated_value_p += coordinates_p[j]*info_p[j];
        interpolated_value_fi += coordinates_fi[j]*info_fi[j];
    }
    std::cout << "Interpolated B_p = " << interpolated_value_p << std::endl;
    std::cout << "Interpolated B_fi = " << interpolated_value_fi << std::endl;
    return Vector3d{-sin(fi)*interpolated_value_fi, cos(fi)*interpolated_value_fi, interpolated_value_p};
}