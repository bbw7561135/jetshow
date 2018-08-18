#include <math.h>
#include <Eigen/Eigen>
#include <utils.h>
#include "VField.h"


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

using Eigen::Vector3d;


ConstFlatVField::ConstFlatVField(double gamma) : gamma_(gamma) {};

Vector3d ConstFlatVField::v(const Vector3d &point) const {
    return Vector3d(0, 0, c*sqrt(1. - 1./(gamma_*gamma_)));
}


ShearedFlatVField::ShearedFlatVField(double gamma_axis, double gamma_border,
                                     double r) : gamma_axis_(gamma_axis),
                                                 gamma_border_(gamma_border),
                                                 r_(r) {}

Vector3d ShearedFlatVField::v(const Vector3d &point) const {
	double x = point[0];
	double y = point[1];
	double r = sqrt(x*x + y*y);
	double gamma = gamma_axis_-(gamma_axis_-gamma_border_)*r/r_;
	return Vector3d(0, 0, c*sqrt(1. - 1./(gamma*gamma)));
}


SheathFlatVField::SheathFlatVField(double gamma_spine, double gamma_sheath,
                                   double r_sheath) :
		gamma_spine_(gamma_spine), gamma_sheath_(gamma_sheath),
		r_sheath_(r_sheath) {}

Vector3d SheathFlatVField::v(const Vector3d &point) const {
	double x = point[0];
	double y = point[1];
	double r = sqrt(x*x + y*y);
	double gamma;
	if (r < r_sheath_) {
		gamma = gamma_spine_;
	} else {
		gamma = gamma_sheath_;
	}
	return Vector3d(0, 0, c*sqrt(1. - 1./(gamma*gamma)));
}


ConstCentralVField::ConstCentralVField(double gamma, Vector3d origin) :
gamma_(gamma), origin_(origin) {}

Vector3d ConstCentralVField::v(const Vector3d &point) const {
    double x = point[0];
    double y = point[1];
    double z = point[2];
    double r = (point + origin_).norm();
    double v_r = c*sqrt(1. - 1./(gamma_*gamma_));
    return Vector3d(v_r*(x+origin_[0])/r, v_r*(y+origin_[1])/r, v_r*(z+origin_[2])/r);

};


ShearedCentralVField::ShearedCentralVField(double gamma_axis,
                                           double gamma_border, double theta,
                                           Vector3d origin) :
		gamma_axis_(gamma_axis), gamma_border_(gamma_border), theta_(theta),
		origin_(origin) {}

Vector3d ShearedCentralVField::v(const Vector3d &point) const {
    double x = point[0];
    double y = point[1];
    double z = point[2];
    double r = (point + origin_).norm();
	double theta = acos(z + origin_[2]/r);
	double gamma = gamma_axis_-(gamma_axis_-gamma_border_)*theta/theta_;
	double v_r = c*sqrt(1. - 1./(gamma*gamma));
	return Vector3d(v_r*(x+origin_[0])/r, v_r*(y+origin_[1])/r, v_r*(z+origin_[2])/r);
}


SheathCentralVField::SheathCentralVField(double gamma_spine,
                                         double gamma_sheath,
                                         double theta_sheath,
                                         Vector3d origin) :
		gamma_spine_(gamma_spine), gamma_sheath_(gamma_sheath),
		theta_sheath_(theta_sheath), origin_(origin) {}

Vector3d SheathCentralVField::v(const Vector3d &point) const {
    double x = point[0];
    double y = point[1];
    double z = point[2];
    double r = (point + origin_).norm();
	double theta = acos(z/r);
	double gamma;
	if (theta < theta_sheath_) {
		gamma = gamma_spine_;
	} else {
		gamma = gamma_sheath_;
	}
	double v_r = c*sqrt(1. - 1./(gamma*gamma));
    return Vector3d(v_r*(x+origin_[0])/r, v_r*(y+origin_[1])/r, v_r*(z+origin_[2])/r);

}

// R0 - radius at z=1pc
ConstParabolicVField::ConstParabolicVField(double gamma, double Rz0) : gamma_(gamma), Rz0_(Rz0) {}

Vector3d ConstParabolicVField::v(const Vector3d &point) const {
    double x = point[0];
    double y = point[1];
    double z = point[2];
    double r = std::hypot(x, y);
    double v_r = c*sqrt(1. - 1./(gamma_*gamma_));
    // Radius of the enclosing parabaloid at current z
    double R = sqrt(Rz0_*Rz0_*z/pc);
    // Radius of the parabaloid enclosing current point at z=1pc
    double r0 = Rz0_*r/R;
    double alpha = atan(2*pc*r/r0/r0);
    Vector3d result;
    if (r == 0) {
        result = Vector3d(0, 0, v_r);
    }
    else {
        result = Vector3d(v_r*cos(alpha)*x/r, v_r*cos(alpha)*y/r, v_r*sin(alpha));
    }
    return result;
};


AccParabolicVField::AccParabolicVField(double gamma0, double R0, double Rz0) : gamma0_(gamma0), R0_(R0), Rz0_(Rz0) {}

Vector3d AccParabolicVField::v(const Vector3d &point) const {
    double x = point[0];
    double y = point[1];
    double z = point[2];
    double r = std::hypot(x, y);

    // Radius of the enclosing parabaloid at current z
    double R = sqrt(Rz0_*Rz0_*z/pc);
    // Radius of the parabaloid enclosing current point at z=1pc
    double r0 = Rz0_*r/R;
    double alpha = atan(2*pc*r/r0/r0);

    double gamma = 1 + (gamma0_-1)/sqrt(R0_)*sqrt(z);
    double v_r = c*sqrt(1. - 1./(gamma*gamma));

    Vector3d result;
    if (r == 0) {
        result = Vector3d(0, 0, v_r);
    }
    else {
        result = Vector3d(v_r*cos(alpha)*x/r, v_r*cos(alpha)*y/r, v_r*sin(alpha));
    }
    return result;
};


// gamma(z) = 1+a*sqrt(z)
// gamma_axis0 - speed at z=R0 at radius=0
// gamma_border0 - speed at z=R0 at radius=R0_border
// Rz0 - radius of parabaloid at z=1pc
ShearedAccParabolicVField::ShearedAccParabolicVField(double gamma_axis0, double gamma_border0, double R0, double Rz0) :
        gamma_axis0_(gamma_axis0), gamma_border0_(gamma_border0), R0_(R0), Rz0_(Rz0) {}

Vector3d ShearedAccParabolicVField::v(const Vector3d &point) const {
    double x = point[0];
    double y = point[1];
    double z = point[2];
    double r = std::hypot(x, y);

    // Need to convert current radius r to radius at z=R0
    // Radius of the enclosing paraboloid at z=R0
    double R0_border = sqrt(Rz0_*Rz0_*R0_/pc);
    // Radius of the enclosing paraboloid at current z
    double R = sqrt(Rz0_*Rz0_*z/pc);
    // Radius of the paraboloid enclosing current point at z=R0
    double rR0 = R0_border*r/R;

//    std::cout << "r/R at z=R0 is " << rR0/R0_border << std::endl;

    double gamma0 = gamma_axis0_-(gamma_axis0_-gamma_border0_)*pow(rR0/R0_border, 0.25);
//    std::cout << "Thus gamma at z=R0 for given r/R is " << gamma0 << std::endl;
    if (gamma0 < 1) {gamma0 = 1;}
    double gamma = 1 + (gamma0-1)/sqrt(R0_)*sqrt(z);
//    std::cout << "Due to acceleration gamma at point is " << gamma << std::endl;
    double v_r = c*sqrt(1. - 1./(gamma*gamma));

    // Radius of the paraboloid enclosing current point at z=1pc
    double r0 = Rz0_*r/R;
    double alpha = atan(2*pc*r/r0/r0);

    Vector3d result;
    if (r == 0) {
        result = Vector3d(0, 0, v_r);
    }
    else {
        result = Vector3d(v_r*cos(alpha)*x/r, v_r*cos(alpha)*y/r, v_r*sin(alpha));
    }
    return result;
}
//ConstParabolicConstConeVField::ConstParabolicConstConeVField(double gamma, double R0, double z0) :
//z0_(z0), conev(gamma), parav(gamma, R0) {}
//
//Vector3d ConstParabolicConstConeVField::v(const Vector3d &point) const {
//    if (point[2] > z0_) {
//        return conev.v(point);
//    }
//    else {
//        return parav.v(point);
//    }
//};


AccParabolicConstConeVField::AccParabolicConstConeVField(double gamma0, double Rz0, double z0) :
        z0_(z0),
        conev(gamma0, Vector3d(0, 0, (Rz0*sqrt(z0/pc)/(0.5*Rz0*pow(z0*pc, -0.5)) - z0))),
        parav(gamma0, z0, Rz0) {
    std::cout << "Cone origin at z= " << (Rz0*sqrt(z0/pc)/(0.5*Rz0*pow(z0*pc, -0.5)) - z0)/pc << std::endl;
    std::cout << "Cone theta = " << atan(sqrt(Rz0*Rz0*z0/pc)/z0) << std::endl;

}

Vector3d AccParabolicConstConeVField::v(const Vector3d &point) const {
    if (point[2] > z0_) {
        return conev.v(point);
    }
    else {
        return parav.v(point);
    }
};


// gamma(z) = 1+a*sqrt(z)
// gamma_axis0 - maximal speed at z=z0 at radius=0
// gamma_border0 - speed at z=z0 at the border of paraboloid (at R=sqrt(Rz0*Rz0*z0/pc))
// Rz0 - radius of parabaloid at z=1pc
// z0 - z-coordinate where paraboloid goes into cone
ShearedAccParabolicConstConeVField::ShearedAccParabolicConstConeVField(double gamma_axis0, double gamma_border0,
        double Rz0, double z0) :
        z0_(z0),
        conev(gamma_axis0, gamma_border0, atan(sqrt(Rz0*Rz0*z0/pc)/z0),
                Vector3d(0, 0, (Rz0*sqrt(z0/pc)/(0.5*Rz0*pow(z0*pc, -0.5)) - z0))),
        parav(gamma_axis0, gamma_border0, z0, Rz0) {
//    std::cout << "z0 = " << z0/pc << std::endl;
//    std::cout << "Rz0 = " << Rz0/pc << std::endl;
//    std::cout << "sqrt = " << sqrt(Rz0*Rz0*z0/pc) << std::endl;
    std::cout << "Cone origin at z= " << (Rz0*sqrt(z0/pc)/(0.5*Rz0*pow(z0*pc, -0.5)) - z0)/pc << std::endl;
    std::cout << "Cone theta = " << atan(sqrt(Rz0*Rz0*z0/pc)/z0) << std::endl;
}

Vector3d ShearedAccParabolicConstConeVField::v(const Vector3d &point) const {
    if (point[2] > z0_) {
        return conev.v(point);
    }
    else {
        return parav.v(point);
    }
};


SimulationVField::SimulationVField(Delaunay_triangulation *tr) : interp_(tr) {}

Vector3d SimulationVField::v(const Vector3d &point) const {
    double gamma = interp_.interpolated_value(point);
    return Vector3d{0, 0, c * sqrt(1. - 1. / (gamma * gamma))};
};