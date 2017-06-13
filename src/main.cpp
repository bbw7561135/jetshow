#include <iostream>
#include <boost/numeric/odeint.hpp>
#include "BField.h"
#include "VField.h"
//#include "Ray.h"
#include "Cone.h"
#include "Jet.h"
#include "utils.h"
#include <math.h>
#include <linspace.h>
#include <Cell.h>
#include <NField.h>
#include "Cylinder.h"



using Eigen::Vector3d;
using std::vector;
using namespace boost::numeric::odeint;


class Rhs {
 public:
  Rhs(Jet* newjet, Vector3d &newpoint_in, Vector3d &newray_direction,
      double newnu) {
    jet = newjet;
    point_in = newpoint_in;
    ray_direction = newray_direction;
    nu = newnu;
  }
  void operator() (const double &x, double &dxdt, const double t) {
    Vector3d point = point_in + t * ray_direction;
    std::cout << "calculating k_I at point " << point << std::endl;
    dxdt = jet->getKI(point, ray_direction, nu);
  }
  void setPointIn(Vector3d &newpoint_in) {point_in = newpoint_in;}
  void setRayDirection(Vector3d &newray_direction) {ray_direction = newray_direction;}
 private:
  Jet* jet;
  Vector3d point_in;
  Vector3d ray_direction;
  double nu;
};

void write_cout(const double &x, const double t) {
  std::cout << t << '\t' << x << std::endl;
}

void test_jet() {
    // Create cone geometry of jet
    Vector3d origin = {0., 0., 0.};
    Vector3d direction = {0., 0., 1.};
    double angle = pi/4.;
    double scale = 10.;
    Cone geometry(origin, direction, angle, scale);

    RadialConicalBField bfield(0.1, 2.0);
    ConstCenralVField vfield(10.);
    BKNField nfield(1., 10.);

    Jet bkjet(&geometry, &vfield, &bfield, &nfield);
    // Test just one point
//    Vector3d test_point{0., 0.3, 3.};
//    std::cout << bkjet.getV(test_point) << std::endl;
//    std::cout << bkjet.getN(test_point) << std::endl;
//    std::cout << bkjet.getB(test_point) << std::endl;


    // Get a ray
    Vector3d ray_direction(0., 1., 0.);
    Vector3d ray_origin(0., 1., 1.);
    Ray ray(ray_origin, ray_direction);
    std::list<Intersection> list_intersect = geometry.hit(ray);
    Intersection intersect = list_intersect.front();
    std::pair<Vector3d,Vector3d> borders = intersect.get_path();
    Vector3d point_in = borders.first;
    Vector3d point_out = borders.second;
    double length = (point_out - point_in).norm();
    double dt = length/100.;
    double nu = 5.*pow(10., 9.);

    Rhs rhs(&bkjet, point_in, ray_direction, nu);


    typedef runge_kutta_dopri5< double > stepper_type;

    // Here x - optical depth \tau
    double x = 0.0;
    integrate_adaptive(make_controlled(1E-18, 1E-18, stepper_type()), rhs, x,
                       0.0 , length, dt, write_cout);

}


int main() {
    test_jet();

//  Intersection intersection{};
//  std::cout << intersection.direction() << std::endl;
//  auto borders = intersection.get_path();
//  std::cout << intersection.direction() << std::endl;
//  std::cout << borders.first << std::endl;
//    Vector3d cone_origin(0., 0., 0.);
//    Vector3d cone_direction(0., 0., 1.);
//    double cone_angle = pi/4.;
//    double scale = 10.0;
//    Cone cone = Cone(cone_origin, cone_direction, cone_angle, scale);
  // There's two intersections
//    Vector3d ray_direction(0., 1., 0.);
//    Vector3d ray_origin(0., 1., 1.);
  // No intersections
//    Vector3d ray_direction(0., 1., 0.);
//    Vector3d ray_origin(1., 1., 0.);
  // Two intersection (one far away)
//    Vector3d ray_direction(0., 1., 1.);
//    Vector3d ray_origin(0., 1., 0.);
  // One intersection
//    Vector3d ray_direction(0., 1., 0.);
//    Vector3d ray_origin(0., 1., 0.);
    // Along border
//    Vector3d ray_origin(0., 0., 0.);
//    Vector3d ray_direction(0., 1., 1.);
//
//    Ray ray(ray_origin, ray_direction);
//    std::list<Intersection> list_intersect = cone.hit(ray);
//    std::cout << "Did ray traverse volume ?" << std::endl << !list_intersect.empty() << std::endl;
//    if (!list_intersect.empty()) {
//        Intersection intersect = list_intersect.front();
//        std::pair<Vector3d,Vector3d> borders = intersect.get_path();
//        std::cout << "Borders : " << borders.first << " and " << borders.second << std::endl;
//    }


      // Test ``is_within``
//    Vector3d point{-2., 1., -5.};
//    bool is_within = cone.is_within(point);
//    std::cout << "Is within : " << is_within << std::endl;

    // This tests cylinder-ray intersections
//    Vector3d cylinder_origin(0., 0., 0.);
//    Vector3d cylinder_direction(0., 0., 1.);
//    double cylinder_radius = 1.;
//    Cylinder cylinder(cylinder_origin, cylinder_direction, cylinder_radius);
//    Vector3d point{0., 0., 2.};
//    bool is_within = cylinder.is_within(point);
//    std::cout << "Is within : " << is_within << std::endl;

    // No intersections
//    Vector3d ray_direction(1., 0., 1.);
//    Vector3d ray_origin(0., 2., 0.);

    // Along border
//    Vector3d ray_direction(0., 0., 1.);
//    Vector3d ray_origin(0., 1., 0.);

    // No interception. Along border internally
//    Vector3d ray_direction(0., 0., 1.);
//    Vector3d ray_origin(-0.1, 0.2, 0.);

    // No interception. Along border externally
//    Vector3d ray_direction(0., 0., 1.);
//    Vector3d ray_origin(10., -0.2, 0.);

    // Single interception
//    Vector3d ray_direction(1., 0., 0.);
//    Vector3d ray_origin(0., 1., 0.);

//    // Two intersections
//    Vector3d ray_direction(0., 1., 0.);
//    Vector3d ray_origin(0., 1., 0.);
//
//    Ray ray(ray_origin, ray_direction);
//    std::list<Intersection> list_intersect = cylinder.hit(ray);
//    std::cout << "Did ray traverse volume ?" << std::endl << !list_intersect.empty() << std::endl;
//    if (!list_intersect.empty()) {
//        Intersection intersect = list_intersect.front();
//        std::pair<Vector3d,Vector3d> borders = intersect.get_path();
//        std::cout << "Borders : " << borders.first << " and " << borders.second << std::endl;
//        std::cout << "Direction " << intersect.direction() << std::endl;
//    }
    // End of cylinder-ray intersections



//    Matrix3d eye_matrix;
//    eye_matrix << 1, 0, 0,
//                  0, 1, 0,
//                  0, 0, 1;
//    // DP
//    Vector3d delta = ray_origin - cone_origin;
//    std::cout << delta << std::endl;
//    // M
//    Matrix3d M = cone_direction * cone_direction.transpose() - cos(cone_angle)*cos(cone_angle)*eye_matrix;
//    std::cout << M << std::endl;
//    double c2 = ray_direction.transpose() * M * ray_direction;
//    double c1 = ray_direction.transpose() * M * delta;
//    double c0 = delta.transpose() * M * delta;
//    std::cout << "c2 = " << c2 << " c1 = " << c1 << " c0 = " << c0 << std::endl;
//    double d = c1*c1 - c0*c2;
//    std::cout << "d = " << d << std::endl;
//    std::cout << "-c1-sqrt(d) = " << -c1 - sqrt(d) << std::endl;
//    double t1 = (-c1 + sqrt(d)) / (c2);
//    double t2 = (-c1 - sqrt(d)) / (c2);
//    std::cout << "t1= " << t1 << " t2= " << t2 << std::endl;

}

//int main() {
//    Vector3d v1(0, 0, 0);
//    Vector3d v2(10, 10, -10);
//    std::vector<Vector3d> vec_1 = linspace(v1, v2, 11);
//    print_vector(vec_1);
//    return 0;
//}

//int main() {
//    Vector3d v1(0, 0, 0);
//    Vector3d v2(10, 10, 10);
//    std::vector<Vector3d> v = {v1, v2};
//    Cell cell = Cell(v);
//    std::vector<Cell> cells = cell.split_n(11);
//    for (size_t i=0; i != cells.size();++i) {
//        std::cout << "Cell number" << i << std::endl;
//        Print(cells[i].points());
//    }
//    return 0;
//
//}



//int main() {
//    Vector3d origin = {0., 0., 0.};
//    Vector3d direction = {1., 1., 1.};
//    direction.normalize();
//    Ray ray(origin, direction);
//    Vector3d point = ray.point(10);
//    std::cout << "Origin:" << std::endl;
//    std::cout << ray.origin() << std::endl;
//    std::cout << "Direction:" << std::endl;
//    std::cout << ray.direction() << std::endl;
//
//    std::cout << "Point:" << std::endl;
//    std::cout << point << std::endl;
//}


/*
int main() {
    BField bfield(1., 10., 10.);
    std::cout << "Bz-Field at z0 : " << bfield.getZ0() << std::endl;
    Vector3d p = {0.1, 0.1, 10};
    Vector3d b;
    b = bfield.bf(p);
    std::cout << "B-Field at z = 10 : " << b.norm();
    return 0;
}*/

