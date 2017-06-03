#include <iostream>
//#include <bfields.h>
//#include "Ray.h"
#include "Cone.h"
#include <math.h>
#include <boost/math/constants/constants.hpp>
#include <linspace.h>
#include <Cell.h>


using Eigen::Vector3d;
using std::vector;

const double pi = boost::math::constants::pi<double>();


//int main() {
//    Vector3d origin(0., 0., 0.);
//    Vector3d direction(0., 0., 1.);
//    double angle = pi/4.;
//    std::cout << "Angle: " << angle << std::endl;
//    Cone cone = Cone(origin, direction, angle);
//    // There's two intersections
////    Vector3d ray_direction(0., 1., 0.);
////    Vector3d ray_origin(0., 1., 1.);
//    // No intersections
////    Vector3d ray_direction(0., 1., 0.);
////    Vector3d ray_origin(1., 1., 0.);
//    // One intersection
//    Vector3d ray_direction(0., 1., 1.);
//    Vector3d ray_origin(0., 1., 0.);
//    Ray ray(ray_origin, ray_direction);
//    Intersection intersect = cone.hit(ray);
//    std::cout << "Has intersection?: " << intersect.has_intersection() << std::endl;
////    for(auto i=0;i < intersect.points().size();i++) {
////        std::cout << intersect.points()[i] << std::endl;
////    }
//    Print(intersect.points());
//}

//int main() {
//    Vector3d v1(0, 0, 0);
//    Vector3d v2(10, 10, -10);
//    std::vector<Vector3d> vec_1 = linspace(v1, v2, 11);
//    print_vector(vec_1);
//    return 0;
//}

int main() {
    Vector3d v1(0, 0, 0);
    Vector3d v2(10, 10, 10);
    std::vector<Vector3d> v = {v1, v2};
    Cell cell = Cell(v);
    std::vector<Cell> cells = cell.split_n(11);
    for (size_t i=0; i != cells.size();++i) {
        std::cout << "Cell number" << i << std::endl;
        Print(cells[i].points());
    }
    return 0;

}

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

