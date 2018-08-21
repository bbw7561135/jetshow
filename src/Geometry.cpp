#include "Geometry.h"
#include "utils.h"
#include <math.h>
#include <iostream>
#include <cnpy.h>
#include <list>
#include "Intersection.h"
#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>

using Eigen::Vector3d;

typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_3 Point_3;
typedef K::Ray_3 Ray_3;
typedef K::Vector_3 Vector_3;
typedef CGAL::Polyhedron_3<K> Polyhedron;
typedef CGAL::AABB_face_graph_triangle_primitive<Polyhedron> Primitive;
typedef CGAL::AABB_traits<K, Primitive> Traits;
typedef CGAL::AABB_tree<Traits> Tree;
typedef Tree::Primitive_id Primitive_id;
typedef boost::optional<Tree::Intersection_and_primitive_id<Ray_3>::Type> Ray_intersection;


std::list<double>
intersection(Vector3d R0, Vector3d Rd, double A, double B, double C, double D, double E, double F, double G, double H,
             double I, double J) {
    std::list<double> result;

    double x0 = R0[0];
    double y0 = R0[1];
    double z0 = R0[2];

    double xd = Rd[0];
    double yd = Rd[1];
    double zd = Rd[2];

    double Aq = A*xd*xd + B*yd*yd + C*zd*zd + D*xd*yd + E*xd*zd + F*yd*zd;

    double Bq = 2*A*x0*xd + 2*B*y0*yd + 2*C*z0*zd + D*(x0*yd + y0*xd) + E*x0*zd + F*(y0*zd + yd*z0) + G*xd + H*yd +
            I*zd;
    double Cq = A*x0*x0 + B*y0*y0 + C*z0*z0 + D*x0*y0 + E*x0*z0 + F*y0*z0 + G*x0 + H*y0 + I*z0 + J;
    double Dscr = Bq*Bq - 4*Aq*Cq;

    if (Aq == 0) {
        result = std::list<double>{-Cq/Bq};
    }

    else if (Dscr < 0) {
        result = std::list<double>{};
    }

    else {
        double t0 = (-Bq - sqrt(Dscr))/(2*Aq);

        if (t0 > 0) {
            result = std::list<double>{t0};
        }
        else {
            double t1 = (-Bq + sqrt(Dscr))/(2*Aq);
            result = std::list<double>{t0, t1};
        }
    }

    return result;
}


std::pair<Vector3d, Vector3d> Geometry::full_infinite_path(Ray &ray) const {
  Vector3d diff = ray.origin() - origin();
  Vector3d ref_point = diff - diff.dot(ray.direction()) * ray.direction();
  Vector3d point_in = ref_point - big_scale() * ray.direction();
  Vector3d point_out = ref_point + big_scale() * ray.direction();
  return std::pair<Vector3d,Vector3d>{point_in, point_out};
}

std::pair<Vector3d, Vector3d> Geometry::half_infinite_path(Ray &ray,
                                                           const Vector3d &point) const {
  Vector3d check_point = point-ray.direction();
  if (is_within(check_point)) {
    return std::pair<Vector3d,Vector3d>{point-big_scale()*ray.direction(), point};

  } else {
    return std::pair<Vector3d,Vector3d>{point, point+big_scale()*ray.direction()};
  }
}


SimulationGeometry::SimulationGeometry(Tree *tree)  {
    tree_ = tree;
}

//void SimulationGeometry::create_tree(double *loaded_data, int nrows) {
//
//    std::vector<Point> points;
//    for (int i=0; i<nrows; i++) {
//        double z = loaded_data[i*3]/pc;
//        double r_p = loaded_data[i*3 + 1]/pc;
//        for (int j=0; j<48; j++) {
//            double x = r_p*sin(j*2*pi/48);
//            double y = r_p*cos(j*2*pi/48);
//            points.emplace_back(Point(x, y, z));
//        }
//    }
//    Polyhedron P;
//    CGAL::convex_hull_3(points.begin(), points.end(), P);
//    tree_ = new Tree(faces(P).first, faces(P).second, P);
//}

//void SimulationGeometry::add_simulation_data(std::string fn) {
//    std::cout << "Loading data from file" << std::endl;
//    cnpy::NpyArray arr = cnpy::npy_load(fn);
//    double* loaded_data = arr.data<double>();
//    size_t nrows = arr.shape[0];
//    create_tree(loaded_data, nrows);
//    std::cout << "Created tree" << std::endl;
//
//
//    Vector3d origin = Vector3d(0, 3, 600);
//    Vector3d direction = Vector3d(0, -1, 0);
//    Vector_3 v(direction[0], direction[1], direction[2]);
//    Point_3 p(origin[0], origin[1], origin[2]);
//    Ray_3 ray_query(p, v);
//
//
//    bool do_inter = tree_->do_intersect(ray_query);
//    std::cout << "Bool intersection inside create_tree= " << do_inter << std::endl;
//}

std::list<Intersection> SimulationGeometry::hit(Ray &ray) const {
//    std::cout << "Hit ray with origin = " << ray.origin() << std::endl;
    Vector3d origin = ray.origin()/pc;
    Vector3d direction = ray.direction();
//    std::cout << "Ray origin/pc = " << origin << std::endl;
//    std::cout << "Ray direction = " << direction << std::endl;
    Vector_3 v(direction[0], direction[1], direction[2]);
    Point_3 p(origin[0], origin[1], origin[2]);
//    Ray_3 ray_query(p, v);
    Line_3 line_query(p, v);
//    std::cout << "Ray query = " << ray_query << std::endl;

    // tests intersections with ray query
    if(tree_->do_intersect(line_query)) {
//        std::cout << "Got intersection!===============================" << std::endl;
        // computes all intersections with ray query (as pairs object - primitive_id)
        std::list<Ray_intersection> intersections;
        std::vector<Vector3d> points;
        tree_->all_intersections(line_query, std::back_inserter(intersections));
        for (auto const& intersection : intersections)
        {
            const Point_3* pt = boost::get<Point_3>(&(intersection->first));
            double x = CGAL::to_double(pt->x());
            double y = CGAL::to_double(pt->y());
            double z = CGAL::to_double(pt->z());
            points.emplace_back(Vector3d(x, y, z));
        }
        return std::list<Intersection>{Intersection(ray, points[1]*pc, points[0]*pc)};
    }
    else
//        std::cout << "No intersections" << std::endl;
        return std::list<Intersection>{};

}