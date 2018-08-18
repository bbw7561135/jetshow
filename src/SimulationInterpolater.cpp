//
// Created by ilya on 18.08.18.
//

#include "SimulationInterpolater.h"
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


SimulationInterpolater::SimulationInterpolater(Delaunay_triangulation *tr) {
    tr_ = tr;
}

double SimulationInterpolater::interpolated_value(Vector3d point) const {
    // Conver 3D point (Vector3d) to (r, r_p) coordinates (Point_)
    double x = point[0];
    double y = point[1];
    double z = point[2];
    double r_p = hypot(x, y);
    Point_ pt(z, r_p);

    Delaunay_triangulation::Face_handle fh = tr_->locate(pt);

    std::vector<Point_ > vertexes;
    std::vector<double> info;

    for (int i=0; i<3; i++) {
        vertexes.push_back(fh->vertex(i)->point());
        info.push_back(fh->vertex(i)->info());

        std::cout << "Triangle:\t" << tr_->triangle(fh) << std::endl;
        std::cout << "Vertex 0:\t" << tr_->triangle(fh)[i] << std::endl;
        std::cout << "Value:\t" << fh->vertex(i)->info() << std::endl;
    }

    // Create an std::vector to store coordinates.
    Scalar_vector coordinates;
    // Instantiate the class Triangle_coordinates_2 for the triangle defined above.
    Triangle_coordinates triangle_coordinates(vertexes[0], vertexes[1], vertexes[2]);

    triangle_coordinates(pt, std::inserter(coordinates, coordinates.end()));

    double interpolated_value = 0;
    for(int j = 0; j < 3; ++j) {
        std::cout << "coordinate " << j + 1 << " = " << coordinates[j] << "; ";
        interpolated_value += coordinates[j]*info[j];
    }
    std::cout << "Interpolated value = " << interpolated_value << std::endl;
    return interpolated_value;
}