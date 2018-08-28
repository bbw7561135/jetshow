//
// Created by ilya on 18.08.18.
//

#include "SimulationInterpolater.h"
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Interpolation_traits_2.h>
#include <CGAL/natural_neighbor_coordinates_2.h>
#include <CGAL/interpolation_functions.h>
#include <CGAL/Barycentric_coordinates_2/Triangle_coordinates_2.h>
#include <SimulationInterpolater.h>
#include <cnpy.h>
#include "utils.h"


#define CGAL_HAS_THREADS

typedef CGAL::Simple_cartesian<double>                                   K_;
typedef K_::Point_2                                                Point_;
typedef CGAL::Triangulation_vertex_base_with_info_2<double, K_>      Vb;
typedef CGAL::Triangulation_data_structure_2<Vb>                  Tds;
typedef CGAL::Delaunay_triangulation_2<K_, Tds>                    Delaunay_triangulation;
typedef K_::FT                                               Coord_type;
typedef std::vector<Coord_type >                            Scalar_vector;
typedef CGAL::Barycentric_coordinates::Triangle_coordinates_2<K_> Triangle_coordinates;


SimulationInterpolater::SimulationInterpolater(Delaunay_triangulation *tr, double nan_value) {
    tr_ = tr;
    std::cout << "Initializing interpolater with nan_value = " << nan_value << std::endl;
    nan_value_ = nan_value;
    previous_hit_fh_ = nullptr;

//    for (Delaunay_triangulation::Finite_faces_iterator fit = tr->finite_faces_begin(); fit != tr->finite_faces_end(); ++fit) {
//        Delaunay_triangulation::Face_handle fh = fit;
//        previous_hit_fh_ = std::make_shared<Delaunay_triangulation::Face_handle>(fh);
//        break;
//    }
}

double SimulationInterpolater::interpolated_value(Vector3d point) const {
    // Conver 3D point (Vector3d) to (r, r_p) coordinates (Point_)
    double x = point[0]/pc;
    double y = point[1]/pc;
    double z = point[2]/pc;
    double r_p = hypot(x, y);
    Point_ pt(z, r_p);

    Delaunay_triangulation::Face_handle fh;
    // TODO: Try T.inexact_locate for speed
    if (previous_hit_fh_ != nullptr) {
//        std::cout << "Hit with hint" << std::endl;
        fh = tr_->locate(pt, previous_hit_fh_);
//        std::cout << "Done Hit with hint" << std::endl;
    } else {
//        std::cout << "First time hit" << std::endl;
        fh = tr_->locate(pt);
    }

    // This previous w/o hints
//    fh = tr_->locate(pt);


    if (tr_->is_infinite(fh)) {
        previous_hit_fh_ = nullptr;
        return nan_value_;
    } else {
        previous_hit_fh_ = fh;
    }

    std::vector<Point_ > vertexes;
    std::vector<double> info;

    for (int i=0; i<3; i++) {
        vertexes.push_back(fh->vertex(i)->point());
        info.push_back(fh->vertex(i)->info());

//        std::cout << "Triangle:\t" << tr_->triangle(fh) << std::endl;
//        std::cout << "Vertex 0:\t" << tr_->triangle(fh)[i] << std::endl;
//        std::cout << "Value:\t" << fh->vertex(i)->info() << std::endl;
    }

    // Create an std::vector to store coordinates.
    Scalar_vector coordinates;
    // Instantiate the class Triangle_coordinates_2 for the triangle defined above.
    Triangle_coordinates triangle_coordinates(vertexes[0], vertexes[1], vertexes[2]);

    triangle_coordinates(pt, std::inserter(coordinates, coordinates.end()));

    double interpolated_value = 0;
    for(int j = 0; j < 3; ++j) {
        interpolated_value += coordinates[j]*info[j];
    }

    if (std::isnan(interpolated_value)) {
        interpolated_value = nan_value_;
    }

    return interpolated_value;
}


void create_triangulation(std::string fn, Delaunay_triangulation *tr) {
    std::vector< std::vector<double> > all_points;
    read_from_txt(fn, all_points);
    size_t nrows = all_points.size();


    std::vector< std::pair<Point_,double> > points;
    for (int i=0; i<nrows; i++) {
        Point_ pt(all_points[i][0]/pc, all_points[i][1]/pc);
        points.emplace_back( std::make_pair( pt,  all_points[i][2]) );
    }
    tr->insert(points.begin(), points.end());
}
