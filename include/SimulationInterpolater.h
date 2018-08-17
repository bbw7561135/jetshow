//
// Created by ilya on 18.08.18.
//

#ifndef JETSHOW_SIMULATIONINTERPOLATER_H
#define JETSHOW_SIMULATIONINTERPOLATER_H


#include <Eigen/Eigen>

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


class SimulationInterpolater {
public:
    explicit SimulationInterpolater(Delaunay_triangulation *tr);
    double interpolated_value(Vector3d point);

private:
    Delaunay_triangulation* tr_;
};


#endif //JETSHOW_SIMULATIONINTERPOLATER_H
