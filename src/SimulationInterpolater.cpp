//
// Created by ilya on 18.08.18.
//

#include "SimulationInterpolater.h"


SimulationInterpolater::SimulationInterpolater(Delaunay_triangulation *tr) {
    tr_ = tr;
}

double SimulationInterpolater::interpolated_value(Vector3d point) {
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
        std::cout << "B_p data:\t" << fh->vertex(i)->info() << std::endl;
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