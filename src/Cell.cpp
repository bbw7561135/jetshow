//
// Created by ilya on 6/2/17.
//

#include "Cell.h"
#include "linspace.h"


Cell::Cell(std::vector<Vector3d> &points) {
    points_ = points;
    center_ = 0.5 * (points[0] + points[1]);
    length_ = (points[1] - points[0]).norm();
    values_ = {};
}

std::vector<Vector3d> Cell::points() {
    return points_;
}

std::vector<Cell> Cell::split_n(int n) {
    std::vector<Eigen::Vector3d> new_points = linspace(points_[0], points_[1], n);
    std::vector<Cell> cells;
    for (std::size_t i=1; i != new_points.size(); ++i) {
        std::vector<Vector3d> cell_points = {new_points[i-1], new_points[i]};
        Cell new_cell = Cell(cell_points);
        cells.push_back(new_cell);
    }
    return cells;
}