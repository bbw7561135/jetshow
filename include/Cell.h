//
// Created by ilya on 6/2/17.
//

#ifndef JETSHOW_CELL_H
#define JETSHOW_CELL_H

#include <Eigen/Eigen>
#include <vector>
#include <map>

using Eigen::Vector3d;

class Cell {
public:
    Cell(std::vector<Vector3d> &points);
    std::vector<Cell> split_n(int n);
    std::vector<Vector3d> points();

private:
    std::vector<Vector3d> points_;
    Vector3d center_;
    double length_;
    std::map<std::string,double> values_;

};


#endif //JETSHOW_CELL_H
