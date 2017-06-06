//
// Created by ilya on 6/2/17.
//

#ifndef JETSHOW_LINSPACE_H
#define JETSHOW_LINSPACE_H

#include <iostream>
#include <vector>
#include <Eigen/Eigen>

// https://stackoverflow.com/a/27030598
template<typename T>
std::vector<Eigen::Vector3d> linspace(T start_in, T end_in, int num_in)
{

    std::vector<Eigen::Vector3d> linspaced;

    Eigen::Vector3d start = static_cast<Eigen::Vector3d>(start_in);
    Eigen::Vector3d end = static_cast<Eigen::Vector3d>(end_in);
    double num = static_cast<double>(num_in);

    if (num == 0) { return linspaced; }
    if (num == 1)
    {
        linspaced.push_back(start);
        return linspaced;
    }

    Eigen::Vector3d delta = (end - start) / (num - 1);

    for(int i=0; i < num-1; ++i)
    {
        linspaced.push_back(start + delta * i);
    }
    linspaced.push_back(end); // I want to ensure that start and end
    // are exactly the same as the input
    return linspaced;
}


#endif //JETSHOW_LINSPACE_H
