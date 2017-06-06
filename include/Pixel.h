//
// Created by ilya on 6/5/17.
//

#ifndef JETSHOW_PIXEL_H
#define JETSHOW_PIXEL_H

#include <array>
#include <string>
#include <bits/unordered_map.h>

using std::array;
using std::unordered_map;
using std::string;

class Pixel {
public:
    Pixel(double size_, const array<double, 2> &coordinate_,
          const array<int, 2> &ij_,
          const unordered_map<string, double> &values_);

    void scale_coordinates(double scale_x, double scale_y);
    void scale_coordinates(double scale);

private:
    double size_;
    pair<double,double> coordinate_;
    pair<int,int> ij_;
    unordered_map<string,double> values_;


};


#endif //JETSHOW_PIXEL_H
