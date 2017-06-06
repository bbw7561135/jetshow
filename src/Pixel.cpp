//
// Created by ilya on 6/5/17.
//

#include "Pixel.h"


Pixel::Pixel(double size, const array<double, 2> &coordinate,
             const array<int, 2> &ij,
             const unordered_map<string, double> &values): size_(size),
                                                            coordinate_(coordinate),
                                                            ij_(ij), values_(values){}

void Pixel::scale_coordinates(double scale_x, double scale_y) {
    coordinate_[0] *= scale_x;
    coordinate_[1] *= scale_y;
}

void Pixel::scale_coordinates(double scale) {
    coordinate_[0] *= scale;
    coordinate_[1] *= scale;
}
