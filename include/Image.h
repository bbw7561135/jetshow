//
// Created by ilya on 6/5/17.
//

#ifndef JETSHOW_IMAGE_H
#define JETSHOW_IMAGE_H

#include <array>
#include <string>
#include <Eigen/Eigen>
#include <memory>
#include "Pixel.h"

using std::array;
using std::vector;
using Eigen::Vector3d;


class Image {
public:
//    Image(const Image&) = delete;
//    Image& operator=(const Image&) = delete;
//    ~Image() = default;
//    explicit Image(pair<int, int> image_size, double pixel_size,
//                   double pixel_scale);
		Image(pair<int, int> image_size, double pixel_size, double pixel_scale);
    const int num_of_pixels_;
    const pair<int,int> image_size;
    vector<vector<double>>& getImage(string value);
    Vector3d getCoordinate(int i, int j);
    Vector3d getScaledCoordinate(int i, int j);
    vector<Pixel>& getPixels();

private:
    pair<int,int> image_size_;
    double pixel_size_;
    double pixel_scale_;
    vector<Pixel> pixels_;

};


#endif //JETSHOW_IMAGE_H
