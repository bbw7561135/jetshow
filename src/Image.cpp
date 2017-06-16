#include <array>
#include <string>
#include <memory>
#include <iostream>
#include "Image.h"

using std::array;
using std::vector;
using std::string;
using Eigen::Vector3d;


vector<vector<double>>& Image::getImage(string value) {
    vector<vector<double>> image;
  for (int i = 0; i < image_size_.first; ++i) {
    for (int j = 0; j < image_size_.second; ++j) {
      image[i][j] = pixels_[i*image_size_.first + j].getValue(value);
    }
  }
}

Image::Image(pair<int, int> image_size, double pixel_size,
             double pixel_scale):
        image_size_(image_size),
        image_size(image_size),
        pixel_size_(pixel_size),
        pixel_scale_(pixel_scale),
        num_of_pixels_(image_size_.first*image_size_.second),
        pixels_() {
  for (int i = 0; i < image_size_.first; ++i) {
    for (int j = 0; j < image_size_.second; ++j) {
      Vector3d coordinate = getScaledCoordinate(i, j);
      auto ij = std::make_pair(i, j);
      auto pxl = Pixel(pixel_size, coordinate, ij);
      pixels_.push_back(std::move(pxl));
    }
  }
}

Vector3d Image::getCoordinate(int i, int j) {
  return Vector3d(0., i-image_size_.first/2.+0.5, j-image_size_.second/2.+0.5);
}

Vector3d Image::getScaledCoordinate(int i, int j) {
  return pixel_scale_*getCoordinate(i, j);
}

vector<Pixel> &Image::getPixels() {
  return pixels_;
};
