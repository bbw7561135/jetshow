#include <Eigen/Eigen>
#include <math.h>
#include <memory>
#include "ImagePlane.h"

using Eigen::Vector3d;
using std::pair;


ImagePlane::ImagePlane(pair<int,int> image_size, double pixel_size,
                       double pixel_scale,
                       double los_angle) : image_(image_size, pixel_size,
                                                  pixel_scale),
                                           los_angle_(los_angle),
                                           image_size(image_size) {
  direction_ = Vector3d{sin(los_angle), 0, -cos(los_angle)};
  // Initialize rays
  vector<std::unique_ptr<Pixel>> pixels = std::move(image_.getPixels());
  for (int i = 0; i < image_size.first; ++i) {
    for (int j = 0; j < image_size.second; ++j) {
      Vector3d coordinate = pixels[i*image_size.first+j].get()->getCoordinate();
      auto ij = std::make_pair(i, j);
      auto ptr = std::make_unique<Ray>(coordinate, direction_);
      rays_.push_back(std::move(ptr));
    }
  }

}

vector<std::unique_ptr<Pixel>> &ImagePlane::getPixels() {
  return image_.getPixels();
}

vector<std::unique_ptr<Ray>> &ImagePlane::getRays() {
  return rays_;
}

vector<vector<double>> &ImagePlane::getImage(string value) {
  return image_.getImage(value);
}
