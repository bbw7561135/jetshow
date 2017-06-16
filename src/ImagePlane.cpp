#include <Eigen/Eigen>
#include <memory>
#include <math.h>
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
  vector<Pixel>& pixels = image_.getPixels();
  for (int i = 0; i < image_size.first; ++i) {
    for (int j = 0; j < image_size.second; ++j) {
      Pixel pxl = pixels[i * image_size.first + j];
      Vector3d& coordinate = pxl.getCoordinate();
      auto ij = std::make_pair(i, j);
      auto ray = Ray(coordinate, direction_);
      rays_.push_back(std::move(ray));
    }
  }

}

vector<Pixel> &ImagePlane::getPixels() {
  return image_.getPixels();
}

vector<Ray> &ImagePlane::getRays() {
  return rays_;
}

vector<vector<double>> &ImagePlane::getImage(string value) {
  return image_.getImage(value);
}
