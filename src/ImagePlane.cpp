#include <Eigen/Eigen>
#include <memory>
#include <math.h>
#include <iostream>
#include "ImagePlane.h"

using Eigen::Vector3d;
using std::pair;


ImagePlane::ImagePlane(pair<int,int> image_size, double pixel_size,
                       double pixel_scale,
                       double los_angle) : image_(image_size, pixel_size,
                                                  pixel_scale),
                                           los_angle_(los_angle),
                                           image_size(image_size) {
  direction_ = Vector3d{-sin(los_angle), 0, -cos(los_angle)};
	Vector3d scale(1., 1., 1./sin(los_angle));
  // Initialize rays
  vector<Pixel>& pixels = image_.getPixels();
  for (int i = 0; i < image_size.first; ++i) {
    for (int j = 0; j < image_size.second; ++j) {
      Pixel pxl = pixels[i * image_size.first + j];
			// Pixels of image has coordinates of observer image (rotated on angle
			// (pi/2-alpha) around y-axis of original xyz jet frame.
      Vector3d coordinate = pxl.getCoordinate();
			// Rays have coordinates in (yz)-plane of the jet.
			coordinate = scale.array()*coordinate.array();
      auto ij = std::make_pair(i, j);
      auto ray = Ray(coordinate, direction_);
//      std::cout << "Setting ray i, j = " << i << " " << j << "with coordinate = " << coordinate << std::endl;
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

vector<vector<double>> ImagePlane::getImage(string value) {
  return image_.getImage(value);
}

vector<double> ImagePlane::getStripe(string value) {
	return image_.getStripe(value);
}
