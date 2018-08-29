#include <array>
#include <string>
#include <memory>
#include <iostream>
#include "Image.h"

using std::array;
using std::vector;
using std::string;
using Eigen::Vector3d;


vector<vector<double>> Image::getImage(string value) {
	vector<vector<double>> image;
	image.resize(image_size_.first);
	for (unsigned long int i = 0; i < image_size_.first; ++i) {
		image[i].resize(image_size_.second);
	}
  for (unsigned long int i = 0; i < image_size_.first; ++i) {
    for (unsigned long int j = 0; j < image_size_.second; ++j) {
      image[i][j] = pixels_[i*image_size_.second + j].getValue(value);
    }
  }
	return image;
}


vector<double> Image::getStripe(string value) {
	vector<double> stripe;
	for (unsigned long int i = 0; i < image_size.first; ++i) {
		for (unsigned long int j = image_size.second / 2; j < image_size.second; ++j) {
			if (i == image_size.first / 2) {
				stripe.push_back(pixels_[i*image_size_.second + j].getValue(value));
			}
		}
	}
	return stripe;
};


Image::Image(pair<unsigned long int, unsigned long int> image_size, double pixel_size,
             double pixel_scale):
        image_size_(image_size),
        pixel_size_(pixel_size),
        pixel_scale_(pixel_scale),
        num_of_pixels_(image_size.first*image_size.second),
        pixels_() {
  for (unsigned long int i = 0; i < image_size_.first; ++i) {
    for (unsigned long int j = 0; j < image_size_.second; ++j) {
      Vector3d coordinate = getScaledCoordinate(i, j);
      auto ij = std::make_pair(i, j);
      auto pxl = Pixel(pixel_size, coordinate, ij);
      pixels_.push_back(std::move(pxl));
    }
  }
}

// Original
//Vector3d Image::getCoordinate(unsigned long int i, unsigned long int j) {
//  return Vector3d{0, i-image_size_.first/2.+0.5, j-image_size_.second/2.+0.5};
//}

// Along jet coordinate starts from 0 to image_size.second.
// Addition to ``j``-part equals ``dr_min_pc/pix_size_pc*sin(theta_LOS)``, where
// ``dr_min`` - minimal distance in simulation output
Vector3d Image::getCoordinate(unsigned long int i, unsigned long int j) {
    return Vector3d{0, i-image_size_.first/2.+0.5, j+0.5+90};
}

Vector3d Image::getScaledCoordinate(unsigned long int i, unsigned long int j) {
  return pixel_scale_*getCoordinate(i, j);
}

vector<Pixel> &Image::getPixels() {
  return pixels_;
};
