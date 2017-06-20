#include <iostream>
#include "Pixel.h"

using std::pair;


Pixel::Pixel(double size, Vector3d &coordinate, pair<int,int> &ij): size_(size),
                                                 coordinate_(coordinate),
                                                 ij_(ij) {};

void Pixel::scale_coordinates(double scale) {
    coordinate_ *= scale;
}

Vector3d Pixel::getCoordinate() {
  return coordinate_;
}

double Pixel::getValue(string value) {
//	std::cout << "Getting value from Pixel #..." << ij_.first
//						<< ij_.second << std::endl;
  return values_[value];
}

void Pixel::setValue(string value, double newvalue) {
//	std::cout << "Setting value of Pixel #..." << ij_.first
//						<< ij_.second << " with " << newvalue << std::endl;
  values_[value] = newvalue;
}
