#include "Pixel.h"

using std::pair;


Pixel::Pixel(double size, Vector3d &coordinate, pair<int,int> &ij): size_(size),
                                                 coordinate_(coordinate),
                                                 ij_(ij) {};

void Pixel::scale_coordinates(double scale) {
    coordinate_ *= scale;
}
Vector3d &Pixel::getCoordinate() {
  return coordinate_;
}

double Pixel::getValue(string value) {
  return values_[value];
}
