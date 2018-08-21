#ifndef JETSHOW_PIXEL_H
#define JETSHOW_PIXEL_H

#include <array>
#include <string>
#include <map>
#include <Eigen/Eigen>

using std::array;
using std::map;
using std::string;
using std::pair;
using Eigen::Vector3d;


class Pixel {
public:
    Pixel(double size_, Vector3d &coordinate_, pair<unsigned long int,unsigned long int> &ij_);
    Vector3d getCoordinate();
    double getValue(string value);
    void setValue(string value, double newvalue);
    void scale_coordinates(double scale);

private:
    double size_;
    Vector3d coordinate_;
    pair<unsigned long int,unsigned  long int> ij_;
    map<string,double> values_;


};


#endif //JETSHOW_PIXEL_H
