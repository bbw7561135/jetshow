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

// ``pixel_scale`` should be number of cm in one pixel. If one pixel is 0.1 mas
// and ``mas_to_pc(z)`` parsecs in one mas for given redshift ``z`` then it
// should be ``0.1*mas_to_pc*pc_to_cm`` where pc_to_cm ~ 3*10^18cm.
class Image {
public:
//    Image(const Image&) = delete;
//    Image& operator=(const Image&) = delete;
//    ~Image() = default;
//    explicit Image(pair<int, int> image_size, double pixel_size,
//                   double pixel_scale);
    Image(pair<unsigned long int, unsigned long int> image_size, double lg_pixel_size_start_cm, double lg_pixel_size_stop_cm);

    const unsigned long int num_of_pixels_;
    const pair<unsigned long int,unsigned long int> image_size;
    vector<vector<double>> getImage(string value);
    vector<vector<double>> getPixelSizes();
    vector<double> getStripe(string value);
    Vector3d getCoordinate(unsigned long int i, unsigned long int j);
    vector<Pixel>& getPixels();

private:
    pair<unsigned long int,unsigned long int> image_size_;
    // Log10 of pixel sizes in cm
    double pixel_size_start_;
    double pixel_size_stop_;
    // n_across, n_along array with pixel sizes [in cm!]
    vector<vector<double>> pixel_sizes_;
    // 2 arrays of the coordinates of pixel centers (n_across, n_along) [in cm!]
    vector<vector<double>> pixel_center_coordinates_along_;
    vector<vector<double>> pixel_center_coordinates_across_;
    vector<Pixel> pixels_;

};


#endif //JETSHOW_IMAGE_H
