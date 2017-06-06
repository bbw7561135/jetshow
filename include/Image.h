//
// Created by ilya on 6/5/17.
//

#ifndef JETSHOW_IMAGE_H
#define JETSHOW_IMAGE_H

#include <array>
#include <string>
#include "Pixel.h"

using std::array;

class Image {
public:
    Image(array<int,2> image_size, array<double,2> pixel_size, array<double,2> pixel_scales);
    const int num_of_pixels_;
    array<array<double,image_size_[0]>, image_size_[1]> get_image(string value);

private:
    array<int,2> image_size_;
    array<double,2> pixel_size_;
    array<double,2> pixel_scales_;
    array<Pixel,num_of_pixels_> pixels_;

};


#endif //JETSHOW_IMAGE_H
