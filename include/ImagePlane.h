#ifndef JETSHOW_IMAGEPLANE_H
#define JETSHOW_IMAGEPLANE_H

#include <array>
#include <Eigen/Eigen>
#include <memory>
#include "Image.h"
#include "Ray.h"

using std::array;
using std::pair;
using Eigen::Vector3d;


class ImagePlane {
public:
//    ImagePlane(const ImagePlane&) = delete;
//    ImagePlane& operator=(const ImagePlane&) = delete;
//    ~ImagePlane() = default;
//    explicit ImagePlane(pair<int,int> image_size, double pixel_size,
//                        double pixel_scale, double los_angle);
		ImagePlane(pair<int,int> image_size, double pixel_size, double pixel_scale,
							 double los_angle);
    vector<Pixel>& getPixels();
    vector<Ray>& getRays();
    vector<vector<double>> getImage(string value);
		vector<double> getStripe(string value);
    const pair<int,int> image_size;

private:
    Image image_;
    double los_angle_;
    Vector3d direction_;
    vector<Ray> rays_;
};


#endif //JETSHOW_IMAGEPLANE_H
