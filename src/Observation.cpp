#include "Observation.h"


Observation::Observation(Jet *newjet, ImagePlane *newimagePlane, double newnu) :
    nu(newnu)
{
  jet = newjet;
  imagePlane = newimagePlane;
};


void Observation::run(int n, double tau_max, double dtau_max) {
  auto image_size = getImageSize();
  auto pixels = imagePlane->getPixels();
  auto rays = imagePlane->getRays();
  for (int j = 0; j < image_size.first; ++j) {
    for (int k = 0; k < image_size.second; ++k) {
      auto ray = rays[j*image_size.first+k].get();
      auto pxl = pixels[j*image_size.first+k].get();
      auto ray_direction = ray->direction();
      std::list<Intersection> list_intersect = jet.hit(&ray);
      if (list_intersect.empty()) { continue;
      } else {
        // Do transfer here
        double background = 0.;
        for (auto it = list_intersect.begin();
             it != list_intersect.end(); ++it) {
          auto borders = intersect.get_path();
          Vector3d point_in = borders.first;
          Vector3d point_out = borders.second;
          double length = (point_out - point_in).norm();
          // FIXME: double cast
          double dt = length/n;

          // This is integration part
          // First integrate till some ``tau_max`` using s
          Tau tau(jet, point_in, ray_direction, nu);
          typedef runge_kutta_dopri5< double > stepper_type;
          // Here x - optical depth \tau
          double x = 0.0;
          integrate_adaptive(make_controlled(1E-21, 1E-18, dtau_max, stepper_type()),
                             tau, x, 0.0 , length, dt, write_cout);

          Vector3d inv_direction = -1.*ray_direction;
          I stokesI(&bkjet, point_out, inv_direction, nu);
          // Here x - Stokes I intensity.
          typedef runge_kutta_dopri5< double > stepper_type;

          double stI = 0.0;
          integrate_adaptive(make_controlled(1E-12, 1E-12, stepper_type()), stokesI,
                             stI, 0.0 , length, dt, write_cout);

          // Write values to pixel
          pxl->setValue('tau', x);
          pxl->setValue('I', stI);

        }
      }
    }
  }
}


vector<vector<double>> &Observation::getImage(string value) {
  return imagePlane->getImage(value);
}


pair<int, int> Observation::getImageSize() {
  return imagePlane->image_size;
}
