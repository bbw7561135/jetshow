#include "Observation.h"
#include <boost/range/algorithm.hpp>


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
  // Cycle for each pixel+ray and make transfer for it. THIS CYCLE HAS TO BE
  // PARALLELIZED!
  for (int j = 0; j < image_size.first; ++j) {
    for (int k = 0; k < image_size.second; ++k) {
      // This is debug printing out begins
      int n_pix = image_size.first*j + k + 1;
      std::cout << "Running on pixel # " << n_pix << std::endl;
      // This is debug printing out ends
      auto ray = rays[j*image_size.first+k].get();
      auto pxl = pixels[j*image_size.first+k].get();
      auto ray_direction = ray->direction();
      std::list<Intersection> list_intersect = jet.hit(&ray);
      if (list_intersect.empty()) {
        // This is debug printing out begins
        std::cout << "No intersections" << std::endl;
        // This is debug printing out ends
        continue;
      } else {
        // Do transfer here
        // Write final values here inside for-cycle
        double background_tau = 0.;

        // On some ``it`` we can stop because of tau > tau_max. Then we should
        // go from that ``it`` back. How to implement it using iterators?
        for (auto it = list_intersect.begin();
             it != list_intersect.end(); ++it) {
          auto borders = it.get_path();
          Vector3d point_in = borders.first;
          Vector3d point_out = borders.second;
          double length = (point_out - point_in).norm();
          // FIXME: cast to double?
          double dt = length/n;

          // First integrate till some ``tau_max`` using s
          Tau tau(jet, point_in, ray_direction, nu);
          double optDepth = background_tau;
          typedef runge_kutta_dopri5< double > stepper_type;
          auto stepper = make_controlled(1E-21, 1E-18, dtau_max,
                                         stepper_type());
          auto iter = boost::find_if(make_adaptive_range(stepper, tau, optDepth,
                                                         0.0 , length, dt));
          std::cout<<"integration stopped at"<<optDepth<<std::endl;
//          integrate_adaptive(stepper, tau, optDepth, 0.0 , length, dt,
//                             write_cout);

          // TODO: If resulting ``optDepth`` > ``tau_max`` then find point where
          // it equals and set ``Intersection`` object's ``point_out`` to this
          // point. Also if ``list_intersect`` contains other elements - delete
          // them using ``erase`` method (``list.erase(++it, list.end)``).


          // Update background value (or write final values if this is last
          // cycle)
          background_tau = optDepth;
          // End for-cycle if optical depth is too high
          if (optDepth >= tau_max) {
            list_intersect.erase(++it, list_intersect.end());
            break;
          }
        }

        // Write final values here inside for-cycle
        double background_I = 0.;
        for (auto it = list_intersect.rend();
             it != list_intersect.rbegin(); --it) {

          Vector3d inv_direction = -1.*ray_direction;
          I stokesI(&bkjet, point_out, inv_direction, nu);
          typedef runge_kutta_dopri5< double > stepper_type;

          // Here stI - Stokes I intensity.
          double stI = background_I;
          integrate_adaptive(make_controlled(1E-12, 1E-12, stepper_type()), stokesI,
                             stI, 0.0 , length, dt, write_cout);
          // Update background value (or write final values if this is last
          // cycle)
          background_I = stI;

        }
        // Write values to pixel
        pxl->setValue('tau', background_tau);
        pxl->setValue('I', background_I);
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
