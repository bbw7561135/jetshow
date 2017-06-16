#include <string>
#include <iostream>
#include "Observation.h"
#include "System.h"
#include <boost/range/algorithm.hpp>
#include <boost/numeric/odeint.hpp>

using namespace boost::numeric::odeint;


Observation::Observation(Jet *newjet, ImagePlane *newimagePlane, double newnu) :
    nu(newnu)
{
  jet = newjet;
  imagePlane = newimagePlane;
};


void Observation::run(int n, double tau_max, double dtau_max) {
  auto image_size = getImageSize();
	vector<std::unique_ptr<Pixel>> pixels = std::move(imagePlane->getPixels());
	vector<std::unique_ptr<Ray>> rays = std::move(imagePlane->getRays());
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
      std::list<Intersection> list_intersect = jet->hit(*ray);
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
          auto borders = (*it).get_path();
          Vector3d point_in = borders.first;
          Vector3d point_out = borders.second;
          double length = (point_out - point_in).norm();
          // FIXME: cast to double?
          double dt = length/n;

          // First integrate till some ``tau_max`` using s
          Tau tau(jet, point_in, ray_direction, nu, tau_max);
          double optDepth = background_tau;
          typedef runge_kutta_dopri5< double > stepper_type;
          auto stepper = make_dense_output(1E-21, 1E-18, dtau_max,
                                           stepper_type());
					using namespace std::placeholders;
					auto is_done = std::bind(check_opt_depth, tau_max, _1);
					auto ode_range = make_adaptive_range(std::ref(stepper), tau, optDepth,
																							 0.0, length, dt);
					auto found_iter = std::find_if(ode_range.first, ode_range.second,
																				 is_done);
					if (found_iter == ode_range.second) {
						// This should prevale for optically thin parts
						std::cout << "Tau less then tau_max" << std::endl;
					} else {
						double t_tau_max = stepper.current_time();
						Ray ray_tau_max(point_in, ray_direction);
						Vector3d point_out_tau_max = ray_tau_max.point(t_tau_max);
						it.operator*().set_point_out(point_out_tau_max);
						// Delete all other furtherst  intersections if any
						// FIXME: Should i check if any other intersections left?
						it = list_intersect.erase(it++, list_intersect.end());
					}
//					std::cout << "Length " << length << std::endl;
//					std::cout << found_iter.get_state() << std::endl;
					// This works only for make_dense_output()
//					std::cout << stepper.current_time() << std::endl;
          // TODO: If resulting ``optDepth`` > ``tau_max`` then find point where
          // it equals and set ``Intersection`` object's ``point_out`` to this
          // point. Also if ``list_intersect`` contains other elements - delete
          // them using ``erase`` method (``list.erase(++it, list.end)``).


          // Update background value (or write final values if this is last
          // cycle)
          background_tau += found_iter.get_state();
        }

        // Write final values here inside for-cycle
        double background_I = 0.;
        for (auto it = list_intersect.rend();
             it != list_intersect.rbegin(); --it) {
					auto borders = (*it).get_path();
					Vector3d point_in = borders.first;
					Vector3d point_out = borders.second;
					double length = (point_out - point_in).norm();
					// FIXME: cast to double?
					double dt = length/n;

          Vector3d inv_direction = -1.*ray_direction;
          I stokesI(jet, point_out, inv_direction, nu);
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
				std::string value ("tau");
        pxl->setValue(value, background_tau);
				value = "I";
				pxl->setValue(value, background_I);
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
