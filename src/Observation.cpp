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

// ``dt_max`` - max. step in pc.
void Observation::run(int n, double tau_max, double dt_max, double tau_min) {
	dt_max *= pc;
  auto image_size = getImageSize();
	vector<Pixel>& pixels = imagePlane->getPixels();
	vector<Ray>& rays = imagePlane->getRays();
	// Comment out for easy debug printing
	#pragma omp parallel for schedule(dynamic) collapse(2)
	// Actually, the best user-time:
	// #pragma omp parallel for num_threads(4) schedule(dynamic) collapse(2)
  for (int j = 0; j < image_size.first; ++j) {
    for (int k = 0; k < image_size.second; ++k) {
      int n_pix = image_size.first*j + k + 1;
      std::cout << "Running on pixel # " << n_pix << std::endl;
      auto &ray = rays[j*image_size.first+k];
      auto &pxl = pixels[j*image_size.first+k];

			auto ray_direction = ray.direction();
      std::list<Intersection> list_intersect = jet->hit(ray);
      if (list_intersect.empty()) {
//        std::cout << "No intersections" << std::endl;
        continue;
      } else {
        // Do transfer here
        // Write final values in this variables inside for-cycle
        double background_tau = 0.;
				double thickness = 0.;

        // On some ``it`` we can stop because of tau > tau_max. Then we should
        // go from that ``it`` back.
        for (auto it = list_intersect.begin();
             it != list_intersect.end(); ++it) {
          auto borders = (*it).get_path();
          Vector3d point_in = borders.first;
          Vector3d point_out = borders.second;


	        double length = (point_out - point_in).norm();
          double dt = length/n;

          // First integrate till some ``tau_max``
	        // Directed to observer (as we integrate ``k`` only)
					Vector3d ray_direction_ = -1. * ray_direction;
					Tau tau(jet, point_in, ray_direction_, nu, tau_max);
					// This is out State
	        double optDepth = 0.0;
          typedef runge_kutta_dopri5< double > stepper_type;
          auto stepper = make_dense_output(1E-3, 1E-14, dt_max,
                                           stepper_type());
					using namespace std::placeholders;
					auto is_done = std::bind(check_opt_depth, tau_max, _1);
					auto ode_range = make_adaptive_range(std::ref(stepper), tau, optDepth,
																							 0.0, length, dt);
					auto found_iter = std::find_if(ode_range.first, ode_range.second,
																				 is_done);

	        // If tau>tau_max for some iteration -> get point where it occurs
					if (found_iter != ode_range.second) {
						// This is t[pc] where tau = ``tau_max11
//						std::cout << "Tau larger then tau_max" << std::endl;
						double t_tau_max = stepper.current_time();

						// TODO: Use https://github.com/headmyshoulder/odeint-v2/blob/master/examples/find_crossing.cpp
						// to find more precise value of t(tau_max)
//						double t_tau_previous = stepper.previous_time();

						Ray ray_tau_max(point_in, ray_direction);
						Vector3d point_out_tau_max = ray_tau_max.point(t_tau_max);
						// Set point of tau=tau_max as outer border for further integrations
						it.operator*().set_point_out(point_out_tau_max);

						// Delete all other furthest intersections if any
						if (list_intersect.size() > 1) {
							++it;
							// This now is container.end()
							it = list_intersect.erase(it, list_intersect.end());
							// Move it before container.end() to make it end() on the new
							// cycle
							--it;
						}
					}
          // Update background value (or write final values if this is last
          // cycle)
          background_tau += found_iter.get_state();
	        // FIXME: Here ``optDepth`` is 0 because out of scope or what? Not
	        // because of scope as ``length`` is nonzero.
					thickness += length;
        }
				double background_I = 0.;
				// Calculate I only if optical depth is high enough
				if (background_tau > tau_min) {
					background_I= 0.;
					// Write final values here inside for-cycle
					for (auto it = list_intersect.rbegin();
							 it != list_intersect.rend(); ++it) {
						auto borders = (*it).get_path();
						Vector3d point_in = borders.first;
						Vector3d point_out = borders.second;

						double length = (point_out - point_in).norm();
						double dt = length / n;

						Vector3d inv_direction = -1. * ray_direction;
//						std::cout << " Direction for integrating I " << inv_direction << std::endl;
						I stokesI(jet, point_out, inv_direction, nu);
						typedef runge_kutta_dopri5<double> stepper_type;

						// Here stI - Stokes I intensity.
						double stI = background_I;
						// TODO: Add ``dt_max`` constrains
//						integrate_adaptive(make_controlled(1E-9, 1E-9, stepper_type()),
//															 stokesI,
//															 stI, 0.0, length, dt, write_cout);
						integrate_adaptive(make_controlled(1E-9, 1E-9, stepper_type()),
						                   stokesI,
						                   stI, 0.0, length, dt);
						// Update background value (or write final values if this is last
						// cycle)
						background_I = stI;
					}
				} else {
//					std::cout << "Too small optical depth..." << std::endl;
				}
        // Write values to pixel
				std::string value ("tau");
//	      std::cout << "Writing to pixel value of tau = " << background_tau << std::endl;
        pxl.setValue(value, background_tau);
				value = "I";
				pxl.setValue(value, background_I);
				value = "l";
//	      std::cout << "Writing to pixel value of length = " << thickness << std::endl;
				pxl.setValue(value, thickness);
      }
    }
  }
}


vector<vector<double>> Observation::getImage(string value) {
  return imagePlane->getImage(value);
}


pair<int, int> Observation::getImageSize() {
  return imagePlane->image_size;
}
