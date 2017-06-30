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
void Observation::run(int n, double tau_max, double dt_max, double tau_min,
                      string type, int n_max, double tau_n_min,
                      double tau_n_max) {
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
        continue;
      } else {

	      std::pair<double,double> tau_l_end;
	      if (type == "constant") {
	        tau_l_end = integrate_tau(list_intersect, ray_direction, nu, tau_max,
	                                  n);
	      } else if (type == "adaptive") {
		      tau_l_end = integrate_tau_adaptive(list_intersect, ray_direction, nu,
		                                         tau_max, n, dt_max);
	      }
	      double background_tau = tau_l_end.first;
	      double thickness = tau_l_end.second;

	      // Write final values here inside integrate_i
				double background_I = 0.;
				// Calculate I only if optical depth is high enough
				if (background_tau > tau_min) {
					string local_type;
					// In high optical depth pixels use adaptive steps
					if (background_tau > 0.01*tau_max) {
						local_type = "adaptive";
					} else {
						local_type = "constant";
					}
					// If adaptive everywhere => then adaptive here too
					if (type == "adaptive") {
						local_type = "adaptive";
					}
					if (local_type == "constant") {
						integrate_i(list_intersect, ray_direction, nu, n, background_tau,
					              tau_n_min, tau_n_max, background_I);
					} else if (local_type == "adaptive") {
						integrate_i_adaptive(list_intersect, ray_direction, nu, n,
						                     background_tau, background_I);
					}
				} else {
//					std::cout << "Too small optical depth..." << std::endl;
				}
        // Write values to pixel
				std::string value ("tau");
        pxl.setValue(value, background_tau);
				value = "I";
				pxl.setValue(value, background_I);
				value = "l";
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

pair<double, double>
Observation::integrate_tau(std::list<Intersection>& list_intersect,
                           Vector3d ray_direction, const double nu,
                           double tau_max, int n) {
	// Write final values in this variables inside for-cycle
	double background_tau = 0.;
	double thickness = 0.;

	for (auto it = list_intersect.begin();
	     it != list_intersect.end(); ++it) {
		auto borders = (*it).get_path();

		Vector3d point_in = borders.first;
		Vector3d point_out = borders.second;

		double length = (point_out - point_in).norm();
		double dt = length/n;

		// Directed to observer (as we integrate ``k`` only)
		Vector3d ray_direction_ = -1. * ray_direction;

		// First integrate till some ``tau_max``
		Tau tau(jet, point_in, ray_direction_, nu, tau_max);
		// This is out State
		double optDepth = 0.0;
		typedef runge_kutta4< double > stepper_type;
		auto stepper = stepper_type();
		using namespace std::placeholders;
		auto is_done = std::bind(check_opt_depth_xt, tau_max, _1);
		auto ode_range = make_const_step_time_range(stepper, tau, optDepth,
		                                            0.0, length, dt);
		auto found_iter = std::find_if(ode_range.first, ode_range.second,
		                               is_done);
		if (found_iter != ode_range.second) {
			auto found_state = found_iter.operator*();
			double t_tau_max = found_state.second;

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
	return std::make_pair(background_tau, thickness);
}


pair<double, double>
Observation::integrate_tau_adaptive(std::list<Intersection> &list_intersect,
                                    Vector3d ray_direction, const double nu,
                                    double tau_max, int n, double dt_max) {

	// Write final values in this variables inside for-cycle
	double background_tau = 0.;
	double thickness = 0.;

	for (auto it = list_intersect.begin();
	     it != list_intersect.end(); ++it) {
		std::pair<double,double> tau_t_end;
		auto borders = (*it).get_path();

		Vector3d point_in = borders.first;
		Vector3d point_out = borders.second;


		double length = (point_out - point_in).norm();
		double dt = length/n;

		// Directed to observer (as we integrate ``k`` only)
		Vector3d ray_direction_ = -1. * ray_direction;

		// First integrate till some ``tau_max``
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
		if (found_iter != ode_range.second) {
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
	return std::make_pair(background_tau, thickness);
}


void Observation::integrate_i(std::list<Intersection> &list_intersect,
                              Vector3d ray_direction, const double nu,
                              int n, double background_tau, double tau_n_min,
                              double tau_n_max, double& background_I) {

	for (auto it = list_intersect.rbegin();
	     it != list_intersect.rend(); ++it) {
		auto borders = (*it).get_path();
		Vector3d point_in = borders.first;
		Vector3d point_out = borders.second;

		double length = (point_out - point_in).norm();
		auto auto_n = steps_schedule(background_tau, n, 10*n, tau_n_min, tau_n_max);
		double dt = length / auto_n;

		Vector3d inv_direction = -1. * ray_direction;
		I stokesI(jet, point_out, inv_direction, nu);
		typedef runge_kutta4< double > stepper_type;
		auto stepper = stepper_type();
		double stI = background_I;
		integrate_const(stepper, stokesI, stI, 0.0, length, dt);
		background_I = stI;
	}
}

void Observation::integrate_i_adaptive(std::list<Intersection> &list_intersect,
                                       Vector3d ray_direction, const double nu,
                                       int n, double background_tau,
                                       double& background_I) {

	for (auto it = list_intersect.rbegin();
	     it != list_intersect.rend(); ++it) {
		auto borders = (*it).get_path();
		Vector3d point_in = borders.first;
		Vector3d point_out = borders.second;

		double length = (point_out - point_in).norm();
		double auto_n = n;
		if (background_tau > 0.1) {
			auto_n = steps_schedule(background_tau, n, 10*n);
		}
		double dt = length / auto_n;

		Vector3d inv_direction = -1. * ray_direction;
		I stokesI(jet, point_out, inv_direction, nu);
		typedef runge_kutta_dopri5<double> stepper_type;
		auto stepper = stepper_type();

		double stI = background_I;
		// TODO: Add ``dt_max`` constrains (using ``make_dense_output``)
		// One can add observer function at the end of the argument list.
		integrate_adaptive(make_controlled(1E-9, 1E-9, stepper_type()),
		                   stokesI,
		                   stI, 0.0, length, dt);
		background_I = stI;
	}
}
