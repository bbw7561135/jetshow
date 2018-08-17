#ifndef JETSHOW_OBSERVATION_H
#define JETSHOW_OBSERVATION_H

#include <string>
#include "Jet.h"
#include "ImagePlane.h"

using std::pair;
typedef std::vector<double> state_type;


class Observation {
  public:
    Observation(Jet* newjet, ImagePlane* imagePlane, double nu);
    void run(int n, double tau_max, double dt_max, double tau_min,
             string integration_type, string output_type,
             int n_max, double tau_n_min, double tau_n_max);
		void run_stripe(int n, double tau_max, double tau_min);
    vector<vector<double>> getImage(string value);
		vector<double> getStripe(string value);
    pair<int,int> getImageSize();
    const double nu;
  private:
    Jet* jet;
    ImagePlane* imagePlane;

		pair<double, double>
		integrate_tau(std::list<Intersection>& list_intersect,
		              Vector3d ray_direction, const double nu,
		              double tau_max, int n);

		pair<double, double>
		integrate_tau_adaptive(std::list<Intersection>& list_intersect,
		                       Vector3d ray_direction, const double nu,
		                       double tau_max, int n, double dt_max);

		void integrate_i(std::list<Intersection>& list_intersect,
		                 Vector3d ray_direction, const double nu,
		                 int n, double tau, double tau_n_min, double tau_n_max,
		                 double& background_I);

		void integrate_i_adaptive(std::list<Intersection>& list_intersect,
		                          Vector3d ray_direction, const double nu,
		                          int n, double tau, double& background_I);

//		void integrate_full_stokes_adaptive(std::list<Intersection>& list_intersect,
//		                                    Vector3d ray_direction, const double nu,
//		                                    int n, double tau,
//		                                    state_type& background);
};

#endif //JETSHOW_OBSERVATION_H
