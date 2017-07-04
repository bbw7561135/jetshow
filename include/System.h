#ifndef JETSHOW_SYSTEMS_H
#define JETSHOW_SYSTEMS_H

#include <iostream>
#include "Jet.h"


typedef std::vector<double> state_type;


class System {
  public:
    System(Jet* newjet, Vector3d &newpoint_in, Vector3d &newray_direction,
           double newnu);

  virtual void operator() (const double &x, double &dxdt, const double t) = 0;

	protected:
    Jet* jet;
    Vector3d point_in;
    Vector3d ray_direction;
    double nu;
};


class Tau : public System {
  public:
    Tau(Jet* newjet, Vector3d &newpoint_in, Vector3d &newray_direction,
        double newnu, double newtau_max) : System(newjet, newpoint_in,
                                               newray_direction, newnu) {};
    void operator() (const double &x, double &dxdt, const double t) override;
};


class I : public System {
 public:
  I(Jet* newjet, Vector3d &newpoint_in, Vector3d &newray_direction,
    double newnu) : System(newjet, newpoint_in, newray_direction, newnu) {};

  void operator() (const double &x, double &dxdt, const double t) override;
};


class FullStokes {
public:
		FullStokes(Jet* newjet, Vector3d &newpoint_in, Vector3d &newray_direction,
		double newnu);

		void operator() (const state_type &x, state_type &dxdt, const double t);

protected:
		Jet* jet;
		Vector3d point_in;
		Vector3d ray_direction;
		double nu;
};


// This debug-print integrands and ``time`` variable
void write_cout(const double &x, const double t);


bool check_opt_depth(double tau_max, const double &x);


bool check_opt_depth_xt(double tau_max, std::pair<double,double> xt);


#endif //JETSHOW_SYSTEMS_H
