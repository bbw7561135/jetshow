#ifndef JETSHOW_SYSTEMS_H
#define JETSHOW_SYSTEMS_H

#include <iostream>
#include "Jet.h"


class System {
 public:
  // Default ctor
  // System();
    System(Jet* newjet, Vector3d &newpoint_in, Vector3d &newray_direction,
           double newnu);

  // This is what differs
  virtual void operator() (const double &x, double &dxdt, const double t) = 0;

  void setPointIn(Vector3d &newpoint_in) {
    point_in = newpoint_in;
  }
  void setRayDirection(Vector3d &newray_direction) {
    ray_direction = newray_direction;
  }

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


//class OptDepthObserver {
//  public:
//    OptDepthObserver(double tau_max, double &t);
//    void operator()( const double &x , double t ) const;
//
//  private:
//    double time_;
//    double value_;
//};


// This debug-print integrands and ``time`` variable
void write_cout(const double &x, const double t);


bool check_opt_depth(double tau_max, const double &x);
//
//class OptDepthMax {
//public:
//    OptDepthMax(double tau_max);
//    bool is_done(const double &x);
//
//private:
//    double tau_max;
//};

#endif //JETSHOW_SYSTEMS_H
