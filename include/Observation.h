#ifndef JETSHOW_OBSERVATION_H
#define JETSHOW_OBSERVATION_H

#include <string>
#include "Jet.h"
#include "ImagePlane.h"

using std::pair;


class Observation {
 public:
  Observation(Jet* newjet, ImagePlane* imagePlane, double nu);
  void run(int n, double tau_max, double dt_max, double tau_min);
  vector<vector<double>>& getImage(string value);
  pair<int,int> getImageSize();
  const double nu;
 private:
  Jet* jet;
  ImagePlane* imagePlane;

};

#endif //JETSHOW_OBSERVATION_H
