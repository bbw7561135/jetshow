//
// Created by ilya on 5/31/17.
//

#ifndef JETSHOW_INTERCEPTION_H
#define JETSHOW_INTERCEPTION_H

#include <Eigen/Eigen>
#include "Ray.h"
#include "Geometry.h"

using Eigen::Vector3d;
class Geometry;
class Ray;


class Intersection {
  public:
  // Ctor for finite intersections
  Intersection(const Ray &ray, const Vector3d &point_in,
               const Vector3d &point_out);
  // Ctor for full infinite intersections
  Intersection(const Ray &ray, const Geometry &geo);
  // Ctor for half infinite intersections
  Intersection(const Ray &ray, const Vector3d &point, const Geometry &geo);

    const Vector3d& direction() const;
    void set_direction(const Vector3d &direction);
    void set_borders(const Vector3d &point_in, const Vector3d &point_out);
    bool is_finite() const;
    const std::pair<Vector3d, Vector3d> & get_path() const;
  private:
    Vector3d direction_;
    bool is_finite_;
    std::pair<Vector3d, Vector3d> borders_;
};

#endif //JETSHOW_INTERCEPTION_H
