#include "Intersection.h"


const Vector3d& Intersection::direction() const {
  return direction_;
}

void Intersection::set_direction(const Vector3d &direction) {
  direction_ = direction;
}


void Intersection::set_borders(const Vector3d &point_in,
                               const Vector3d &point_out) {
  borders_ = std::pair<Vector3d, Vector3d>{point_in, point_out};
}

const std::pair<Vector3d, Vector3d> & Intersection::get_path() const {
  return borders_;
}

// Ctor for finite intersections
Intersection::Intersection(const Ray &ray,
                           const Vector3d &new_point_in,
                           const Vector3d &new_point_out) {
  set_direction(ray.direction());
  set_borders(new_point_in, new_point_out);
  is_finite_ = 1;
}

// Ctor for full infinite intersections
Intersection::Intersection(const Ray &ray, const Geometry &geo) {
  Vector3d diff = ray.origin() - geo.origin();
  Vector3d ref_point = diff - diff.dot(ray.direction()) * ray.direction();
  Intersection(ray,
               ref_point - geo.big_scale() * ray.direction(),
               ref_point + geo.big_scale() * ray.direction());
  is_finite_ = 0;
}

// Ctor for half infinite intersections
Intersection::Intersection(const Ray &ray, const Vector3d &point,
                           const Geometry &geometry) {
  // Check if ``ray`` was inside ``geometry`` before coming to ``point``.
  Vector3d check_point = point-ray.direction();
  if (geometry.is_within(check_point)) {
      Intersection(ray, point-geometry.big_scale()*ray.direction(),
                   point);
  } else {
    Intersection(ray, point, point+geometry.big_scale()*ray.direction());
  }
}


bool Intersection::is_finite() const {
  return is_finite_;
}
