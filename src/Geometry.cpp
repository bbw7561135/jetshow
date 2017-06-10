#include "Geometry.h"


std::pair<Vector3d, Vector3d> Geometry::full_infinite_path(Ray &ray) const {
  Vector3d diff = ray.origin() - origin();
  Vector3d ref_point = diff - diff.dot(ray.direction()) * ray.direction();
  Vector3d point_in = ref_point - big_scale() * ray.direction();
  Vector3d point_out = ref_point + big_scale() * ray.direction();
  return std::pair<Vector3d,Vector3d>{point_in, point_out};
}

std::pair<Vector3d, Vector3d> Geometry::half_infinite_path(Ray &ray,
                                                           const Vector3d &point) const {
  Vector3d check_point = point-ray.direction();
  if (is_within(check_point)) {
    return std::pair<Vector3d,Vector3d>{point-big_scale()*ray.direction(), point};

  } else {
    return std::pair<Vector3d,Vector3d>{point, point+big_scale()*ray.direction()};
  }
}
