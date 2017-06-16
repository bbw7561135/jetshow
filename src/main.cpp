#include <iostream>
#include <boost/numeric/odeint.hpp>
#include "ImagePlane.h"
#include "Observation.h"
#include "Image.h"
#include "System.h"
#include "BField.h"
#include "VField.h"
#include "Ray.h"
#include "Cone.h"
#include "Jet.h"
#include "utils.h"
#include "math.h"
#include "linspace.h"
#include "Cell.h"
#include "NField.h"
#include "Cylinder.h"
#include "Pixel.h"
#include <boost/range/algorithm.hpp>
#include <memory>
#include <mpi.h>
#include <ctime>
#include <chrono>





using Eigen::Vector3d;
using std::vector;
using std::pair;
using namespace boost::numeric::odeint;
typedef std::chrono::high_resolution_clock Clock;


void test_pixel() {
  auto ij = std::make_pair(1, 10);
  Vector3d coordinate{1., 2., 3.};
  double size = 0.1;
  Pixel pxl(0.1, coordinate, ij);
  std::cout << "Pixel coordinate = " << pxl.getCoordinate() << std::endl;
  pxl.scale_coordinates(100.);
  std::cout << "New Pixel coordinate = " << pxl.getCoordinate() << std::endl;
}

void test_image() {
  auto image_size = std::make_pair(10, 10);
  double pixel_size = 0.1;
  double pixel_scale = 100.;
  Image image(image_size, pixel_size, pixel_scale);
  vector<Pixel>& pixels = image.getPixels();
  for (int j = 0; j < image_size.first; ++j) {
    std::cout << "j = " << j << std::endl;
    for (int k = 0; k < image_size.second; ++k) {
      std::cout << "k = " << k << std::endl;
      auto pxl = pixels[j*image_size.first+k];
      std::cout << "coordinate " << pxl.getCoordinate() << std::endl;
    }
  }

	vector<Pixel>& pixels2 = image.getPixels();
	for (int j = 0; j < image_size.first; ++j) {
		std::cout << "j = " << j << std::endl;
		for (int k = 0; k < image_size.second; ++k) {
			std::cout << "k = " << k << std::endl;
			auto pxl = pixels2[j*image_size.first+k];
			std::cout << "coordinate2 " << pxl.getCoordinate() << std::endl;
		}
	}
}

void test_image_plane() {
  auto image_size = std::make_pair(10, 10);
  double pixel_size = 0.1;
  double pixel_scale = 100.;
  double los_angle = pi/6.;
  ImagePlane imagePlane(image_size, pixel_size, pixel_scale, los_angle);
  vector<Ray>& rays = imagePlane.getRays();
  for (int j = 0; j < image_size.first; ++j) {
    std::cout << "j = " << j << std::endl;
    for (int k = 0; k < image_size.second; ++k) {
      std::cout << "k = " << k << std::endl;
      auto ray = rays[j*image_size.first+k];
      std::cout << "Ray origin " << ray.origin() << std::endl;
      std::cout << "Ray direction " << ray.direction() << std::endl;
    }
  }

	vector<Pixel>& pxls = imagePlane.getPixels();
	for (int j = 0; j < image_size.first; ++j) {
		std::cout << "j = " << j << std::endl;
		for (int k = 0; k < image_size.second; ++k) {
			std::cout << "k = " << k << std::endl;
			auto pxl = pxls[j*image_size.first+k];
			std::cout << "Pixel Coordinate " << pxl.getCoordinate() << std::endl;
		}
	}
}


void test_jet() {
    // Create cone geometry of jet
    Vector3d origin = {0., 0., 0.};
    Vector3d direction = {0., 0., 1.};
    double angle = pi/4.;
    double scale = 10.;
    Cone geometry(origin, direction, angle, scale);

    RadialConicalBField bfield(0.1, 2.0);
    ConstCenralVField vfield(10.*c);
    BKNField nfield(1., 10.);

    Jet bkjet(&geometry, &vfield, &bfield, &nfield);
    // Test just one point
//    Vector3d test_point{0., 0.3, 3.};
//    std::cout << bkjet.getV(test_point) << std::endl;
//    std::cout << bkjet.getN(test_point) << std::endl;
//    std::cout << bkjet.getB(test_point) << std::endl;


    // Get a ray
    Vector3d ray_direction(0., 1., 0.);
    Vector3d ray_origin(0., 1., 1.);
    Ray ray(ray_origin, ray_direction);
    std::list<Intersection> list_intersect = geometry.hit(ray);
    Intersection intersect = list_intersect.front();
    std::pair<Vector3d,Vector3d> borders = intersect.get_path();
    Vector3d point_in = borders.first;
    Vector3d point_out = borders.second;
    double length = (point_out - point_in).norm();
    double dt = length/100.;
    std::cout << "Dt " << dt << std::endl;
    double nu = 5.*pow(10., 9.);
    double tau_max = 1E-18;
    double dt_max = 0.01;

    // This is integration part
    Tau tau(&bkjet, point_in, ray_direction, nu, tau_max);
//    OptDepthMax stopTau(tau_max);
    typedef runge_kutta_dopri5< double > stepper_type;
    // Here x - optical depth \tau
    double x = 0.0;
//    std::unique_ptr<double>


    // New part
//    auto stepper = make_controlled(1E-21, 1E-18, dt_max,
//                                   stepper_type());
  // This stepper doesn't hold value in x (use iter.get_state())
    auto stepper = make_dense_output(1E-21, 1E-18, dt_max,
                                     stepper_type());
    using namespace std::placeholders;
    auto is_done = std::bind(check_opt_depth, tau_max, _1);
    auto ode_range = make_adaptive_range(std::ref(stepper), tau, x, 0.0, length,
                                       dt);
    auto iter = std::find_if(ode_range.first, ode_range.second, is_done);

    if (iter == ode_range.second) {
      std::cout << "Tau less then tau_max" << std::endl;
    // no threshold crossing -> return time after t_end and ic
    }
    std::cout << "Length " << length << std::endl;
    std::cout << iter.get_state() << std::endl;
    // This works only for make_dense_output()
    std::cout << stepper.current_time() << std::endl;



// Using instead of this:
//    integrate_adaptive(make_controlled(1E-21, 1E-18, 0.01, stepper_type()), tau, x,
//                       0.0 , length, dt, write_cout);


//    Vector3d inv_direction = -1.*ray_direction;
//    I stokesI(&bkjet, point_out, inv_direction, nu);
//     Here x - Stokes I intensity.
//    typedef runge_kutta_dopri5< double > stepper_type;
//
//    double stI = 0.0;
//    integrate_adaptive(make_controlled(1E-12, 1E-12, stepper_type()), stokesI,
//                       stI, 0.0 , length, dt, write_cout);

}


void test_observations() {
    Vector3d origin = {0., 0., 0.};
    Vector3d direction = {0., 0., 1.};
    double angle = pi/6.;
    double scale = 10.;
    Cone geometry(origin, direction, angle, scale);

    RadialConicalBField bfield(0.1, 2.0);
    ConstCenralVField vfield(10.*c);
    BKNField nfield(1., 10.);

    Jet bkjet(&geometry, &vfield, &bfield, &nfield);

    auto image_size = std::make_pair(500, 500);
    double pixel_size = 0.1;
    double pixel_scale = 100.;
    double los_angle = pi/3.;
    ImagePlane imagePlane(image_size, pixel_size, pixel_scale, los_angle);

    double nu = 5.*pow(10., 9.);
    Observation observation(&bkjet, &imagePlane, nu);
    double tau_max = 100.;
    double dt_max = 10.;
    int n = 100;
    observation.run(n, tau_max, dt_max);
}

void test_erase() {
  std::list<double> a{10., 11., 12., 13., 14.};

  std::cout << "Before";
  for (auto it = a.begin(); it != a.end(); it++) {
    std::cout << " " << *it;}
  std::cout << std::endl;

  for (auto it = a.begin(); it != a.end(); ++it) {
    if (*it > 11.5) {
      it = a.erase(it++, a.end());
    }
  }

  std::cout << "After";
  for (auto it = a.begin(); it != a.end(); it++) {
    std::cout << " " << *it;}
  std::cout << std::endl;
  }

void test_mpi() {
		const int size = 12;
		double sinTable[size];

		#pragma omp parallel
		{
				// Code inside this region runs in parallel.
				printf("Hello!\n");
		}

		#pragma omp parallel for
		for(int n=0; n<size; ++n) {
			std::cout <<  n << std::endl;
			sinTable[n] = std::sin(2 * M_PI * n / size);
		}
};

int main() {
	auto t1 = Clock::now();
	std::clock_t start;
	start = std::clock();
//	test_mpi();
  test_observations();
//  test_erase();
//    test_jet();
//  test_pixel();
//    test_image();
//    test_image_plane();
	std::cout << "Time: " << (std::clock() - start) / (double)(CLOCKS_PER_SEC / 1000) << " ms" << std::endl;
	auto t2 = Clock::now();
	std::cout << "Delta t2-t1: "
						<< std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count()
						<< " milliseconds" << std::endl;

}

//  Intersection intersection{};
//  std::cout << intersection.direction() << std::endl;
//  auto borders = intersection.get_path();
//  std::cout << intersection.direction() << std::endl;
//  std::cout << borders.first << std::endl;
//    Vector3d cone_origin(0., 0., 0.);
//    Vector3d cone_direction(0., 0., 1.);
//    double cone_angle = pi/4.;
//    double scale = 10.0;
//    Cone cone = Cone(cone_origin, cone_direction, cone_angle, scale);
  // There's two intersections
//    Vector3d ray_direction(0., 1., 0.);
//    Vector3d ray_origin(0., 1., 1.);
  // No intersections
//    Vector3d ray_direction(0., 1., 0.);
//    Vector3d ray_origin(1., 1., 0.);
  // Two intersection (one far away)
//    Vector3d ray_direction(0., 1., 1.);
//    Vector3d ray_origin(0., 1., 0.);
  // One intersection
//    Vector3d ray_direction(0., 1., 0.);
//    Vector3d ray_origin(0., 1., 0.);
    // Along border
//    Vector3d ray_origin(0., 0., 0.);
//    Vector3d ray_direction(0., 1., 1.);
//
//    Ray ray(ray_origin, ray_direction);
//    std::list<Intersection> list_intersect = cone.hit(ray);
//    std::cout << "Did ray traverse volume ?" << std::endl << !list_intersect.empty() << std::endl;
//    if (!list_intersect.empty()) {
//        Intersection intersect = list_intersect.front();
//        std::pair<Vector3d,Vector3d> borders = intersect.get_path();
//        std::cout << "Borders : " << borders.first << " and " << borders.second << std::endl;
//    }


      // Test ``is_within``
//    Vector3d point{-2., 1., -5.};
//    bool is_within = cone.is_within(point);
//    std::cout << "Is within : " << is_within << std::endl;

    // This tests cylinder-ray intersections
//    Vector3d cylinder_origin(0., 0., 0.);
//    Vector3d cylinder_direction(0., 0., 1.);
//    double cylinder_radius = 1.;
//    Cylinder cylinder(cylinder_origin, cylinder_direction, cylinder_radius);
//    Vector3d point{0., 0., 2.};
//    bool is_within = cylinder.is_within(point);
//    std::cout << "Is within : " << is_within << std::endl;

    // No intersections
//    Vector3d ray_direction(1., 0., 1.);
//    Vector3d ray_origin(0., 2., 0.);

    // Along border
//    Vector3d ray_direction(0., 0., 1.);
//    Vector3d ray_origin(0., 1., 0.);

    // No interception. Along border internally
//    Vector3d ray_direction(0., 0., 1.);
//    Vector3d ray_origin(-0.1, 0.2, 0.);

    // No interception. Along border externally
//    Vector3d ray_direction(0., 0., 1.);
//    Vector3d ray_origin(10., -0.2, 0.);

    // Single interception
//    Vector3d ray_direction(1., 0., 0.);
//    Vector3d ray_origin(0., 1., 0.);

//    // Two intersections
//    Vector3d ray_direction(0., 1., 0.);
//    Vector3d ray_origin(0., 1., 0.);
//
//    Ray ray(ray_origin, ray_direction);
//    std::list<Intersection> list_intersect = cylinder.hit(ray);
//    std::cout << "Did ray traverse volume ?" << std::endl << !list_intersect.empty() << std::endl;
//    if (!list_intersect.empty()) {
//        Intersection intersect = list_intersect.front();
//        std::pair<Vector3d,Vector3d> borders = intersect.get_path();
//        std::cout << "Borders : " << borders.first << " and " << borders.second << std::endl;
//        std::cout << "Direction " << intersect.direction() << std::endl;
//    }
    // End of cylinder-ray intersections



//    Matrix3d eye_matrix;
//    eye_matrix << 1, 0, 0,
//                  0, 1, 0,
//                  0, 0, 1;
//    // DP
//    Vector3d delta = ray_origin - cone_origin;
//    std::cout << delta << std::endl;
//    // M
//    Matrix3d M = cone_direction * cone_direction.transpose() - cos(cone_angle)*cos(cone_angle)*eye_matrix;
//    std::cout << M << std::endl;
//    double c2 = ray_direction.transpose() * M * ray_direction;
//    double c1 = ray_direction.transpose() * M * delta;
//    double c0 = delta.transpose() * M * delta;
//    std::cout << "c2 = " << c2 << " c1 = " << c1 << " c0 = " << c0 << std::endl;
//    double d = c1*c1 - c0*c2;
//    std::cout << "d = " << d << std::endl;
//    std::cout << "-c1-sqrt(d) = " << -c1 - sqrt(d) << std::endl;
//    double t1 = (-c1 + sqrt(d)) / (c2);
//    double t2 = (-c1 - sqrt(d)) / (c2);
//    std::cout << "t1= " << t1 << " t2= " << t2 << std::endl;


//int main() {
//    Vector3d v1(0, 0, 0);
//    Vector3d v2(10, 10, -10);
//    std::vector<Vector3d> vec_1 = linspace(v1, v2, 11);
//    print_vector(vec_1);
//    return 0;
//}

//int main() {
//    Vector3d v1(0, 0, 0);
//    Vector3d v2(10, 10, 10);
//    std::vector<Vector3d> v = {v1, v2};
//    Cell cell = Cell(v);
//    std::vector<Cell> cells = cell.split_n(11);
//    for (size_t i=0; i != cells.size();++i) {
//        std::cout << "Cell number" << i << std::endl;
//        Print(cells[i].points());
//    }
//    return 0;
//
//}



//int main() {
//    Vector3d origin = {0., 0., 0.};
//    Vector3d direction = {1., 1., 1.};
//    direction.normalize();
//    Ray ray(origin, direction);
//    Vector3d point = ray.point(10);
//    std::cout << "Origin:" << std::endl;
//    std::cout << ray.origin() << std::endl;
//    std::cout << "Direction:" << std::endl;
//    std::cout << ray.direction() << std::endl;
//
//    std::cout << "Point:" << std::endl;
//    std::cout << point << std::endl;
//}


/*
int main() {
    BField bfield(1., 10., 10.);
    std::cout << "Bz-Field at z0 : " << bfield.getZ0() << std::endl;
    Vector3d p = {0.1, 0.1, 10};
    Vector3d b;
    b = bfield.bf(p);
    std::cout << "B-Field at z = 10 : " << b.norm();
    return 0;
}*/

