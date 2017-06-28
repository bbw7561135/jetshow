#include <utils.h>
#include "Cells.h"


RandomCells::RandomCells(int N, double n_exponent,
                         unsigned int seed) : seed_(seed),
                                              exponent_(n_exponent),
                                              n(N) {
	cells_.resize(3, N);
	directions_.resize(N);
}

Vector3d RandomCells::getID(const Vector3d &point) {
	MatrixXd::Index index;
	(cells_.colwise() - point).colwise().squaredNorm().minCoeff(&index);
	Vector3d result = directions_[index];
//	std::cout << "In RandomCells.getID n = " << result << std::endl;
	return result;
}

MatrixXd &RandomCells::getCells() {
	return cells_;
}

void RandomCells::setSeed(unsigned int seed) {
	seed_ = seed;
}


RandomCellsInSphere::RandomCellsInSphere(int N, double n_exponent,
                                         unsigned int seed) :
		RandomCells(N, n_exponent, seed), r_min_(), r_max_(), theta_min_(),
		theta_max_() {};

void RandomCellsInSphere::create() {
		std::vector<Vector3d> points = generate_random_points_sphere(n, r_max_,
	                                                               exponent_,
	                                                               seed_,
	                                                               r_min_,
	                                                               theta_max_);
	for (int j = 0; j < points.size()/2; ++j) {
		cells_.col(j) = points[j];
	}

	directions_ = std::move(generate_random_directions(n, seed_));
}

void RandomCellsInSphere::setGeometry(double coord1_max, double coord2_max,
                                      double coord1_min, double coord2_min) {
	r_max_ = coord1_max;
	theta_max_ = coord2_max;
	r_min_ = coord1_min;
	theta_min_ = coord2_min;
}


GeneralGeometryRandomCells::GeneralGeometryRandomCells(int N, Geometry *geo,
                                                       double n_exponent,
                                                       unsigned int seed) :
		RandomCells(N, n_exponent, seed) {
	geo_ = geo;
}

void GeneralGeometryRandomCells::create() {
	std::vector<Vector3d> points = generate_random_points_general(n, r_max_,
	                                                              geo_,
	                                                              exponent_,
	                                                              seed_);
	for (int j = 0; j < points.size(); ++j) {
		cells_.row(j) = points[j];
	}
}

void
GeneralGeometryRandomCells::setGeometry(double coord1_max, double coord2_max,
                                        double coord1_min, double coord2_min) {
r_max_ = coord1_max;
};