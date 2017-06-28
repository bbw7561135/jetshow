#ifndef JETSHOW_CELLS_H
#define JETSHOW_CELLS_H

#include <Eigen/Eigen>
#include "Geometry.h"

using Eigen::Vector3d;
using Eigen::MatrixXd;


class Cells {
public:
		// Get ID of cell (=> direction and other possible properties) for given
		// point.
		virtual Vector3d getID(const Vector3d &point) = 0 ;
		// Restricts space for generation of cells.
		virtual void setGeometry(double coord1_max, double coord2_max,
		                         double coord1_min=0.0, double coord2_min=0.0) = 0;
		virtual void create() = 0 ;

};


class RandomCells : public Cells {
public:
		RandomCells();

		RandomCells(int N, double n_exponent, unsigned int seed=0);

		void setSeed(unsigned int seed);
		MatrixXd& getCells();

		// For random cells we just look for the nearest neighbor cell center.
		Vector3d getID(const Vector3d &point) override ;
		const int n;

protected:
		MatrixXd cells_;
		std::vector<Vector3d> directions_;
		unsigned int seed_;
		double exponent_;
};


class RandomCellsInSphere : public RandomCells {
public:
		// ``N`` - number of cells, ``n_exponent`` - power exponent for density
		// radial dependence, ``geo`` - ``Geometry`` instance to restrict possible
		// cells position.
		RandomCellsInSphere(int N, double n_exponent, unsigned int seed=0);

		void setGeometry(double coord1_max, double coord2_max,
		                 double coord1_min=0.0, double coord2_min=0.0) override ;

		// Here we implement geometry-specific generation
		void create() override ;
private:
		double r_min_;
		double r_max_;
		double theta_min_;
		double theta_max_;

};


// Thid class creates cells randomly using ``Geometry`` instance with class
// memeber ``is_within`` as restriction.
class GeneralGeometryRandomCells : public RandomCells {
public:
		GeneralGeometryRandomCells(int N, Geometry* geo, double n_exponent,
		                           unsigned int seed=0);
		// Using ``Geometry`` instance with class memeber ``is_within`` as
		// restriction.
		void create() override ;
		void setGeometry(double coord1_max, double coord2_max=0,
		                 double coord1_min=0.0, double coord2_min=0.0) override ;

private:
		Geometry* geo_;
		double r_max_;
};

#endif //JETSHOW_CELLS_H
