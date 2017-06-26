#ifndef JETSHOW_CELLS_H
#define JETSHOW_CELLS_H

#include <Eigen/Eigen>

using Eigen::Vector3d;


class Cells {
public:
		virtual Vector3d getID(const Vector3d &point) const = 0 ;
		virtual void populate() = 0 ;

};


class RandomCells : public Cells {
		Vector3d getID(const Vector3d &point) const override ;
		void populate() override ;
};

#endif //JETSHOW_CELLS_H
