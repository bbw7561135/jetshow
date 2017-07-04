#ifndef JETSHOW_JET_H
#define JETSHOW_JET_H


#include "Geometry.h"
#include "VField.h"
#include "BField.h"
#include "NField.h"
#include "utils.h"


class Jet {
public:
    Jet(Geometry* geo, VField* vfield, BField* bField, NField* nField);

    // Vector of the magnetic field in the lab frame at point ``point``.
    const Vector3d getB(const Vector3d& point);

    // Vector of the bulk motion speed in the lab frame at point ``point``.
    const Vector3d getV(const Vector3d& point);

    // Particle density in the lab frame at point ``point``.
    const double getN(const Vector3d& point);


    // Absorption coefficient in ``point`` of the jet in the observer (lab)
    // frame. ``n`` is unit LOS vector in the observer frame.
    double getKI(Vector3d &point, Vector3d &n_los, double nu);

		double getKQ(Vector3d &point, Vector3d &n_los, double nu);

		double getKU(Vector3d &point, Vector3d &n_los, double nu);

		// TODO: For e+ must be "+" sign
		double getKV(Vector3d &point, Vector3d &n_los, double nu);

		double getKF(Vector3d &point, Vector3d &n_los, double nu);

		double getKC(Vector3d &point, Vector3d &n_los, double nu);

		double gethQ(Vector3d &point, Vector3d &n_los, double nu);

    // Emission coefficient in ``point`` of the jet in the observer (lab)
    // frame. ``n`` is unit LOS vector in the observer frame.
    double getEtaI(Vector3d &point, Vector3d &n_los, double nu);

		double getEtaQ(Vector3d &point, Vector3d &n_los, double nu);

		double getEtaU(Vector3d &point, Vector3d &n_los, double nu);

		// TODO: For e+ must be "+" sign
		double getEtaV(Vector3d &point, Vector3d &n_los, double nu);

    std::list<Intersection> hit(Ray& ray);


private:
    Geometry* geometry_;
    VField* vfield_;
    BField* bfield_;
    NField* nfield_;
};


#endif //JETSHOW_JET_H
