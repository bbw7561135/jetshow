//
// Created by ilya on 6/6/17.
//

#ifndef JETSHOW_JET_H
#define JETSHOW_JET_H


#include "Geometry.h"
#include "VField.h"
#include "BField.h"
#include "NField.h"

class Jet {
public:
    Jet(Geometry* geo, VField* vfield, BField* bField, NField* nField);
    // Absorption coefficient in ``point`` of the jet in the observer (lab)
    // frame. ``n`` is direction vector in the observer frame.
    double k_I(Vector3d &point, Vector3d &n);
    // Emission coefficient in ``point`` of the jet in the observer (lab)
    // frame. ``n`` is direction vector in the observer frame.
    double eta_I(Vector3d &point, Vector3d &n);

};


#endif //JETSHOW_JET_H
