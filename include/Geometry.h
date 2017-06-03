//
// Created by ilya on 5/31/17.
//

#ifndef JETSHOW_GEOMETRY_H
#define JETSHOW_GEOMETRY_H

#include "Ray.h"
#include "Intersection.h"

class Ray;

class Geometry {
        public:
        virtual Intersection  hit(Ray &ray) = 0;

};


#endif //JETSHOW_GEOMETRY_H
