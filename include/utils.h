//
// Created by ilya on 6/6/17.
//

#ifndef JETSHOW_UTILS_H
#define JETSHOW_UTILS_H

#include <math.h>
#include <boost/>

const double mas_to_rad = 4.8481368*pow(10., -9.);
const double rad_to_mas = 1./mas_to_rad;
// Parsec [cm]
const double pc = 3.0857*pow(10., 18.);
// Mass of electron [g]
const double m_e = 9.109382*pow(10., -28.);
// Mass of proton [g]
const double m_p = 1.672621*pow(10., -24.);
// Charge of electron [C]
// const double q_e = 1.602176*pow(10.,-19.);
const double q_e = 4.8*pow(10.,-10.);
// Charge of proton [C]
const double q_p = 4.8*pow(10.,-10.);
// Speed of light [cm / s]
const double c = 3.*pow(10.,10.);
// Jy in cgc
const double to_jy = pow(10.,23.);



#endif //JETSHOW_UTILS_H
