#include <iostream>
#include <bfields.h>

using Eigen::Vector3d;

int main() {
    BField bfield(1., 10., 10.);
    std::cout << "Bz-Field at z0 : " << bfield.getZ0() << " \n";
    Vector3d p;
    p(0) = 0.1;
    p(1) = 0.1;
    p(2) = 10.;
    Vector3d b;
    b = bfield.bf(p);
    std::cout << "B-Field at z = 10 : " << b.norm();
    return 0;
}