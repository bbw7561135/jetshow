#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <boost/numeric/odeint.hpp>
//#include <cnpy.h>
#include "Geometry.h"
#include "Cone.h"
#include "BField.h"
#include "VField.h"
#include "NField.h"
#include "Jet.h"

#include "ImagePlane.h"
#include "Observation.h"
#include "Image.h"
#include "System.h"
#include "Ray.h"
#include "Parabaloid.h"
#include "ParabaloidCone.h"
#include "utils.h"
#include "math.h"
#include "linspace.h"
#include "Cell.h"
#include "Pixel.h"
#include <boost/range/algorithm.hpp>
#include <memory>
#include <mpi.h>
#include <ctime>
#include <chrono>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <Cells.h>


#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Interpolation_traits_2.h>
#include <CGAL/natural_neighbor_coordinates_2.h>
#include <CGAL/interpolation_functions.h>
#include <CGAL/Barycentric_coordinates_2/Triangle_coordinates_2.h>


#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

typedef CGAL::Simple_cartesian<double>                                   K_;
typedef K_::Point_2                                                Point_;
typedef CGAL::Triangulation_vertex_base_with_info_2<double, K_>      Vb;
typedef CGAL::Triangulation_data_structure_2<Vb>                  Tds;
typedef CGAL::Delaunay_triangulation_2<K_, Tds>                    Delaunay_triangulation;
typedef K_::FT                                               Coord_type;
typedef std::vector<Coord_type >                            Scalar_vector;
typedef CGAL::Barycentric_coordinates::Triangle_coordinates_2<K_> Triangle_coordinates;


using Eigen::Vector3d;
using Eigen::Matrix3Xd;
using std::vector;
using namespace boost::numeric::odeint;
namespace ph = std::placeholders;


// Pybind11
namespace py = pybind11;

//std::pair<vector<vector<double>>, vector<vector<double>>>
vector<vector<double>>
get_i_image(double los_angle, double redshift, unsigned long int number_of_pixels_along,
            unsigned long int number_of_pixels_across, double pixel_size_mas_start,
            double pixel_size_mas_stop, double cone_half_angle, double B_1, double m, double K_1,
            double n, double Gamma, double nu_observed_ghz, double tau_max, bool central_vfield) {

    // Setting geometry
    Vector3d origin = {0., 0., 0.};
    Vector3d direction = {0., 0., 1.};
    double big_scale = 100*pc;
    Cone geometry(origin, direction, cone_half_angle, 100);

    // Setting BField
    RandomScalarBFieldZ bfield(B_1, m);
    //ToroidalBField bfield(B_1, m, true);

    // Setting N_field
    BKNFieldZ nfield(K_1, n, true);

    // Setting V-field
    VField* vfield;
    if (central_vfield) {
        vfield = new ConstCentralVField(Gamma);
    } else {
        vfield = new ConstFlatVField(Gamma);
    }

    Jet bkjet(&geometry, vfield, &bfield, &nfield);

    auto image_size = std::make_pair(number_of_pixels_across, number_of_pixels_along);
    auto pc_in_mas = mas_to_pc(redshift);
    std::cout << "pc_in_mas " << pc_in_mas << std::endl;

    // Log10 of pixel size in cm
    auto lg_pixel_size_start = log10(pixel_size_mas_start*pc_in_mas*pc);
    auto lg_pixel_size_stop = log10(pixel_size_mas_stop*pc_in_mas*pc);


    ImagePlane imagePlane(image_size, lg_pixel_size_start, lg_pixel_size_stop, los_angle);
    std::cout << "Setting pixel size pc from " << pow(10.0, lg_pixel_size_start)/pc << " to " << pow(10.0, lg_pixel_size_stop)/pc << std::endl;

    // Array of pixel sizes in cm
    auto pixel_sizes = imagePlane.getPixelSizes();
    // Array of pixel solid angles in rad*rad
    std::vector<std::vector<double>> pixel_solid_angles;
    pixel_solid_angles.resize(pixel_sizes.size());

    for(unsigned long i=0; i < pixel_sizes.size(); i++) {
        pixel_solid_angles[i].resize(pixel_sizes[0].size());
        for(unsigned long j=0; j < pixel_sizes[0].size(); j++) {
            // Divide by ``pc_in_mas*pc`` to bring ``cm`` to ``mas`` at source redshift
            pixel_solid_angles[i][j] = (pixel_sizes[i][j]/(pc_in_mas*pc))*(pixel_sizes[i][j]/(pc_in_mas*pc))*mas_to_rad*mas_to_rad;
        }
    }

    // Array of scale factors. Divide resulting image on this to obtain flux density in Jy. Accounts for cosmological
    // scaling of intensity
    std::vector<std::vector<double>> scales;
    scales.resize(pixel_sizes.size());
    for(unsigned long i=0; i < pixel_sizes.size(); i++) {
        scales[i].resize(pixel_sizes[0].size());
        for(unsigned long j=0; j < pixel_sizes[0].size(); j++) {
            scales[i][j] = 1E-23*(1.+redshift)*(1.+redshift)*(1.+redshift)/pixel_solid_angles[i][j];
        }
    }


    // Setting observed frequency to the BH frame in Hz
    nu_observed_ghz *= 1E+09;
    nu_observed_ghz *= (1.+redshift);

    Observation observation(&bkjet, &imagePlane, nu_observed_ghz);

    // Transfer-specific parameters
    string step_type = "adaptive";
    double tau_n_min = 0.1;
    double dt_max_pc = 0.001;
    double dt_max = pc*dt_max_pc;
    double tau_min_log10 = -10.0;
    double tau_min = pow(10.,tau_min_log10);
    int n_ = 100;
    int n_tau_max = 2000;
    std::string calculate = "I";


    observation.run(n_, tau_max, dt_max, tau_min, step_type, calculate,
                    n_tau_max, tau_n_min, tau_max);
    string value = "tau";
    auto image_tau = observation.getImage(value);

    value = "I";
    auto image_i = observation.getImage(value);

    for (unsigned long int i = 0; i < image_i.size(); ++i)
    {
        for (unsigned long int j = 0; j < image_i[i].size(); ++j)
        {
            image_i[i][j] = image_i[i][j]/scales[i][j];
        }
    }

    //auto result = std::make_pair(image_tau, image_i);
    //return result;
    return image_i;
}


// Compile with:
// c++ -O3 -Wall -shared -std=c++14 -lCGAL -L/usr/lib/openmpi -lmpi_cxx -fopenmp -fPIC `python3 -m pybind11
// --includes` -o pyjetshow`python3-config --extension-suffix` src/pyjetshow.cpp src/BField.cpp src/Cone.cpp
// src/Geometry.cpp src/Ray.cpp src/Intersection.cpp src/Cell.cpp src/Pixel.cpp src/Image.cpp
// src/ImagePlane.cpp src/NField.cpp src/VField.cpp src/Jet.cpp src/Cylinder.cpp src/utils.cpp src/System.cpp src/Observation.cpp src/SimulationInterpolater.cpp -I/home/ilya/github/jetshow/include  -I/usr/include/eigen3
PYBIND11_MODULE(pyjetshow, m) {
using namespace pybind11::literals; // for _a literal to define arguments
m.doc() = "Radiative transfer for BK models"; // optional module docstring

m.def("get_i_image", &get_i_image, "Obtain Stokes I image with random B-field", "los_angle"_a, "redshift"_a,
"number_of_pixels_along"_a, "number_of_pixels_across"_a, "pixel_size_mas_start"_a, "pixel_size_mas_stop"_a,
"cone_half_angle"_a, "B_1"_a, "m"_a, "K_1"_a, "n"_a, "Gamma"_a, "nu_observed_ghz"_a, "tau_max"_a=10000.,
"central_vfield"_a=false);
}
