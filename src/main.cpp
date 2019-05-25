#include <string>
#include <iostream>
#include <fstream>
#include <boost/numeric/odeint.hpp>
//#include <cnpy.h>
#include "ImagePlane.h"
#include "Observation.h"
#include "Geometry.h"
#include "Image.h"
#include "System.h"
#include "BField.h"
#include "VField.h"
#include "Ray.h"
#include "Cone.h"
#include "Parabaloid.h"
#include "ParabaloidCone.h"
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
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <Cells.h>
//#include "logspace.h"


// FIXME: Do not need this in analytical runs
//#include <CGAL/Simple_cartesian.h>
//#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
//#include <CGAL/Delaunay_triangulation_2.h>
//#include <CGAL/Triangulation_vertex_base_with_info_2.h>
//#include <CGAL/Interpolation_traits_2.h>
//#include <CGAL/natural_neighbor_coordinates_2.h>
//#include <CGAL/interpolation_functions.h>
//#include <CGAL/Barycentric_coordinates_2/Triangle_coordinates_2.h>
//
//typedef CGAL::Simple_cartesian<double>                                   K_;
//typedef K_::Point_2                                                Point_;
//typedef CGAL::Triangulation_vertex_base_with_info_2<double, K_>      Vb;
//typedef CGAL::Triangulation_data_structure_2<Vb>                  Tds;
//typedef CGAL::Delaunay_triangulation_2<K_, Tds>                    Delaunay_triangulation;
//typedef K_::FT                                               Coord_type;
//typedef std::vector<Coord_type >                            Scalar_vector;
//typedef CGAL::Barycentric_coordinates::Triangle_coordinates_2<K_> Triangle_coordinates;


using Eigen::Vector3d;
using Eigen::Matrix3Xd;
using std::vector;
using std::pair;
using namespace boost::numeric::odeint;
typedef std::chrono::high_resolution_clock Clock;
//namespace pt = boost::property_tree;
namespace ph = std::placeholders;





//void test_observations_rnd_bfield() {
//	// Create a root
//	pt::ptree root;
//
//	Geometry* geometry;
//	VectorBField* bfield;
//	VField* vfield;
//	NField* nfield;
//
//	// Load the json file in this ptree
//	pt::read_json("../config.json", root);
//	// Read values
//	double los_angle = root.get<double>("observation.los_angle");
//	double z = root.get<double>("observation.redshift");
//	int number_of_pixels = root.get<int>("image.number_of_pixels");
//	double pixel_size_mas = root.get<double>("image.pixel_size_mas");
//
//	// Setting geometry
//	Vector3d origin = {0., 0., 0.};
//	Vector3d direction = {0., 0., 1.};
//	std::string geotype = root.get<std::string>("jet.geometry.type");
//	std::cout << "Geometry type " << geotype << std::endl;
//	if (geotype == "cone") {
//		// Create Cone with parameters
//		double cone_angle = root.get<double>("jet.geometry.parameters.angle");
//		double scale_pc = root.get<double>("jet.geometry.parameters.scale_pc");
//		double scale = pc * scale_pc;
//		geometry = new Cone(origin, direction, cone_angle, scale);
////				geometry = geometry_;
//	}
//
//	std::string btype = root.get<std::string>("jet.bfield.type");
//	std::cout << "B-field type : " << btype << std::endl;
//	if (btype == "radial_conical") {
//		double b_1 = root.get<double>("jet.bfield.parameters.b_1");
//		double n_b = root.get<double>("jet.bfield.parameters.n_b");
//		bfield = new RadialConicalBField(b_1, n_b);
//	} else if (btype == "spiral_conical") {
//		double b_1 = root.get<double>("jet.bfield.parameters.b_1");
//		double pitch_angle = root.get<double>("jet.bfield.parameters.pitch_angle");
//		bfield = new SpiralConicalBField(b_1, pitch_angle);
//	};
//
//	// Create ``Cells`` instance
////	int N = 1000;
////	RandomCellsInSphere cells(N, 2.0);
////	double cone_angle = root.get<double>("jet.geometry.parameters.angle");
////	double scale_pc = root.get<double>("jet.geometry.parameters.scale_pc");
////	cells.setGeometry(scale_pc*pc, cone_angle);
////	cells.create();
//	// Create ``RandomVectorBField`` instance
//	double rnd_fraction = root.get<double>("jet.bfield.parameters.random_fraction");
////	RandomCellsBField rnd_bfield(&cells, bfield, rnd_fraction);
//	unsigned int seed = 123;
//	if (rnd_fraction > 0) {
//		bfield = new PointsRandomVectorBField(bfield, rnd_fraction, seed);
//	}
//
//	std::string vtype = root.get<std::string>("jet.vfield.type");
//	std::cout << "Velocity type : " << vtype << std::endl;
//	if (vtype == "const_central") {
//		double gamma = root.get<double>("jet.vfield.parameters.gamma");
//		vfield = new ConstCentralVField(gamma);
//	}
//
//	std::string ntype = root.get<std::string>("jet.nfield.type");
//	std::cout << "Density type : " << ntype << std::endl;
//	if (ntype == "bk") {
//		double n_1 = root.get<double>("jet.nfield.parameters.n_1");
//		double n_n = root.get<double>("jet.nfield.parameters.n_n");
//		nfield = new BKNField(n_1, n_n);
//	}
//
//
//	Jet bkjet(geometry, vfield, bfield, nfield);
//
//	auto image_size = std::make_pair(number_of_pixels, number_of_pixels);
//	auto pc_in_mas = mas_to_pc(z);
////		auto cm_in_mas = pc * pc_in_mas;
//	auto pixel_size = pixel_size_mas*0.25*pc_in_mas*pc;
//	auto pix_solid_angle = pixel_solid_angle(pixel_size_mas, z);
//
//	ImagePlane imagePlane(image_size, pixel_size, pixel_size, los_angle);
//
//	double nu = root.get<double>("observation.frequency_ghz");
//	nu *= 1E+09;
//	// Redshifting to SBH frame
//	nu *= (1.+z);
//	Observation observation(&bkjet, &imagePlane, nu);
//	string step_type = root.get<string>("integration.step_type");
//	double tau_max = root.get<double>("integration.parameters.tau_max");
//	double tau_n_min = root.get<double>("integration.parameters.tau_n_min");
//	double dt_max_pc = root.get<double>("integration.parameters.dl_max_pc");
//	double dt_max = pc*dt_max_pc;
//	double tau_min_log10 = root.get<double>("integration.parameters.log10_tau_min");
//	double tau_min = pow(10.,tau_min_log10);
//	int n = root.get<int>("integration.parameters.n");
//	int n_tau_max = root.get<int>("integration.parameters.n_tau_max");
//
//	std::cout << "Integrating using max. opt. depth = " << tau_max << std::endl;
//	std::cout << "Integrating using min. lg(opt.depth) = " << tau_min_log10 << std::endl;
//	std::cout << "Integrating using max. step [pc] = " << dt_max_pc << std::endl;
//	std::cout << "Integrating using default number of steps = " << n << std::endl;
//
//
//	observation.run(n, tau_max, dt_max, tau_min, step_type, n_tau_max, tau_n_min,
//	                tau_max);
//
//	string value = "tau";
//	auto image = observation.getImage(value);
//	std::fstream fs;
//	std::string file_tau = root.get<std::string>("output.file_tau");
//	fs.open(file_tau, std::ios::out | std::ios::app);
//
//	if (fs.is_open())
//	{
//		write_2dvector(fs, image);
//		fs.close();
//	}
//
//	value = "I";
//	image = observation.getImage(value);
//	std::string file_i = root.get<std::string>("output.file_i");
//	fs.open(file_i, std::ios::out | std::ios::app);
//	double scale = 1E-23/pix_solid_angle;
//	std::cout << "Scaling Stokes I by " << scale << std::endl;
//
//	if (fs.is_open())
//	{
//		// Scaling to Jy
//		write_2dvector(fs, image, scale);
//		// Just to show how it can be used
//		// write_2dvector(std::cout, image);
//		fs.close();
//	}
//
//	value = "l";
//	image = observation.getImage(value);
//	std::string file_length = root.get<std::string>("output.file_length");
//	fs.open(file_length, std::ios::out | std::ios::app);
//
//	if (fs.is_open())
//	{
//		write_2dvector(fs, image, pc);
//		fs.close();
//	}
//}

//
//void test_observations_full() {
//	// Create a root
//	pt::ptree root;
//
//	Geometry* geometry;
//	VectorBField* bfield;
//	VField* vfield;
//	NField* nfield;
//
//	// Load the json file in this ptree
//	pt::read_json("../config.json", root);
//	// Read values
//	double los_angle = root.get<double>("observation.los_angle");
//	double z = root.get<double>("observation.redshift");
//	int number_of_pixels = root.get<int>("image.number_of_pixels");
//	double pixel_size_mas = root.get<double>("image.pixel_size_mas");
//
//	// Setting geometry
//	Vector3d origin = {0., 0., 0.};
//	Vector3d direction = {0., 0., 1.};
//	std::string geotype = root.get<std::string>("jet.geometry.type");
//	std::cout << "Geometry type " << geotype << std::endl;
//	if (geotype == "cone") {
//		// Create Cone with parameters
//		double cone_angle = root.get<double>("jet.geometry.parameters.angle");
//		double scale_pc = root.get<double>("jet.geometry.parameters.scale_pc");
//		double scale = pc * scale_pc;
//		geometry = new Cone(origin, direction, cone_angle, scale);
////				geometry = geometry_;
//	}
//
//	std::string btype = root.get<std::string>("jet.bfield.type");
//	std::cout << "B-field type : " << btype << std::endl;
//	if (btype == "radial_conical") {
//		double b_1 = root.get<double>("jet.bfield.parameters.b_1");
//		double n_b = root.get<double>("jet.bfield.parameters.n_b");
//		bfield = new RadialConicalBField(b_1, n_b);
//	} else if (btype == "spiral_conical") {
//		double b_1 = root.get<double>("jet.bfield.parameters.b_1");
//		double pitch_angle = root.get<double>("jet.bfield.parameters.pitch_angle");
//		bfield = new SpiralConicalBField(b_1, pitch_angle);
//	} else if (btype == "toroidal") {
//		double b_1 = root.get<double>("jet.bfield.parameters.b_1");
//		double n_b = root.get<double>("jet.bfield.parameters.n_b");
//		bfield = new ToroidalBField(b_1, n_b);
//	};
//
//	// Create ``Cells`` instance
////	int N = 1000;
////	RandomCellsInSphere cells(N, 2.0);
////	double cone_angle = root.get<double>("jet.geometry.parameters.angle");
////	double scale_pc = root.get<double>("jet.geometry.parameters.scale_pc");
////	cells.setGeometry(scale_pc*pc, cone_angle);
////	cells.create();
//	// Create ``RandomVectorBField`` instance
//	double rnd_fraction = root.get<double>("jet.bfield.parameters.random_fraction");
////	RandomCellsBField rnd_bfield(&cells, bfield, rnd_fraction);
//	unsigned int seed = 123;
//	if (rnd_fraction > 0) {
//		bfield = new PointsRandomVectorBField(bfield, rnd_fraction, seed);
//	}
//
//	std::string vtype = root.get<std::string>("jet.vfield.type");
//	std::cout << "Velocity type : " << vtype << std::endl;
//	if (vtype == "const_central") {
//		double gamma = root.get<double>("jet.vfield.parameters.gamma");
//		vfield = new ConstCentralVField(gamma);
//	}
//
//	std::string ntype = root.get<std::string>("jet.nfield.type");
//	std::cout << "Density type : " << ntype << std::endl;
//	if (ntype == "bk") {
//		double n_1 = root.get<double>("jet.nfield.parameters.n_1");
//		double n_n = root.get<double>("jet.nfield.parameters.n_n");
//		nfield = new BKNField(n_1, n_n);
//	}
//
//
//	Jet bkjet(geometry, vfield, bfield, nfield);
//
//	auto image_size = std::make_pair(number_of_pixels, number_of_pixels);
//	auto pc_in_mas = mas_to_pc(z);
////		auto cm_in_mas = pc * pc_in_mas;
//	auto pixel_size = pixel_size_mas*0.25*pc_in_mas*pc;
//	auto pix_solid_angle = pixel_solid_angle(pixel_size_mas, z);
//
//	ImagePlane imagePlane(image_size, pixel_size, pixel_size, los_angle);
//
//	double nu = root.get<double>("observation.frequency_ghz");
//	nu *= 1E+09;
//	// Redshifting to SBH frame
//	nu *= (1.+z);
//	Observation observation(&bkjet, &imagePlane, nu);
//	string step_type = root.get<string>("integration.step_type");
//	double tau_max = root.get<double>("integration.parameters.tau_max");
//	double tau_n_min = root.get<double>("integration.parameters.tau_n_min");
//	double dt_max_pc = root.get<double>("integration.parameters.dl_max_pc");
//	double dt_max = pc*dt_max_pc;
//	double tau_min_log10 = root.get<double>("integration.parameters.log10_tau_min");
//	double tau_min = pow(10.,tau_min_log10);
//	int n = root.get<int>("integration.parameters.n");
//	int n_tau_max = root.get<int>("integration.parameters.n_tau_max");
//
//	std::cout << "Integrating using max. opt. depth = " << tau_max << std::endl;
//	std::cout << "Integrating using min. lg(opt.depth) = " << tau_min_log10 << std::endl;
//	std::cout << "Integrating using max. step [pc] = " << dt_max_pc << std::endl;
//	std::cout << "Integrating using default number of steps = " << n << std::endl;
//
//
//	observation.run(n, tau_max, dt_max, tau_min, step_type, n_tau_max, tau_n_min,
//	                tau_max);
//
//	string value = "tau";
//	auto image = observation.getImage(value);
//	std::fstream fs;
//	std::string file_tau = root.get<std::string>("output.file_tau");
//	fs.open(file_tau, std::ios::out | std::ios::app);
//
//	if (fs.is_open())
//	{
//		write_2dvector(fs, image);
//		fs.close();
//	}
//
//	value = "I";
//	image = observation.getImage(value);
//	std::string file_i = root.get<std::string>("output.file_i");
//	fs.open(file_i, std::ios::out | std::ios::app);
//	double scale = 1E-23/pix_solid_angle;
//	std::cout << "Scaling Stokes I by " << scale << std::endl;
//
//	if (fs.is_open())
//	{
//		// Scaling to Jy
//		write_2dvector(fs, image, scale);
//		// Just to show how it can be used
//		// write_2dvector(std::cout, image);
//		fs.close();
//	}
//
//	value = "Q";
//	image = observation.getImage(value);
//	std::string file_q = root.get<std::string>("output.file_q");
//	fs.open(file_q, std::ios::out | std::ios::app);
//
//	if (fs.is_open())
//	{
//		// Scaling to Jy
//		write_2dvector(fs, image, scale);
//		// Just to show how it can be used
//		// write_2dvector(std::cout, image);
//		fs.close();
//	}
//
//	value = "U";
//	image = observation.getImage(value);
//	std::string file_u = root.get<std::string>("output.file_u");
//	fs.open(file_u, std::ios::out | std::ios::app);
//
//	if (fs.is_open())
//	{
//		// Scaling to Jy
//		write_2dvector(fs, image, scale);
//		// Just to show how it can be used
//		// write_2dvector(std::cout, image);
//		fs.close();
//	}
//
//	value = "V";
//	image = observation.getImage(value);
//	std::string file_v = root.get<std::string>("output.file_v");
//	fs.open(file_v, std::ios::out | std::ios::app);
//
//	if (fs.is_open())
//	{
//		// Scaling to Jy
//		write_2dvector(fs, image, scale);
//		// Just to show how it can be used
//		// write_2dvector(std::cout, image);
//		fs.close();
//	}
//
//	value = "l";
//	image = observation.getImage(value);
//	std::string file_length = root.get<std::string>("output.file_length");
//	fs.open(file_length, std::ios::out | std::ios::app);
//
//	if (fs.is_open())
//	{
//		write_2dvector(fs, image, pc);
//		fs.close();
//	}
//}



// THIS IS MASTER FUNCTION
//void test_stripe() {
//	// Create a root
//	pt::ptree root;
//
//	Geometry* geometry;
////	RandomScalarBField* bfield;
//	VectorBField* bfield;
//	VField* vfield;
//	NField* nfield;
//
//	// Load the json file in this ptree
//	pt::read_json("../config.json", root);
//	// Read values
//	double los_angle = root.get<double>("observation.los_angle");
//	double z = root.get<double>("observation.redshift");
//	int number_of_pixels = root.get<int>("image.number_of_pixels");
//	double pixel_size_mas = root.get<double>("image.pixel_size_mas");
//
//	// Setting geometry
//	Vector3d origin = {0., 0., 0.};
//	Vector3d direction = {0., 0., 1.};
//	std::string geotype = root.get<std::string>("jet.geometry.type");
//	std::cout << "Geometry type " << geotype << std::endl;
//	if (geotype == "cone") {
//		// Create Cone with parameters
//		double cone_angle = root.get<double>("jet.geometry.parameters.angle");
//		double scale_pc = root.get<double>("jet.geometry.parameters.scale_pc");
//		double scale = pc * scale_pc;
//		geometry = new Cone(origin, direction, cone_angle, scale);
////				geometry = geometry_;
//	}
//
//
//	std::string btype = root.get<std::string>("jet.bfield.type");
//	std::cout << "B-field type : " << btype << std::endl;
////	if (btype == "radial_conical") {
////		double b_1 = root.get<double>("jet.bfield.parameters.b_1");
////		double n_b = root.get<double>("jet.bfield.parameters.n_b");
////		bfield = new RadialConicalBField(b_1, n_b);
////	} else if (btype == "spiral_conical") {
////		double b_1 = root.get<double>("jet.bfield.parameters.b_1");
////		double pitch_angle = root.get<double>("jet.bfield.parameters.pitch_angle");
////		bfield = new SpiralConicalBField(b_1, pitch_angle);
////	} else if (btype == "toroidal") {
////		double b_1 = root.get<double>("jet.bfield.parameters.b_1");
////		double n_b = root.get<double>("jet.bfield.parameters.n_b");
////		bfield = new ToroidalBField(b_1, n_b);
////	};
//
//	double b_1 = root.get<double>("jet.bfield.parameters.b_1");
//	double n_b = root.get<double>("jet.bfield.parameters.n_b");
//	std::cout << "B-field value, exponent : " << b_1 << " " << n_b << std::endl;
//	bfield = new RadialConicalBField(b_1, n_b);
//
////	double rnd_fraction = root.get<double>("jet.bfield.parameters.random_fraction");
////	if (rnd_fraction > 0.01) {
////		// Create ``Cells`` instance
////		int N = root.get<int>("jet.bfield.parameters.n_cells_1pc");
////		RandomCellsInSphere cells(N, 2.0);
////		double cone_angle = root.get<double>("jet.geometry.parameters.angle");
////		double scale_pc = root.get<double>("jet.geometry.parameters.scale_pc");
////		cells.setGeometry(scale_pc * pc, cone_angle);
////		cells.create();
////		// Create ``RandomBField`` instance
////		RandomCellsBField rnd_bfield(&cells, bfield, rnd_fraction);
////		unsigned int seed = 123;
////		if (rnd_fraction > 0) {
////			bfield = new RandomPointBField(bfield, rnd_fraction, seed);
////		}
////	}
//
//	std::string vtype = root.get<std::string>("jet.vfield.type");
//	std::cout << "Velocity type : " << vtype << std::endl;
//	if (vtype == "const_central") {
//		double gamma = root.get<double>("jet.vfield.parameters.gamma");
//		vfield = new ConstCentralVField(gamma);
//	} else if (vtype == "const_flat") {
//		double gamma = root.get<double>("jet.vfield.parameters.gamma");
//		vfield = new ConstFlatVField(gamma);
//	};
//
//	std::string ntype = root.get<std::string>("jet.nfield.type");
//	std::cout << "Density type : " << ntype << std::endl;
//	if (ntype == "bk") {
//		double n_1 = root.get<double>("jet.nfield.parameters.n_1");
//		double n_n = root.get<double>("jet.nfield.parameters.n_n");
//		nfield = new BKNField(n_1, n_n);
//	}
//
//
//	Jet bkjet(geometry, vfield, bfield, nfield);
//
//	auto image_size = std::make_pair(number_of_pixels, number_of_pixels);
//	auto pc_in_mas = mas_to_pc(z);
////		auto cm_in_mas = pc * pc_in_mas;
//	std::cout << "Setting pixel size to " << pixel_size_mas << " mas" << std::endl;
//	auto pixel_size = pixel_size_mas*pc_in_mas*pc;
//	auto pix_solid_angle = pixel_solid_angle(pixel_size_mas, z);
//
//	ImagePlane imagePlane(image_size, pixel_size, pixel_size, los_angle);
//	std::cout << "Setting pixel size pc " << pixel_size/pc << std::endl;
//
//	double nu = root.get<double>("observation.frequency_ghz");
//	nu *= 1E+09;
//	// Redshifting to SBH frame
//	nu *= (1.+z);
//	Observation observation(&bkjet, &imagePlane, nu);
//	string step_type = root.get<string>("integration.step_type");
//	double tau_max = root.get<double>("integration.parameters.tau_max");
//	double tau_n_min = root.get<double>("integration.parameters.tau_n_min");
//	double dt_max_pc = root.get<double>("integration.parameters.dl_max_pc");
//	double dt_max = pc*dt_max_pc;
//	double tau_min_log10 = root.get<double>("integration.parameters.log10_tau_min");
//	double tau_min = pow(10.,tau_min_log10);
//	int n = root.get<int>("integration.parameters.n");
//	int n_tau_max = root.get<int>("integration.parameters.n_tau_max");
//	std::string calculate = root.get<std::string>("calculate");
//
//	std::cout << "Integrating using max. opt. depth = " << tau_max << std::endl;
//	std::cout << "Integrating using min. lg(opt.depth) = " << tau_min_log10 << std::endl;
//	std::cout << "Integrating using max. step [pc] = " << dt_max_pc << std::endl;
//	std::cout << "Integrating using default number of steps = " << n << std::endl;
//	std::cout << "Calculate mode = " << calculate << std::endl;
//
//
//
//	if (calculate == "tau") {
//		observation.run_stripe(n, tau_max, tau_min);
//
//		string value = "tau";
//		auto stripe = observation.getStripe(value);
//		std::fstream fs;
//		std::string file_tau = root.get<std::string>("output.file_tau_stripe");
//		fs.open(file_tau, std::ios::out | std::ios::app);
//
//		if (fs.is_open()) {
//			write_vector(fs, stripe);
//			fs.close();
//		}
//	}
//
//
//	else {
//		observation.run(n, tau_max, dt_max, tau_min, step_type, calculate,
//		                n_tau_max, tau_n_min, tau_max);
//		string value = "tau";
//		auto image = observation.getImage(value);
//		std::fstream fs;
//		std::string file_tau = root.get<std::string>("output.file_tau");
//		fs.open(file_tau, std::ios::out | std::ios::app);
//
//		if (fs.is_open())
//		{
//			write_2dvector(fs, image);
//			fs.close();
//		}
//
//		value = "I";
//		image = observation.getImage(value);
//		std::string file_i = root.get<std::string>("output.file_i");
//		fs.open(file_i, std::ios::out | std::ios::app);
//    // double scale = 1E-23/pix_solid_angle;
//		// I_{\nu}/\nu^3 = inv => to get Flux in the observer frame in Jy:
//		double scale = 1E-23*(1.+z)*(1.+z)*(1.+z)/pix_solid_angle;
//
//		// Now convert to observed flux using luminosity distance
////		double da = comoving_transfer_distance(z)/(1.+z);
////		double da = da_old(z);
////		double scale = 1E-23*da*da*(1.+z)*(1.+z)*(1.+z)/(pixel_size_mas*pc_in_mas*pixel_size_mas*pc_in_mas);
//
//		std::cout << "Scaling Stokes I by " << scale << std::endl;
//
//		if (fs.is_open())
//		{
//			// Scaling to Jy
//			write_2dvector(fs, image, scale);
//			// Just to show how it can be used
//			// write_2dvector(std::cout, image);
//			fs.close();
//		}
//
//		// TODO: Don't save l
////		value = "l";
////		image = observation.getImage(value);
////		std::string file_length = root.get<std::string>("output.file_length");
////		fs.open(file_length, std::ios::out | std::ios::app);
////
////		if (fs.is_open()) {
////			write_2dvector(fs, image, pc);
////			fs.close();
////		}
//
//
//		if (calculate == "full") {
//			value = "Q";
//			image = observation.getImage(value);
//			std::string file_q = root.get<std::string>("output.file_q");
//			fs.open(file_q, std::ios::out | std::ios::app);
//
//			if (fs.is_open()) {
//				// Scaling to Jy
//				write_2dvector(fs, image, scale);
//				// Just to show how it can be used
//				// write_2dvector(std::cout, image);
//				fs.close();
//			}
//
//			value = "U";
//			image = observation.getImage(value);
//			std::string file_u = root.get<std::string>("output.file_u");
//			fs.open(file_u, std::ios::out | std::ios::app);
//
//			if (fs.is_open()) {
//				// Scaling to Jy
//				write_2dvector(fs, image, scale);
//				// Just to show how it can be used
//				// write_2dvector(std::cout, image);
//				fs.close();
//			}
//
//			value = "V";
//			image = observation.getImage(value);
//			std::string file_v = root.get<std::string>("output.file_v");
//			fs.open(file_v, std::ios::out | std::ios::app);
//
//			if (fs.is_open()) {
//				// Scaling to Jy
//				write_2dvector(fs, image, scale);
//				// Just to show how it can be used
//				// write_2dvector(std::cout, image);
//				fs.close();
//			}
//
//			// TODO: Don't save l
////			value = "l";
////			image = observation.getImage(value);
////			std::string file_length = root.get<std::string>("output.file_length");
////			fs.open(file_length, std::ios::out | std::ios::app);
////
////			if (fs.is_open()) {
////				write_2dvector(fs, image, pc);
////				fs.close();
////			}
//		}
//	}
//}
//


//void test_velocity() {
//    VField* vfield;
//    // Radius of parabaloid at z0=1pc
//    double r0 = 0.1*pc;
//    double gamma_axis0 = 30;
//    double gamma_border0 = 1;
//    // Z at which the speed on the axis is gamma0_axis
//    double z_at_gamma0 = 1*pc;
//    // Radius of parabaloid at z = 1pc
//    double Rz0 = r0;
//    vfield = new ShearedAccParabolicVField(gamma_axis0, gamma_border0, z_at_gamma0, Rz0);
//    Vector3d point = {0, 0.025*pc, 1*pc};
//    std::cout << vfield->vf(point) << std::endl;
//
//}


// Testing geometries, composite fields
//void test_collimations() {
//    Geometry* geometry;
////    BField* bfield;
//    RandomScalarBFieldZ* bfield;
//    VField* vfield;
//    NField* nfield;
//
//    double los_angle = pi/2;
//    double z = 0.00436;
//    int number_of_pixels = 1024;
//    double pixel_size_mas = 0.005;
//
//    // Setting geometry
//    Vector3d origin = {0., 0., 0.};
//    Vector3d direction = {0., 0., 1.};
////    double cone_half_angle = 0.1;
//    double big_scale = 100*pc;
//    // Radius of parabaloid at z0=1pc
//    double r0 = 0.05*pc;
//    // Distance where collimation stops (100 pix for 0.005mas pixel)
//    double z0 = 0.05*pc;
////    geometry = new Cone(origin, direction, cone_half_angle, big_scale);
//    geometry = new Parabaloid(origin, direction, r0, big_scale);
////    geometry = new ParabaloidCone(origin, direction, r0, z0, big_scale);
//
//    double b_1 = 0.01;
//    double n_b = 1.35;
////    bfield = new RadialConicalBField(b_1, n_b);
//    bfield = new RandomScalarBFieldZ(b_1, n_b);
//
////    double gamma = 10;
////    double z_at_gamma = 1*pc;
////    vfield = new ConstCentralVField(gamma);
////    vfield = new AccParabolicVField(gamma, z_at_gamma);
//
////    double gamma_axis0 = 10;
////    double gamma_border0 = 1;
////    // Z at which the speed on the axis is gamma_0_axis
////    double z_at_gamma0 = 0.1*pc;
////    // Radius of parabaloid at z = 1pc
////    double Rz0 = r0;
////
////    vfield = new ShearedAccParabolicVField(gamma_axis0, gamma_border0, z_at_gamma0, Rz0);
//
////    double gamma0 = 10;
////    vfield = new AccParabolicConstConeVField(gamma0, r0, z0);
//
//
//    double gamma_axis0 = 10;
//    double gamma_border0 = 1;
//    // Z at which the speed on the axis is gamma_0_axis
//    double z_at_gamma0 = 0.05*pc;
//    // Radius of parabaloid at z = 1pc
//    double Rz0 = r0;
//
//    vfield = new ShearedAccParabolicConstConeVField(gamma_axis0, gamma_border0, Rz0, z_at_gamma0);
//
//    double n_1 = 10;
//    double n_n = 2.0;
//    nfield = new BKNField(n_1, n_n);
//
//    Jet bkjet(geometry, vfield, bfield, nfield);
//
//    auto image_size = std::make_pair(number_of_pixels, number_of_pixels);
//    auto pc_in_mas = mas_to_pc(z);
//    auto pixel_size = pixel_size_mas*pc_in_mas*pc;
//    auto pix_solid_angle = pixel_solid_angle(pixel_size_mas, z);
//
//    ImagePlane imagePlane(image_size, pixel_size, pixel_size, los_angle);
//    std::cout << "Setting pixel size pc " << pixel_size/pc << std::endl;
//
//    // Setting observed frequency
//    double nu = 15.4;
//    nu *= 1E+09;
//    nu *= (1.+z);
//
//    Observation observation(&bkjet, &imagePlane, nu);
//
//    string step_type = "adaptive";
//    double tau_max = 100.0;
//    double tau_n_min = 0.1;
//    double dt_max_pc = 0.001;
//    double dt_max = pc*dt_max_pc;
//    double tau_min_log10 = -8.0;
//    double tau_min = pow(10.,tau_min_log10);
//    int n = 100;
//    int n_tau_max = 2000;
//    std::string calculate = "I";
//
//
//    observation.run(n, tau_max, dt_max, tau_min, step_type, calculate,
//                    n_tau_max, tau_n_min, tau_max);
//    string value = "tau";
//    auto image = observation.getImage(value);
//    std::fstream fs;
//    std::string file_tau = "map_tau.txt";
//    fs.open(file_tau, std::ios::out | std::ios::app);
//
//    if (fs.is_open())
//    {
//        write_2dvector(fs, image);
//        fs.close();
//    }
//
//    value = "I";
//    image = observation.getImage(value);
//    std::string file_i = "map_i.txt";
//    fs.open(file_i, std::ios::out | std::ios::app);
//    double scale = 1E-23*(1.+z)*(1.+z)*(1.+z)/pix_solid_angle;
//
//    if (fs.is_open())
//    {
//        // Scaling to Jy
//        write_2dvector(fs, image, scale);
//        fs.close();
//    }
//
//    value = "l";
//    image = observation.getImage(value);
//    std::string file_length = "map_l.txt";
//    fs.open(file_length, std::ios::out | std::ios::app);
//
//    if (fs.is_open()) {
//        write_2dvector(fs, image, pc);
//        fs.close();
//    }
//}


//void test_reading_npy() {
//    cnpy::NpyArray arr = cnpy::npy_load("gamma_full.npy");
//    double* loaded_data = arr.data<double>();
//    size_t nrows = arr.shape[0];
//    size_t ncols = arr.shape[1];
//    std::cout << "nrows=" << nrows << ", ncols=" << ncols << std::endl;
//    size_t row = 310000;
//    double z = loaded_data[row*nrows + 0]/pc;
//    double r_p = loaded_data[row*nrows + 1]/pc;
//    std::cout << z << " " << r_p << std::endl;
//}


//void test_interpolation() {
//    Delaunay_triangulation tr_v;
//    create_triangulation("vfield_10.txt", &tr_v);
//    SimulationVField vfield(&tr_v);
//
//    size_t n_points = 100;
//    // Grid of r and r_p where to interpolate
//    std::vector<double> r_grid = MyLinearSpacedArray(0., 4*pow(10, 21), 10*n_points);
//    std::vector<double> r_p_grid = MyLinearSpacedArray(0., pow(10, 19), n_points);
//    std::vector<vector<double >> result;
//    for (auto & r: r_grid) {
//        for (auto & r_p : r_p_grid) {
//            Vector3d point1(0, r_p, r);
//            Vector3d point2(0, -r_p, r);
//            Vector3d v1 = vfield.vf(point1);
//            Vector3d v2 = vfield.vf(point2);
//            double gamma1 = sqrt(1./(1.- (v1/c).squaredNorm()));
//            double gamma2 = sqrt(1./(1.- (v2/c).squaredNorm()));
//            result.emplace_back(std::vector<double>{r, r_p, gamma1});
//            result.emplace_back(std::vector<double>{r, -r_p, gamma2});
//        }
//    }
//
//    std::fstream fs;
//    std::string file_length = "interpolated.txt";
//    fs.open(file_length, std::ios::out | std::ios::app);
//
//    if (fs.is_open()) {
//        write_2dvector(fs, result);
//        fs.close();
//    }
//}


//void run_on_simulations() {
//    // M87
////    double los_angle = 0.314;
//    // Some blazar -- 2.5 deg
//    double los_angle = 0.0436;
//    // Some blazar -- 5 deg
////    double los_angle = 0.0872;
//    double redshift = 0.00436;
////    unsigned long int number_of_pixels_along = 500;
////    unsigned long int number_of_pixels_across = 100;
//    unsigned long int number_of_pixels_along = 200;
//    unsigned long int number_of_pixels_across = 200;
//    double pixel_size_mas = 0.5;
//
//
//    // Setting geometry
//    // Using simulations
//    vector< vector<double> > all_points;
////    read_from_txt("vfield_10.txt", all_points);
//    read_from_txt("vfield_10_v2.txt", all_points);
//    size_t nrows = all_points.size();
//
//    std::vector<Point_3> points;
//    int n_circle = 36;
//    std::cout << "Reading geometry file with #rows = " << nrows << std::endl;
//    for (size_t i=0; i<nrows; i++) {
//        double z = all_points[i][0]/pc;
//        double r_p = all_points[i][1]/pc;
//        for (int j=0; j<n_circle; j++) {
//            double x = r_p*sin(j*2*pi/n_circle);
//            double y = r_p*cos(j*2*pi/n_circle);
//            double length_ = sqrt(x*x + y*y + z*z);
//            points.emplace_back(Point_3(x, y, z));
//        }
//    }
//
//    Polyhedron P;
//    CGAL::convex_hull_3(points.begin(), points.end(), P);
//    Tree tree(faces(P).first, faces(P).second, P);
//    SimulationGeometry geometry(&tree);
//
////    // Using analytical shapes
////    Vector3d origin = {0., 0., 0.};
////	Vector3d direction = {0., 0., 1.};
////    double cone_half_angle = 0.01;
////    double big_scale = 100*pc;
////    // Radius of parabaloid at z0=1pc
////    double r0 = 0.1*pc;
////    // Distance where collimation stops
////    double z0 = 3.0*pc;
//////    Cone geometry(origin, direction, cone_half_angle, 100);
////    ParabaloidCone geometry(origin, direction, r0, z0, big_scale);
//
//
//    // Setting VectorBField
//    // Using simulations
//    Delaunay_triangulation tr_p;
//    Delaunay_triangulation tr_fi;
////    create_triangulation("bfield_p_10.txt", &tr_p);
////    create_triangulation("bfield_fi_10.txt", &tr_fi);
//    create_triangulation("bfield_p_10_v2.txt", &tr_p);
//    create_triangulation("bfield_fi_10_v2.txt", &tr_fi);
//    SimulationBField bfield(&tr_p, &tr_fi, false);
////    // Using analytical expressions
//////    RadialConicalBField bfield(0.1, 1, true);
//////    RandomScalarBFieldZ bfield(1, 1.35);
////    CompositeRandomScalarBFieldZ bfield(1, 0.5, 1, z0);
//
//
//    // Setting N_field
//    // Using simulations
//    Delaunay_triangulation tr_n;
////    create_triangulation("nfield_10.txt", &tr_n);
//    create_triangulation("nfield_10_v2.txt", &tr_n);
//    SimulationNField nfield(&tr_n, true);
////    // Using analytical expressions
////    BKNField nfield(1000, 2, true);
////    CompositeBKNField nfield(1000, 1, 2, z0, true);
//
//
//    // Setting V-field
//    // Using simulations
//    Delaunay_triangulation tr_v;
////    create_triangulation("vfield_10.txt", &tr_v);
//    create_triangulation("vfield_10_v2.txt", &tr_v);
//    SimulationVField vfield(&tr_v);
////    // Using analytical expressions
//////    ConstCentralVField vfield(10);
////    AccParabolicConstConeVField vfield(10, r0, z0);
//
//    Jet bkjet(&geometry, &vfield, &bfield, &nfield);
//
//    auto image_size = std::make_pair(number_of_pixels_across, number_of_pixels_along);
//    auto pc_in_mas = mas_to_pc(redshift);
//    auto pixel_size = pixel_size_mas*pc_in_mas*pc;
//    auto pix_solid_angle = pixel_solid_angle(pixel_size_mas, redshift);
//
//    ImagePlane imagePlane(image_size, pixel_size, pixel_size, los_angle);
//    std::cout << "Setting pixel size pc " << pixel_size/pc << std::endl;
//
//    // Setting observed frequency
//    double nu = 15.4;
//    nu *= 1E+09;
//    nu *= (1.+redshift);
//
//    Observation observation(&bkjet, &imagePlane, nu);
//
//    string step_type = "adaptive";
//    double tau_max = 100.0;
//    double tau_n_min = 0.1;
//    double dt_max_pc = 0.001;
//    double dt_max = pc*dt_max_pc;
//    double tau_min_log10 = -10.0;
//    double tau_min = pow(10.,tau_min_log10);
//    int n = 100;
//    int n_tau_max = 2000;
//    std::string calculate = "I";
//
//
//    observation.run(n, tau_max, dt_max, tau_min, step_type, calculate,
//                    n_tau_max, tau_n_min, tau_max);
//    string value = "tau";
//    auto image = observation.getImage(value);
//    std::fstream fs;
//    std::string file_tau = "map_tau.txt";
//    fs.open(file_tau, std::ios::out | std::ios::app);
//
//    if (fs.is_open())
//    {
//        write_2dvector(fs, image);
//        fs.close();
//    }
//
//    value = "I";
//    image = observation.getImage(value);
//    std::string file_i = "map_i.txt";
//    fs.open(file_i, std::ios::out | std::ios::app);
//    double scale = 1E-23*(1.+redshift)*(1.+redshift)*(1.+redshift)/pix_solid_angle;
//
//    if (fs.is_open())
//    {
//        // Scaling to Jy
//        write_2dvector(fs, image, scale);
//        fs.close();
//    }
//
//    value = "l";
//    image = observation.getImage(value);
//    std::string file_length = "map_l.txt";
//    fs.open(file_length, std::ios::out | std::ios::app);
//
//    if (fs.is_open()) {
//        write_2dvector(fs, image, pc);
//        fs.close();
//    }
//}


//void run_analytical() {
//    // 38 deg
//    //double los_angle = 0.66323;
//    double los_angle = pi/2;
//    //double redshift = 0.0165;
//    double redshift = 0.5;
//    unsigned long int number_of_pixels_along = 5000;
//    unsigned long int number_of_pixels_across = 1000;
//    double pixel_size_mas = 0.0005;
//
//    // Setting geometry
//    Vector3d origin = {0., 0., 0.};
//	Vector3d direction = {0., 0., 1.};
//    // 3.7 deg
//    double cone_half_angle = 0.064577;
//    double big_scale = 100*pc;
//    Cone geometry(origin, direction, cone_half_angle, 100);
//
//    // Setting BField
//    RandomScalarBFieldZ bfield(10, 1.0);
//
//    // Setting N_field
//    BKNFieldZ nfield(10000, 2, true);
//
//    // Setting V-field
//    ConstFlatVField vfield(2.25);
//
//    Jet bkjet(&geometry, &vfield, &bfield, &nfield);
//
//    auto image_size = std::make_pair(number_of_pixels_across, number_of_pixels_along);
//    auto pc_in_mas = mas_to_pc(redshift);
//    std::cout << "pc_in_mas " << pc_in_mas << std::endl;
//    auto pixel_size = pixel_size_mas*pc_in_mas*pc;
//    auto pix_solid_angle = pixel_solid_angle(pixel_size_mas, redshift);
//
//    ImagePlane imagePlane(image_size, pixel_size, pixel_size, los_angle);
//    std::cout << "Setting pixel size pc " << pixel_size/pc << std::endl;
//
//    // Setting observed frequency
//    double nu = 15.4;
//    nu *= 1E+09;
//    nu *= (1.+redshift);
//
//    Observation observation(&bkjet, &imagePlane, nu);
//
//    string step_type = "adaptive";
//    double tau_max = 10000.0;
//    double tau_n_min = 0.1;
//    double dt_max_pc = 0.001;
//    double dt_max = pc*dt_max_pc;
//    double tau_min_log10 = -10.0;
//    double tau_min = pow(10.,tau_min_log10);
//    int n = 100;
//    int n_tau_max = 2000;
//    std::string calculate = "I";
//
//
//    observation.run(n, tau_max, dt_max, tau_min, step_type, calculate,
//                    n_tau_max, tau_n_min, tau_max);
//    string value = "tau";
//    auto image = observation.getImage(value);
//    std::fstream fs;
//    std::string file_tau = "map_tau.txt";
//    fs.open(file_tau, std::ios::out | std::ios::app);
//
//    if (fs.is_open())
//    {
//        write_2dvector(fs, image);
//        fs.close();
//    }
//
//    value = "I";
//    image = observation.getImage(value);
//    std::string file_i = "map_i.txt";
//    fs.open(file_i, std::ios::out | std::ios::app);
//    double scale = 1E-23*(1.+redshift)*(1.+redshift)*(1.+redshift)/pix_solid_angle;
//
//    if (fs.is_open())
//    {
//        // Scaling to Jy
//        write_2dvector(fs, image, scale);
//        fs.close();
//    }
//
//    value = "l";
//    image = observation.getImage(value);
//    std::string file_length = "map_l.txt";
//    fs.open(file_length, std::ios::out | std::ios::app);
//
//    if (fs.is_open()) {
//        write_2dvector(fs, image, pc);
//        fs.close();
//    }
//}


std::pair<vector<vector<double>>, vector<vector<double>>> run_analytical_params(double los_angle, double redshift,
        unsigned long int number_of_pixels_along, unsigned long int number_of_pixels_across,
        double pixel_size_mas_start, double pixel_size_mas_stop, double cone_half_angle, double B_1, double m,
        double K_1, double n, double s, double Gamma, double nu_observed_ghz, double tau_max, bool central_vfield=true);

std::pair<vector<vector<double>>, vector<vector<double>>>
run_analytical_params(double los_angle, double redshift, unsigned long int number_of_pixels_along,
                      unsigned long int number_of_pixels_across, double pixel_size_mas_start,
                      double pixel_size_mas_stop, double cone_half_angle, double B_1, double m, double K_1,
                      double n, double s, double Gamma, double nu_observed_ghz, double tau_max, bool central_vfield) {

    // Setting geometry
    Vector3d origin = {0., 0., 0.};
    Vector3d direction = {0., 0., 1.};
    double big_scale = 100*pc;
    Cone geometry(origin, direction, cone_half_angle, 100);

    // Setting BField
    RandomScalarBFieldZ bfield(B_1, m);
    //ToroidalBField bfield(B_1, m, true);

    // Setting N_field
    BKNFieldZ nfield(K_1, n, true, s);

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

    auto result = std::make_pair(image_tau, image_i);
    return result;
}


//void check_logspace() {
//    std::vector<double> logspace;
//    logspace = pyLogspace(-1, 1, 10);
//    for(auto el: logspace) {
//        std::cout << el << std::endl;
//    }
//}


//void check_image_building() {
//
//    auto image_size = std::make_pair(6, 12);
//    double pixel_size_start = -1;
//    double pixel_size_stop = 1;
//
//    std::vector<std::vector<double>> pixel_sizes_;
//    std::vector<std::vector<double>> pixel_center_coordinates_along_;
//    std::vector<std::vector<double>> pixel_center_coordinates_across_;
//
//
//    // Create array of pixel sizes
//    auto pixel_sizes_along = pyLogspace(pixel_size_start, pixel_size_stop, image_size.second);
//
//    //for(auto el: pixel_sizes_along) {
//    //    std::cout << el << std::endl;
//    //}
//
//
//    pixel_sizes_.reserve(image_size.first);
//    for (int i=0; i < image_size.first; ++i) {
//        pixel_sizes_.push_back(pixel_sizes_along);
//    }
//
//    for (int i = 0; i < image_size.first; i++) {
//        for (int j = 0; j < image_size.second; j++) {
//            std::cout << pixel_sizes_[i][j] << " ";
//            std::cout << std::endl;
//
//        }
//    }
//
//
//    // Create array of pixel center coordinates (n_across, n_along) for direction ALONG the jet
//    std::vector<double> cumsum_along;
//    cumsum_along.reserve(pixel_sizes_along.size());
//    std::partial_sum(pixel_sizes_along.begin(), pixel_sizes_along.end(), cumsum_along.begin());
//    for(int i = 0; i<pixel_sizes_along.size(); ++i) {
//        cumsum_along[i] = cumsum_along[i] - 0.5*pixel_sizes_along[i];
//    }
//
//
//    for (int i=0; i < image_size.second; ++i) {
//        std::cout << cumsum_along[i] << std::endl;
//    }
//
//    pixel_center_coordinates_along_.reserve(image_size.first);
//    for (int i=0; i < image_size.first; ++i) {
//        pixel_center_coordinates_along_[i].reserve(image_size.second);
//        for (int j=0; j < image_size.second; ++j) {
//            pixel_center_coordinates_along_[i][j] = cumsum_along[j];
//        }
//    }
//
//    std::cout << "Coordinates along = " << std::endl;
//    for (int i=0; i < image_size.first; ++i) {
//        std::cout << "Print i=" << i << std::endl;
//        for (int j=0; j < image_size.second; ++j) {
//            std::cout << pixel_center_coordinates_along_[i][j] << " ";
//        }
//        std::cout << std::endl;
//    }
//    std::cout << std::endl;
//
//
//
//    // Create array of pixel center coordinates (n_across, n_along) for direction ACROSS the jet
//    std::vector<std::vector<double>> cumsum_across;
//    // In each of the n_along arrays will be cumsums of the pixel coordinates in transverse direction
//    cumsum_across.reserve(image_size.second);
//    std::vector<double> pixel_sizes_transverse(image_size.first/2, 1.0);
//
//
//    // Going along the jet with larger and larger pixel sizes
//    for (int i=0; i < image_size.second; ++i) {
//        // Array of the constant pixel sizes across  [j] the jet at given distance [i] from center
//        for (int j=0; j < image_size.first/2; ++j) {
//            pixel_sizes_transverse[j] = pixel_sizes_along[i]*pixel_sizes_transverse[j];
//        }
//        cumsum_across[i].reserve(image_size.first/2);
//        std::partial_sum(pixel_sizes_transverse.begin(), pixel_sizes_transverse.end(), cumsum_across[i].begin());
//
//        for(int k = 0; k<pixel_sizes_transverse.size(); ++k) {
//            cumsum_across[i][k] = cumsum_across[i][k] - 0.5*pixel_sizes_transverse[k];
//        }
//
//        // Get ready for next transverse slice
//        for (int j=0; j < image_size.first/2; ++j) {
//            pixel_sizes_transverse[j] = 1.0;
//        }
//    }
//
//    // Flip
//    std::vector<std::vector<double> > trans_vec(image_size.first/2, std::vector<double>(image_size.second));
//    for (int i = 0; i < image_size.second; i++)
//    {
//        for (int j = 0; j < image_size.first/2; j++)
//        {
//            trans_vec[j][i] = cumsum_across[i][j];
//        }
//    }
//
//    // Negative coordinates
//    std::vector<std::vector<double> > trans_vec_neg(image_size.first/2, std::vector<double>(image_size.second));
//    for (int i = 0; i < image_size.second; i++)
//    {
//        for (int j = 0; j < image_size.first/2; j++)
//        {
//            trans_vec_neg[j][i] = -trans_vec[j][i];
//        }
//    }
//
//    // Flip positive coordinates
//    std::vector<std::vector<double> > trans_vec_flip(image_size.first/2, std::vector<double>(image_size.second));
//    for (int i = 0; i < image_size.second; i++)
//    {
//        for (int j = 0; j < image_size.first/2; j++)
//        {
//            trans_vec_flip[j][i] = trans_vec[image_size.first/2-j-1][i];
//        }
//    }
//
//    // Concatenate flipped positive coordinates with negative
//    //std::vector<std::vector<double> > result;
//    pixel_center_coordinates_across_ = trans_vec_flip;
//    pixel_center_coordinates_across_.insert(pixel_center_coordinates_across_.end(), trans_vec_neg.begin(), trans_vec_neg.end() );
//
//
//    std::cout << "Coordinates across = " << std::endl;
//    for (int i=0; i < image_size.first; ++i) {
//        std::cout << "Print i=" << i << std::endl;
//        for (int j=0; j < image_size.second; ++j) {
//            std::cout << pixel_center_coordinates_across_[i][j] << " ";
//        }
//        std::cout << std::endl;
//    }
//
//
//    //// Print out what happened
//    //std::cout << "Across sizes" << std::endl;
//    //for (int i=0; i < image_size.second; ++i) {
//    //    std::cout << "Along slice #" << i << std::endl;
//    //
//    //    // Array of the constant pixel sizes across  [j] the jet at given distance [i] from center
//    //    for (int j = 0; j < image_size.first/2; ++j) {
//    //        std::cout << cumsum_across[i][j] << " ";
//    //    }
//    //    std::cout << std::endl;
//    //}
//
//    //std::cout << "Along sizes" << std::endl;
//    //for (int i=0; i < image_size.second; ++i) {
//    //    std::cout << "Across slice #" << i << std::endl;
//    //    // Array of the constant pixel sizes across  [i] the jet at any transverse distance from axes.
//    //    std::cout << cumsum_along[i] << " ";
//    //    }
//
//    //std::cout << "Along sizes" << std::endl;
//    //for (int i=0; i < image_size.first; ++i) {
//    //    std::cout << "Along slice #" << i << std::endl;
//    //
//    //    // Array of the constant pixel sizes across  [j] the jet at given distance [i] from center
//    //    for (int j = 0; j < image_size.second; ++j) {
//    //        std::cout << pixel_center_coordinates_along_[i][j] << " ";
//    //    }
//    //    std::cout << std::endl;
//    //}
//
//}


//void check_image_plane() {
//
//    double redshift = 0.1;
//    double pixel_size_mas_start = 0.01;
//    double pixel_size_mas_stop = 1;
//    auto pc_in_mas = mas_to_pc(redshift);
//    std::cout << "pc_in_mas " << pc_in_mas << std::endl;
//    auto lg_pixel_size_start = log10(pixel_size_mas_start*pc_in_mas*pc);
//    auto lg_pixel_size_stop = log10(pixel_size_mas_stop*pc_in_mas*pc);
//
//    unsigned long int number_of_pixels_across = 6;
//    unsigned long int number_of_pixels_along = 12;
//    double los_angle = 3.1415926/2.0;
//    auto image_size = std::make_pair(number_of_pixels_across, number_of_pixels_along);
//
//    //Image image(image_size, pixel_size_start, pixel_size_stop);
//    //std::vector<Pixel>& pixels = image.getPixels();
//
//    ImagePlane imagePlane(image_size, lg_pixel_size_start, lg_pixel_size_stop, los_angle);
//    std::vector<Pixel>& pixels = imagePlane.getPixels();
//    std::vector<Ray>& rays = imagePlane.getRays();
//}


//int main() {
//    //check_image_building();
//    //check_image_plane();
//}


int main() {
	auto t1 = Clock::now();
	std::clock_t start;
	start = std::clock();

    //run_on_simulations();
    //run_analytical();
    auto result = run_analytical_params(0.663225, 0.0165, 500, 60, 0.001, 0.1, 0.06832789, 0.1, 1.0, 10000., 2.0, 2.5,
        5.0, 15.35, 100000, false);
    std::cout << "DONE" << std::endl;
    auto image_tau = result.first;
    auto image_i = result.second;

    std::fstream fs;
    std::string file_tau = "map_tau_num_tangled_flat_u.txt";
    fs.open(file_tau, std::ios::out | std::ios::app);

    if (fs.is_open())
    {
        write_2dvector(fs, image_tau);
        fs.close();
    }

    std::string file_i = "map_i_num_tangled_flat_u.txt";
    fs.open(file_i, std::ios::out | std::ios::app);

    if (fs.is_open())
    {
        write_2dvector(fs, image_i);
        fs.close();
    }


	std::cout << "CPU Time: "
						<< (std::clock() - start) / (double) (CLOCKS_PER_SEC)
						<< " s" << std::endl;
	auto t2 = Clock::now();
	std::cout << "User time: "
						<< std::chrono::duration_cast<std::chrono::seconds>(
								t2 - t1).count()
						<< " s" << std::endl;
}