#include <string>
#include <iostream>
#include <fstream>
#include <boost/numeric/odeint.hpp>
#include <cnpy.h>
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


using Eigen::Vector3d;
using Eigen::Matrix3Xd;
using std::vector;
using std::pair;
using namespace boost::numeric::odeint;
typedef std::chrono::high_resolution_clock Clock;
namespace pt = boost::property_tree;
namespace ph = std::placeholders;


//void test_observations_rnd_bfield() {
//	// Create a root
//	pt::ptree root;
//
//	Geometry* geometry;
//	BField* bfield;
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
//	// Create ``RandomBField`` instance
//	double rnd_fraction = root.get<double>("jet.bfield.parameters.random_fraction");
////	RandomCellsBField rnd_bfield(&cells, bfield, rnd_fraction);
//	unsigned int seed = 123;
//	if (rnd_fraction > 0) {
//		bfield = new RandomPointBField(bfield, rnd_fraction, seed);
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
//	BField* bfield;
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
//	// Create ``RandomBField`` instance
//	double rnd_fraction = root.get<double>("jet.bfield.parameters.random_fraction");
////	RandomCellsBField rnd_bfield(&cells, bfield, rnd_fraction);
//	unsigned int seed = 123;
//	if (rnd_fraction > 0) {
//		bfield = new RandomPointBField(bfield, rnd_fraction, seed);
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


//void test_intersection() {
//    Vector3d R0 = {0, 1, 1};
//    Vector3d Rd = {0, 1, 0};
//    Ray ray = Ray(R0, Rd);
//
//
//    std::list<double> result = intersection(R0, Rd, 1., 1., -1.);
//    std::cout << "Number of intersections = " << result.size() << std::endl;
//	for (double t : result) {
//		std::cout << ray.point(t) << '\n';
//	}
//
//}


// THIS IS MASTER FUNCTION
//void test_stripe() {
//	// Create a root
//	pt::ptree root;
//
//	Geometry* geometry;
////	RandomScalarBField* bfield;
//	BField* bfield;
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


void test_velocity() {
    VField* vfield;
    // Radius of parabaloid at z0=1pc
    double r0 = 0.1*pc;
    double gamma_axis0 = 30;
    double gamma_border0 = 1;
    // Z at which the speed on the axis is gamma0_axis
    double z_at_gamma0 = 1*pc;
    // Radius of parabaloid at z = 1pc
    double Rz0 = r0;
    vfield = new ShearedAccParabolicVField(gamma_axis0, gamma_border0, z_at_gamma0, Rz0);
    Vector3d point = {0, 0.025*pc, 1*pc};
    std::cout << vfield->v(point) << std::endl;

}

// Testing geometries, composite fields
void test_collimations() {
    Geometry* geometry;
//    BField* bfield;
    RandomScalarBFieldZ* bfield;
    VField* vfield;
    NField* nfield;

    double los_angle = pi/2;
    double z = 0.00436;
    int number_of_pixels = 1024;
    double pixel_size_mas = 0.005;

    // Setting geometry
    Vector3d origin = {0., 0., 0.};
    Vector3d direction = {0., 0., 1.};
//    double cone_half_angle = 0.1;
    double big_scale = 100*pc;
    // Radius of parabaloid at z0=1pc
    double r0 = 0.05*pc;
    // Distance where collimation stops (100 pix for 0.005mas pixel)
    double z0 = 0.05*pc;
//    geometry = new Cone(origin, direction, cone_half_angle, big_scale);
    geometry = new Parabaloid(origin, direction, r0, big_scale);
//    geometry = new ParabaloidCone(origin, direction, r0, z0, big_scale);

    double b_1 = 0.01;
    double n_b = 1.35;
//    bfield = new RadialConicalBField(b_1, n_b);
    bfield = new RandomScalarBFieldZ(b_1, n_b);

//    double gamma = 10;
//    double z_at_gamma = 1*pc;
//    vfield = new ConstCentralVField(gamma);
//    vfield = new AccParabolicVField(gamma, z_at_gamma);

//    double gamma_axis0 = 10;
//    double gamma_border0 = 1;
//    // Z at which the speed on the axis is gamma_0_axis
//    double z_at_gamma0 = 0.1*pc;
//    // Radius of parabaloid at z = 1pc
//    double Rz0 = r0;
//
//    vfield = new ShearedAccParabolicVField(gamma_axis0, gamma_border0, z_at_gamma0, Rz0);

//    double gamma0 = 10;
//    vfield = new AccParabolicConstConeVField(gamma0, r0, z0);


    double gamma_axis0 = 10;
    double gamma_border0 = 1;
    // Z at which the speed on the axis is gamma_0_axis
    double z_at_gamma0 = 0.05*pc;
    // Radius of parabaloid at z = 1pc
    double Rz0 = r0;

    vfield = new ShearedAccParabolicConstConeVField(gamma_axis0, gamma_border0, Rz0, z_at_gamma0);

    double n_1 = 10;
    double n_n = 2.0;
    nfield = new BKNField(n_1, n_n);

    Jet bkjet(geometry, vfield, bfield, nfield);

    auto image_size = std::make_pair(number_of_pixels, number_of_pixels);
    auto pc_in_mas = mas_to_pc(z);
    auto pixel_size = pixel_size_mas*pc_in_mas*pc;
    auto pix_solid_angle = pixel_solid_angle(pixel_size_mas, z);

    ImagePlane imagePlane(image_size, pixel_size, pixel_size, los_angle);
    std::cout << "Setting pixel size pc " << pixel_size/pc << std::endl;

    // Setting observed frequency
    double nu = 15.4;
    nu *= 1E+09;
    nu *= (1.+z);

    Observation observation(&bkjet, &imagePlane, nu);

    string step_type = "adaptive";
    double tau_max = 100.0;
    double tau_n_min = 0.1;
    double dt_max_pc = 0.001;
    double dt_max = pc*dt_max_pc;
    double tau_min_log10 = -8.0;
    double tau_min = pow(10.,tau_min_log10);
    int n = 100;
    int n_tau_max = 2000;
    std::string calculate = "I";


    observation.run(n, tau_max, dt_max, tau_min, step_type, calculate,
                    n_tau_max, tau_n_min, tau_max);
    string value = "tau";
    auto image = observation.getImage(value);
    std::fstream fs;
    std::string file_tau = "map_tau.txt";
    fs.open(file_tau, std::ios::out | std::ios::app);

    if (fs.is_open())
    {
        write_2dvector(fs, image);
        fs.close();
    }

    value = "I";
    image = observation.getImage(value);
    std::string file_i = "map_i.txt";
    fs.open(file_i, std::ios::out | std::ios::app);
    double scale = 1E-23*(1.+z)*(1.+z)*(1.+z)/pix_solid_angle;

    if (fs.is_open())
    {
        // Scaling to Jy
        write_2dvector(fs, image, scale);
        fs.close();
    }

    value = "l";
    image = observation.getImage(value);
    std::string file_length = "map_l.txt";
    fs.open(file_length, std::ios::out | std::ios::app);

    if (fs.is_open()) {
        write_2dvector(fs, image, pc);
        fs.close();
    }
}


void test_reading_npy() {
    cnpy::NpyArray arr = cnpy::npy_load("gamma.npy");
    double* loaded_data = arr.data<double>();
    size_t nrows = arr.shape[0];
    size_t ncols = arr.shape[1];
    std::cout << "nrows=" << nrows << ", ncols=" << ncols << std::endl;
}


void test_interpolation() {
    cnpy::NpyArray arr = cnpy::npy_load("gamma.npy");
    double* loaded_data = arr.data<double>();
    size_t nrows = arr.shape[0];
    size_t ncols = arr.shape[1];
    int n_points = nrows*ncols;
    std::cout << "nrows=" << nrows << ", ncols=" << ncols << std::endl;

    vector<double> r;
    vector<double> r_p;
    vector<double> gamma;

    for (int i=0; i<nrows; i++) {
        r.push_back(loaded_data[i*nrows]);
        r_p.push_back(loaded_data[i*nrows + 1]);
        gamma.push_back(loaded_data[i*nrows + 2]);
    }


//    double val = interp(10E20, 10E18);
//    std::cout << val << std::endl;
//
//
//    // Grid of r and r_p where to interpolate
//    std::vector<double> r_grid = linspace(0., 4*pow(10, 21), 2*n_points);
//    std::vector<double> r_p_grid = linspace(0., pow(10, 19), 2*n_points);
//
//    std::vector<double> interpolated_vals;
//
//    // Save values to file
//    cnpy::npy_save("interpolated.npy", &interpolated_vals[0], {r_grid.size(), r_p_grid.size()}, "w");
}


int main() {
	auto t1 = Clock::now();
	std::clock_t start;
	start = std::clock();

//	test_observations_rnd_bfield();
//	test_observations_full();
//	test_intersection();
	test_collimations();
//  test_reading_npy();
//    test_interpolation();
//	test_velocity();
//	test_stripe();

	std::cout << "CPU Time: "
						<< (std::clock() - start) / (double) (CLOCKS_PER_SEC)
						<< " s" << std::endl;
	auto t2 = Clock::now();
	std::cout << "User time: "
						<< std::chrono::duration_cast<std::chrono::seconds>(
								t2 - t1).count()
						<< " s" << std::endl;
}