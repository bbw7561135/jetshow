#include <string>
#include <iostream>
#include <fstream>
#include <boost/numeric/odeint.hpp>
#include "ImagePlane.h"
#include "Observation.h"
#include "Image.h"
#include "System.h"
#include "BField.h"
#include "VField.h"
#include "Ray.h"
#include "Cone.h"
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


void test_observations_rnd_bfield() {
	// Create a root
	pt::ptree root;

	Geometry* geometry;
	BField* bfield;
	VField* vfield;
	NField* nfield;

	// Load the json file in this ptree
	pt::read_json("../config.json", root);
	// Read values
	double los_angle = root.get<double>("observation.los_angle");
	double z = root.get<double>("observation.redshift");
	int number_of_pixels = root.get<int>("image.number_of_pixels");
	double pixel_size_mas = root.get<double>("image.pixel_size_mas");

	// Setting geometry
	Vector3d origin = {0., 0., 0.};
	Vector3d direction = {0., 0., 1.};
	std::string geotype = root.get<std::string>("jet.geometry.type");
	std::cout << "Geometry type " << geotype << std::endl;
	if (geotype == "cone") {
		// Create Cone with parameters
		double cone_angle = root.get<double>("jet.geometry.parameters.angle");
		double scale_pc = root.get<double>("jet.geometry.parameters.scale_pc");
		double scale = pc * scale_pc;
		geometry = new Cone(origin, direction, cone_angle, scale);
//				geometry = geometry_;
	}

	std::string btype = root.get<std::string>("jet.bfield.type");
	std::cout << "B-field type : " << btype << std::endl;
	if (btype == "radial_conical") {
		double b_1 = root.get<double>("jet.bfield.parameters.b_1");
		double n_b = root.get<double>("jet.bfield.parameters.n_b");
		bfield = new RadialConicalBField(b_1, n_b);
	} else if (btype == "spiral_conical") {
		double b_1 = root.get<double>("jet.bfield.parameters.b_1");
		double pitch_angle = root.get<double>("jet.bfield.parameters.pitch_angle");
		bfield = new SpiralConicalBField(b_1, pitch_angle);
	};

	// Create ``Cells`` instance
//	int N = 1000;
//	RandomCellsInSphere cells(N, 2.0);
//	double cone_angle = root.get<double>("jet.geometry.parameters.angle");
//	double scale_pc = root.get<double>("jet.geometry.parameters.scale_pc");
//	cells.setGeometry(scale_pc*pc, cone_angle);
//	cells.create();
	// Create ``RandomBField`` instance
	double rnd_fraction = root.get<double>("jet.bfield.parameters.random_fraction");
//	RandomCellsBField rnd_bfield(&cells, bfield, rnd_fraction);
	unsigned int seed = 123;
	if (rnd_fraction > 0) {
		bfield = new RandomPointBField(bfield, rnd_fraction, seed);
	}

	std::string vtype = root.get<std::string>("jet.vfield.type");
	std::cout << "Velocity type : " << vtype << std::endl;
	if (vtype == "const_central") {
		double gamma = root.get<double>("jet.vfield.parameters.gamma");
		vfield = new ConstCentralVField(gamma);
	}

	std::string ntype = root.get<std::string>("jet.nfield.type");
	std::cout << "Density type : " << ntype << std::endl;
	if (ntype == "bk") {
		double n_1 = root.get<double>("jet.nfield.parameters.n_1");
		double n_n = root.get<double>("jet.nfield.parameters.n_n");
		nfield = new BKNField(n_1, n_n);
	}


	Jet bkjet(geometry, vfield, bfield, nfield);

	auto image_size = std::make_pair(number_of_pixels, number_of_pixels);
	auto pc_in_mas = mas_to_pc(z);
//		auto cm_in_mas = pc * pc_in_mas;
	auto pixel_size = pixel_size_mas*0.25*pc_in_mas*pc;
	auto pix_solid_angle = pixel_solid_angle(pixel_size_mas, z);

	ImagePlane imagePlane(image_size, pixel_size, pixel_size, los_angle);

	double nu = root.get<double>("observation.frequency_ghz");
	nu *= 1E+09;
	// Redshifting to SBH frame
	nu *= (1.+z);
	Observation observation(&bkjet, &imagePlane, nu);
	string step_type = root.get<string>("integration.step_type");
	double tau_max = root.get<double>("integration.parameters.tau_max");
	double tau_n_min = root.get<double>("integration.parameters.tau_n_min");
	double dt_max_pc = root.get<double>("integration.parameters.dl_max_pc");
	double dt_max = pc*dt_max_pc;
	double tau_min_log10 = root.get<double>("integration.parameters.log10_tau_min");
	double tau_min = pow(10.,tau_min_log10);
	int n = root.get<int>("integration.parameters.n");
	int n_tau_max = root.get<int>("integration.parameters.n_tau_max");

	std::cout << "Integrating using max. opt. depth = " << tau_max << std::endl;
	std::cout << "Integrating using min. lg(opt.depth) = " << tau_min_log10 << std::endl;
	std::cout << "Integrating using max. step [pc] = " << dt_max_pc << std::endl;
	std::cout << "Integrating using default number of steps = " << n << std::endl;


	observation.run(n, tau_max, dt_max, tau_min, step_type, n_tau_max, tau_n_min,
	                tau_max);

	string value = "tau";
	auto image = observation.getImage(value);
	std::fstream fs;
	std::string file_tau = root.get<std::string>("output.file_tau");
	fs.open(file_tau, std::ios::out | std::ios::app);

	if (fs.is_open())
	{
		write_2dvector(fs, image);
		fs.close();
	}

	value = "I";
	image = observation.getImage(value);
	std::string file_i = root.get<std::string>("output.file_i");
	fs.open(file_i, std::ios::out | std::ios::app);
	double scale = 1E-23/pix_solid_angle;
	std::cout << "Scaling Stokes I by " << scale << std::endl;

	if (fs.is_open())
	{
		// Scaling to Jy
		write_2dvector(fs, image, scale);
		// Just to show how it can be used
		// write_2dvector(std::cout, image);
		fs.close();
	}

	value = "l";
	image = observation.getImage(value);
	std::string file_length = root.get<std::string>("output.file_length");
	fs.open(file_length, std::ios::out | std::ios::app);

	if (fs.is_open())
	{
		write_2dvector(fs, image, pc);
		fs.close();
	}
}


void test_observations_full() {
	// Create a root
	pt::ptree root;

	Geometry* geometry;
	BField* bfield;
	VField* vfield;
	NField* nfield;

	// Load the json file in this ptree
	pt::read_json("../config.json", root);
	// Read values
	double los_angle = root.get<double>("observation.los_angle");
	double z = root.get<double>("observation.redshift");
	int number_of_pixels = root.get<int>("image.number_of_pixels");
	double pixel_size_mas = root.get<double>("image.pixel_size_mas");

	// Setting geometry
	Vector3d origin = {0., 0., 0.};
	Vector3d direction = {0., 0., 1.};
	std::string geotype = root.get<std::string>("jet.geometry.type");
	std::cout << "Geometry type " << geotype << std::endl;
	if (geotype == "cone") {
		// Create Cone with parameters
		double cone_angle = root.get<double>("jet.geometry.parameters.angle");
		double scale_pc = root.get<double>("jet.geometry.parameters.scale_pc");
		double scale = pc * scale_pc;
		geometry = new Cone(origin, direction, cone_angle, scale);
//				geometry = geometry_;
	}

	std::string btype = root.get<std::string>("jet.bfield.type");
	std::cout << "B-field type : " << btype << std::endl;
	if (btype == "radial_conical") {
		double b_1 = root.get<double>("jet.bfield.parameters.b_1");
		double n_b = root.get<double>("jet.bfield.parameters.n_b");
		bfield = new RadialConicalBField(b_1, n_b);
	} else if (btype == "spiral_conical") {
		double b_1 = root.get<double>("jet.bfield.parameters.b_1");
		double pitch_angle = root.get<double>("jet.bfield.parameters.pitch_angle");
		bfield = new SpiralConicalBField(b_1, pitch_angle);
	};

	// Create ``Cells`` instance
//	int N = 1000;
//	RandomCellsInSphere cells(N, 2.0);
//	double cone_angle = root.get<double>("jet.geometry.parameters.angle");
//	double scale_pc = root.get<double>("jet.geometry.parameters.scale_pc");
//	cells.setGeometry(scale_pc*pc, cone_angle);
//	cells.create();
	// Create ``RandomBField`` instance
	double rnd_fraction = root.get<double>("jet.bfield.parameters.random_fraction");
//	RandomCellsBField rnd_bfield(&cells, bfield, rnd_fraction);
	unsigned int seed = 123;
	if (rnd_fraction > 0) {
		bfield = new RandomPointBField(bfield, rnd_fraction, seed);
	}

	std::string vtype = root.get<std::string>("jet.vfield.type");
	std::cout << "Velocity type : " << vtype << std::endl;
	if (vtype == "const_central") {
		double gamma = root.get<double>("jet.vfield.parameters.gamma");
		vfield = new ConstCentralVField(gamma);
	}

	std::string ntype = root.get<std::string>("jet.nfield.type");
	std::cout << "Density type : " << ntype << std::endl;
	if (ntype == "bk") {
		double n_1 = root.get<double>("jet.nfield.parameters.n_1");
		double n_n = root.get<double>("jet.nfield.parameters.n_n");
		nfield = new BKNField(n_1, n_n);
	}


	Jet bkjet(geometry, vfield, bfield, nfield);

	auto image_size = std::make_pair(number_of_pixels, number_of_pixels);
	auto pc_in_mas = mas_to_pc(z);
//		auto cm_in_mas = pc * pc_in_mas;
	auto pixel_size = pixel_size_mas*0.25*pc_in_mas*pc;
	auto pix_solid_angle = pixel_solid_angle(pixel_size_mas, z);

	ImagePlane imagePlane(image_size, pixel_size, pixel_size, los_angle);

	double nu = root.get<double>("observation.frequency_ghz");
	nu *= 1E+09;
	// Redshifting to SBH frame
	nu *= (1.+z);
	Observation observation(&bkjet, &imagePlane, nu);
	string step_type = root.get<string>("integration.step_type");
	double tau_max = root.get<double>("integration.parameters.tau_max");
	double tau_n_min = root.get<double>("integration.parameters.tau_n_min");
	double dt_max_pc = root.get<double>("integration.parameters.dl_max_pc");
	double dt_max = pc*dt_max_pc;
	double tau_min_log10 = root.get<double>("integration.parameters.log10_tau_min");
	double tau_min = pow(10.,tau_min_log10);
	int n = root.get<int>("integration.parameters.n");
	int n_tau_max = root.get<int>("integration.parameters.n_tau_max");

	std::cout << "Integrating using max. opt. depth = " << tau_max << std::endl;
	std::cout << "Integrating using min. lg(opt.depth) = " << tau_min_log10 << std::endl;
	std::cout << "Integrating using max. step [pc] = " << dt_max_pc << std::endl;
	std::cout << "Integrating using default number of steps = " << n << std::endl;


	observation.run(n, tau_max, dt_max, tau_min, step_type, n_tau_max, tau_n_min,
	                tau_max);

	string value = "tau";
	auto image = observation.getImage(value);
	std::fstream fs;
	std::string file_tau = root.get<std::string>("output.file_tau");
	fs.open(file_tau, std::ios::out | std::ios::app);

	if (fs.is_open())
	{
		write_2dvector(fs, image);
		fs.close();
	}

	value = "I";
	image = observation.getImage(value);
	std::string file_i = root.get<std::string>("output.file_i");
	fs.open(file_i, std::ios::out | std::ios::app);
	double scale = 1E-23/pix_solid_angle;
	std::cout << "Scaling Stokes I by " << scale << std::endl;

	if (fs.is_open())
	{
		// Scaling to Jy
		write_2dvector(fs, image, scale);
		// Just to show how it can be used
		// write_2dvector(std::cout, image);
		fs.close();
	}

	value = "Q";
	image = observation.getImage(value);
	std::string file_q = root.get<std::string>("output.file_q");
	fs.open(file_q, std::ios::out | std::ios::app);

	if (fs.is_open())
	{
		// Scaling to Jy
		write_2dvector(fs, image, scale);
		// Just to show how it can be used
		// write_2dvector(std::cout, image);
		fs.close();
	}

	value = "U";
	image = observation.getImage(value);
	std::string file_u = root.get<std::string>("output.file_u");
	fs.open(file_u, std::ios::out | std::ios::app);

	if (fs.is_open())
	{
		// Scaling to Jy
		write_2dvector(fs, image, scale);
		// Just to show how it can be used
		// write_2dvector(std::cout, image);
		fs.close();
	}

	value = "V";
	image = observation.getImage(value);
	std::string file_v = root.get<std::string>("output.file_v");
	fs.open(file_v, std::ios::out | std::ios::app);

	if (fs.is_open())
	{
		// Scaling to Jy
		write_2dvector(fs, image, scale);
		// Just to show how it can be used
		// write_2dvector(std::cout, image);
		fs.close();
	}

	value = "l";
	image = observation.getImage(value);
	std::string file_length = root.get<std::string>("output.file_length");
	fs.open(file_length, std::ios::out | std::ios::app);

	if (fs.is_open())
	{
		write_2dvector(fs, image, pc);
		fs.close();
	}
}


void test_stripe() {
	// Create a root
	pt::ptree root;

	Geometry* geometry;
	BField* bfield;
	VField* vfield;
	NField* nfield;

	// Load the json file in this ptree
	pt::read_json("../config.json", root);
	// Read values
	double los_angle = root.get<double>("observation.los_angle");
	double z = root.get<double>("observation.redshift");
	int number_of_pixels = root.get<int>("image.number_of_pixels");
	double pixel_size_mas = root.get<double>("image.pixel_size_mas");

	// Setting geometry
	Vector3d origin = {0., 0., 0.};
	Vector3d direction = {0., 0., 1.};
	std::string geotype = root.get<std::string>("jet.geometry.type");
	std::cout << "Geometry type " << geotype << std::endl;
	if (geotype == "cone") {
		// Create Cone with parameters
		double cone_angle = root.get<double>("jet.geometry.parameters.angle");
		double scale_pc = root.get<double>("jet.geometry.parameters.scale_pc");
		double scale = pc * scale_pc;
		geometry = new Cone(origin, direction, cone_angle, scale);
//				geometry = geometry_;
	}

	std::string btype = root.get<std::string>("jet.bfield.type");
	std::cout << "B-field type : " << btype << std::endl;
	if (btype == "radial_conical") {
		double b_1 = root.get<double>("jet.bfield.parameters.b_1");
		double n_b = root.get<double>("jet.bfield.parameters.n_b");
		bfield = new RadialConicalBField(b_1, n_b);
	} else if (btype == "spiral_conical") {
		double b_1 = root.get<double>("jet.bfield.parameters.b_1");
		double pitch_angle = root.get<double>("jet.bfield.parameters.pitch_angle");
		bfield = new SpiralConicalBField(b_1, pitch_angle);
	};

	// Create ``Cells`` instance
//	int N = 1000;
//	RandomCellsInSphere cells(N, 2.0);
//	double cone_angle = root.get<double>("jet.geometry.parameters.angle");
//	double scale_pc = root.get<double>("jet.geometry.parameters.scale_pc");
//	cells.setGeometry(scale_pc*pc, cone_angle);
//	cells.create();
	// Create ``RandomBField`` instance
	double rnd_fraction = root.get<double>("jet.bfield.parameters.random_fraction");
//	RandomCellsBField rnd_bfield(&cells, bfield, rnd_fraction);
	unsigned int seed = 123;
	if (rnd_fraction > 0) {
		bfield = new RandomPointBField(bfield, rnd_fraction, seed);
	}

	std::string vtype = root.get<std::string>("jet.vfield.type");
	std::cout << "Velocity type : " << vtype << std::endl;
	if (vtype == "const_central") {
		double gamma = root.get<double>("jet.vfield.parameters.gamma");
		vfield = new ConstCentralVField(gamma);
	}

	std::string ntype = root.get<std::string>("jet.nfield.type");
	std::cout << "Density type : " << ntype << std::endl;
	if (ntype == "bk") {
		double n_1 = root.get<double>("jet.nfield.parameters.n_1");
		double n_n = root.get<double>("jet.nfield.parameters.n_n");
		nfield = new BKNField(n_1, n_n);
	}


	Jet bkjet(geometry, vfield, bfield, nfield);

	auto image_size = std::make_pair(number_of_pixels, number_of_pixels);
	auto pc_in_mas = mas_to_pc(z);
//		auto cm_in_mas = pc * pc_in_mas;
	auto pixel_size = pixel_size_mas*0.25*pc_in_mas*pc;
	auto pix_solid_angle = pixel_solid_angle(pixel_size_mas, z);

	ImagePlane imagePlane(image_size, pixel_size, pixel_size, los_angle);

	double nu = root.get<double>("observation.frequency_ghz");
	nu *= 1E+09;
	// Redshifting to SBH frame
	nu *= (1.+z);
	Observation observation(&bkjet, &imagePlane, nu);
	string step_type = root.get<string>("integration.step_type");
	double tau_max = root.get<double>("integration.parameters.tau_max");
	double tau_n_min = root.get<double>("integration.parameters.tau_n_min");
	double dt_max_pc = root.get<double>("integration.parameters.dl_max_pc");
	double dt_max = pc*dt_max_pc;
	double tau_min_log10 = root.get<double>("integration.parameters.log10_tau_min");
	double tau_min = pow(10.,tau_min_log10);
	int n = root.get<int>("integration.parameters.n");
	int n_tau_max = root.get<int>("integration.parameters.n_tau_max");

	std::cout << "Integrating using max. opt. depth = " << tau_max << std::endl;
	std::cout << "Integrating using min. lg(opt.depth) = " << tau_min_log10 << std::endl;
	std::cout << "Integrating using max. step [pc] = " << dt_max_pc << std::endl;
	std::cout << "Integrating using default number of steps = " << n << std::endl;

	observation.run_stripe(n, tau_max, tau_min);
//	observation.run(n, tau_max, dt_max, tau_min, step_type, n_tau_max, tau_n_min,
//	                tau_max);

	string value = "tau";
	auto image = observation.getImage(value);
	std::fstream fs;
	std::string file_tau = root.get<std::string>("output.file_tau");
	fs.open(file_tau, std::ios::out | std::ios::app);

	if (fs.is_open())
	{
		write_2dvector(fs, image);
		fs.close();
	}
}


int main() {
	auto t1 = Clock::now();
	std::clock_t start;
	start = std::clock();

//	test_observations_rnd_bfield();
//	test_observations_full();
	test_stripe();

	std::cout << "CPU Time: "
						<< (std::clock() - start) / (double) (CLOCKS_PER_SEC)
						<< " s" << std::endl;
	auto t2 = Clock::now();
	std::cout << "User time: "
						<< std::chrono::duration_cast<std::chrono::seconds>(
								t2 - t1).count()
						<< " s" << std::endl;
}