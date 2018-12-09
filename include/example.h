#ifndef JETSHOW_EXAMPLE_H
#define JETSHOW_EXAMPLE_H


#include <vector>


using std::vector;


std::pair<vector<vector<double>>, vector<vector<double>>> run_analytical_params(double los_angle, double redshift,
                                                                                unsigned long int number_of_pixels_along, unsigned long int number_of_pixels_across, double pixel_size_mas,
                                                                                double cone_half_angle, double B_1, double m, double K_1, double n, double Gamma, double nu_observed_ghz,
                                                                                double tau_max, bool central_vfield=true);



#endif //JETSHOW_EXAMPLE_H
