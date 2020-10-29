/*
 * support.hpp
 *
 *  Created on: Oct 29, 2020
 *      Author: d-w-h
 */

#ifndef SUPPORT_HPP_
#define SUPPORT_HPP_

#include <string>

double min(double a, double b);
double psi(double r);
void export_grid_data(std::string grid_data_file, int nx, int ny, int nt);
void export_solver_data(std::string file_prefix, double **f, int nx, int ny, int ts_counter);

#endif /* SUPPORT_HPP_ */
