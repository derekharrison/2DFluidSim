/*
 * support.cpp
 *
 *  Created on: Oct 29, 2020
 *      Author: d-w-h
 */

#include <stdio.h>
#include <string.h>
#include <string>
#include <cstring>
#include <iostream>
#include <fstream>

double min(double a, double b) {
    double min;
    if(a > b) {
        min = b;
    }
    else {
        min = a;
    }

    return min;
}

double psi(double r) {
    double result;
    if(r > 0) {
        result = min(r, 1);
    }
    else {
        result = 0.0;
    }

    return result;
}

void export_grid_data(std::string grid_data_file, int nx, int ny, int nt) {
    std::ofstream myfile(grid_data_file);
    myfile << nx << " " << ny << " " << nt;
    myfile.close();
}

void export_solver_data(std::string file_prefix, double **f, int nx, int ny, int ts_counter) {
    std::string index = std::to_string(ts_counter);
    std::string file = file_prefix + index + ".txt";
    std::ofstream myfile(file);

    for(int i = 0; i < nx; ++i) {
        for(int j = 0; j < ny; ++j) {
            myfile << f[i][j] << " ";
        }
        myfile << "\n";
    }
    myfile.close();
}

void export_max_and_min(std::string file_name, double max, double min) {
    std::ofstream myfile(file_name);
    myfile << max << " " << min;
    myfile.close();
}
