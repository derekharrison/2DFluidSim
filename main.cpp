/*
 * main.cpp
 *
 *  Created on: Oct 24, 2020
 *      Author: d-w-h
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "fluid.hpp"
#include "memory_functions.hpp"
#include "support.hpp"
#include "user_types.hpp"

int main(int argc, char* argv[]) {

    p_params physical_params = {0};
    t_data time_data = {0};
    g_data grid_data = {0};
    s_data solver_data = {0};
    obj_data object_data = {0};

    /* Parameters */
    physical_params.pin = 1.0;    //Inlet pressure
    physical_params.pout = 0.0;   //Outlet pressure
    physical_params.rho = 3.0;    //Density
    physical_params.mu = 0.01;    //Viscosity

    time_data.ti = 0.0;           //Initial simulation time
    time_data.tf = 10.012;        //Final simulation time
    time_data.nt = 225;            //Number of timesteps

    grid_data.L = 1.0;            //Length of domain
    grid_data.W = 0.1;            //Width of domain
    grid_data.n_p_x = 60;         //Number of scalar nodes in x direction
    grid_data.n_p_y = 25;          //Number of scalar nodes in y direction. Note, this should be an odd number.

    object_data.x_loc = grid_data.L/3;
    object_data.y_loc = grid_data.W/2.5;
    object_data.R = grid_data.W/4;

    /* Measure execution time */
    clock_t start_time = clock();

    /* Allocate memory for solver data */
    int n_p_x, n_p_y, n_ux_x, n_ux_y, n_uy_x, n_uy_y;
    n_p_x = grid_data.n_p_x;
    n_p_y = grid_data.n_p_y;
    n_ux_x = n_p_x - 1;
    n_ux_y = n_p_y;
    n_uy_x = n_p_x;
    n_uy_y = n_p_y - 1;

    solver_data.p = array2D(n_p_x, n_p_y);
    solver_data.p_coord_x = array2D(n_p_x, n_p_y);
    solver_data.p_coord_y = array2D(n_p_x, n_p_y);
    solver_data.ux = array2D(n_ux_x, n_ux_y);
    solver_data.ux_coord_x = array2D(n_ux_x, n_ux_y);
    solver_data.ux_coord_y = array2D(n_ux_x, n_ux_y);
    solver_data.uy = array2D(n_uy_x, n_uy_y);
    solver_data.uy_coord_x = array2D(n_uy_x, n_uy_y);
    solver_data.uy_coord_y = array2D(n_uy_x, n_uy_y);

    /* Execute solver */
    fluid_solver(physical_params, time_data, grid_data, object_data, &solver_data);

    /* Print data */
    for(int i = 0; i < n_uy_x; ++i) {
        for(int j = 0; j < n_uy_y; ++j) {
            printf("uy[%i][%i]: %f, ", i, j, solver_data.uy[i][j]);
        }
        printf("\n");
    }

    printf("\n");
    for(int i = 0; i < n_ux_x; ++i) {
        for(int j = 0; j < n_ux_y; ++j) {
            printf("ux[%i][%i]: %f, ", i, j, solver_data.ux[i][j]);
        }
        printf("\n");
    }

    clock_t end_time = clock();
    double time_taken = (double) (end_time - start_time)/CLOCKS_PER_SEC;
    printf("time taken: %f\n", time_taken);
    return 0;
}
