/*
 * user_types.hpp
 *
 *  Created on: Oct 29, 2020
 *      Author: d-w-h
 */

#ifndef USER_TYPES_HPP_
#define USER_TYPES_HPP_

typedef struct physical_params {
	double pin;
	double pout;
	double rho;
	double mu;
} p_params;

typedef struct time_data {
	double ti;
	double tf;
	int nt;
} t_data;

typedef struct grid_data {
	double L;
	double W;
	int n_p_x;
	int n_p_y;
} g_data;

typedef struct solver_data {
    double **ux;
    double **ux_coord_x;
    double **ux_coord_y;
    double **uy;
    double **uy_coord_x;
    double **uy_coord_y;
    double **p;
    double **p_coord_x;
    double **p_coord_y;
} s_data;

typedef struct object_data {
	double x_loc;
	double y_loc;
	double R;
} obj_data;

#endif /* USER_TYPES_HPP_ */
