/*
 * fluid.cpp
 *
 *  Created on: Oct 29, 2020
 *      Author: d-w-h
 */

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <math.h>
#include "memory_functions.hpp"
#include "support.hpp"
#include "user_types.hpp"

void fluid_solver(p_params physical_params, t_data time_data, g_data grid_data, obj_data object_data, s_data* solver_data) {
    double **ux, **uy, **p, L, W, rho, mu, pin, pout, ti, tf;
    int n_p_x, n_p_y, n_ux_x, n_ux_y, n_uy_x, n_uy_y, nt;

    /* Parameters */
    pin = physical_params.pin;
    pout = physical_params.pout;
    rho = physical_params.rho;
    mu = physical_params.mu;

    ti = time_data.ti;
    tf = time_data.tf;
    nt = time_data.nt;

    L = grid_data.L;
    W = grid_data.W;
    n_p_x = grid_data.n_p_x;
    n_p_y = grid_data.n_p_y;

    /* Object parameters */
    double x_loc, y_loc, R;
    x_loc = object_data.x_loc;
    y_loc = object_data.y_loc;
    R = object_data.R;

    /* Start calculations */
    double del_x, del_y, del_z, del_t;
    del_x = L / (n_p_x - 1);
    del_y = W / n_p_y;
    del_z = 1.0;
    del_t = (tf - ti) / nt;
    n_ux_x = n_p_x - 1;
    n_ux_y = n_p_y;
    n_uy_x = n_p_x;
    n_uy_y = n_p_y - 1;

    /* Allocate memory for data */
    double **ux_o, **uy_o, **p_o;
    double **x_coord, **y_coord;
    double **x_coord_uy, **y_coord_uy;
    double **x_p, **y_p, **ux_p, **uy_p;
    ux = array2D(n_ux_x, n_ux_y);
    x_coord = array2D(n_ux_x, n_ux_y);
    y_coord = array2D(n_ux_x, n_ux_y);
    uy = array2D(n_uy_x, n_uy_y);
    x_coord_uy = array2D(n_uy_x, n_uy_y);
    y_coord_uy = array2D(n_uy_x, n_uy_y);
    p = array2D(n_p_x, n_p_y);
    ux_o = array2D(n_ux_x, n_ux_y);
    uy_o = array2D(n_uy_x, n_uy_y);
    p_o = array2D(n_p_x, n_p_y);
    x_p = array2D(n_p_x, n_p_y);
    y_p = array2D(n_p_x, n_p_y);
    ux_p = array2D(n_p_x, n_p_y);
    uy_p = array2D(n_p_x, n_p_y);

    /* Initialize data */
    //Pressure
    for(int i = 0; i < n_p_x; ++i)
        for(int j = 0; j < n_p_y; ++j) {
            p[i][j] = (pout - pin) / (n_p_x - 1) * i + pin;
            p_o[i][j] = (pout - pin) / (n_p_x - 1) * i + pin;
            x_p[i][j] = i*del_x;
            y_p[i][j] = j*del_y + 0.5*del_y;
            ux_p[i][j] = 0.0;
            uy_p[i][j] = 0.0;
        }

    //X velocity
    for(int i = 0; i < n_ux_x; ++i)
        for(int j = 0; j < n_ux_y; ++j) {
            ux[i][j] = 0.0;
            ux_o[i][j] = 0.0;
            x_coord[i][j] = i*del_x + 0.5*del_x;
            y_coord[i][j] = j*del_y + 0.5*del_y;
        }

    //Y velocity
    for(int i = 0; i < n_uy_x; ++i)
        for(int j = 0; j < n_uy_y; ++j) {
            uy[i][j] = 0.0;
            uy_o[i][j] = 0.0;
            x_coord_uy[i][j] = i*del_x;
            y_coord_uy[i][j] = j*del_y + del_y;
        }

    /* Start Gauss-Seidel iterations */
    double max_ux, min_ux;
    double max_uy, min_uy;
    max_ux = -3e+8;
    min_ux = 3e+8;
    max_uy = -3e+8;
    min_uy = 3e+8;
    int max_iter = 1000;
    int ts_counter = 0;
    while(ts_counter < nt) {
        int iter_counter = 0;
        while(iter_counter < max_iter) {
            /* Ux nodes */
            //Upper left node ux
            double ux_w = 1.5 * ux_o[0][n_ux_y-1] - 0.5 * ux_o[1][n_ux_y-1];
            double ux_e = ux_o[0][n_ux_y-1];
            double ux_s = 0.5*(ux_o[0][n_ux_y-1] + ux_o[0][n_ux_y-2]);
            double uy_s = 0.5*(uy_o[0][n_uy_y-1] + uy_o[1][n_uy_y-1]);
            double ux_n = 0.0;
            double uy_n = 0.0;
            double Smx = del_y*del_z*rho*ux_w*ux_w - del_y*del_z*rho*ux_e*ux_e + del_x*del_z*rho*uy_s*ux_s - del_x*del_z*rho*uy_n*ux_n + del_y*del_z*(p_o[0][n_p_y-1] - p_o[1][n_p_y-1]);
            double duy_dx_s = (uy_o[1][n_uy_y-1] - uy_o[0][n_uy_y-1]) / del_x;
            double duy_dx_n = 0.0;
            ux[0][n_ux_y-1] = (Smx + del_x*del_z*-mu*duy_dx_s + mu*del_x*del_z*duy_dx_n + del_x*del_z*mu/del_y*ux[0][n_ux_y-2] + del_x*del_y*del_z*rho/del_t*ux_o[0][n_ux_y-1])/(del_x*del_z*mu/del_y + del_x*del_z*mu/(0.5*del_y) + del_x*del_y*del_z*rho/del_t);

            //Upper left central node ux
            double rw = 0.0;
            ux_w = ux_o[0][n_ux_y-1] + 0.5*(ux_o[1][n_ux_y-1] - ux_o[0][n_ux_y-1]);
            double re = (ux_o[1][n_ux_y-1] - ux_o[0][n_ux_y-1]) / (ux_o[2][n_ux_y-1] - ux_o[1][n_ux_y-1]);
            ux_e = ux_o[1][n_ux_y-1] + 0.5*psi(re)*(ux_o[2][n_ux_y-1] - ux_o[1][n_ux_y-1]);
            ux_s = 0.5*(ux_o[1][n_ux_y-1] + ux_o[1][n_ux_y-2]);
            uy_s = 0.5*(uy_o[1][n_uy_y-1] + uy_o[2][n_uy_y-1]);
            ux_n = 0.0;
            uy_n = 0.0;
            Smx = del_y*del_z*rho*ux_w*ux_w - del_y*del_z*rho*ux_e*ux_e + del_x*del_z*rho*uy_s*ux_s - del_x*del_z*rho*uy_n*ux_n + del_y*del_z*(p_o[1][n_p_y-1] - p_o[2][n_p_y-1]);
            duy_dx_s = (uy_o[2][n_uy_y-1] - uy_o[1][n_uy_y-1]) / del_x;
            duy_dx_n = 0.0;
            ux[1][n_ux_y-1] = (Smx + del_y*del_z*2*mu/del_x*ux[0][n_ux_y-1] + del_y*del_z*2*mu/del_x*ux[2][n_ux_y-1] + del_x*del_z*-mu*duy_dx_s + del_x*del_z*mu/del_y*ux[1][n_ux_y-2] + del_x*del_z*mu*duy_dx_n + del_x*del_y*del_z*rho/del_t*ux_o[1][n_ux_y-1]) / (del_y*del_z*2*mu/del_x + del_y*del_z*2*mu/del_x + del_x*del_z*mu/del_y + del_x*del_z*mu/(0.5*del_y) + del_x*del_y*del_z*rho/del_t);

            //Upper central nodes ux
            for(int i = 2; i < n_ux_x - 1; ++i) {
                rw = (ux_o[i-1][n_ux_y-1] - ux_o[i-2][n_ux_y-1]) / (ux_o[i+1][n_ux_y-1] - ux_o[i][n_ux_y-1]);
                ux_w = ux_o[i-1][n_ux_y-1] + 0.5*psi(rw)*(ux_o[i][n_ux_y-1] - ux_o[i-1][n_ux_y-1]);
                re = (ux_o[i][n_ux_y-1] - ux_o[i-1][n_ux_y-1]) / (ux_o[i+1][n_ux_y-1] - ux_o[i][n_ux_y-1]);
                ux_e = ux_o[i][n_ux_y-1] + 0.5*psi(re)*(ux_o[i+1][n_ux_y-1] - ux_o[i][n_ux_y-1]);
                ux_s = 0.5*(ux_o[i][n_ux_y-1] + ux_o[i][n_ux_y-2]);
                uy_s = 0.5*(uy_o[i][n_uy_y-1] + uy_o[i+1][n_uy_y-1]);
                ux_n = 0.0;
                uy_n = 0.0;
                Smx = del_y*del_z*rho*ux_w*ux_w - del_y*del_z*rho*ux_e*ux_e + del_x*del_z*rho*uy_s*ux_s - del_x*del_z*rho*uy_n*ux_n + del_y*del_z*(p_o[i][n_p_y-1] - p_o[i+1][n_p_y-1]);
                duy_dx_s = (uy_o[i+1][n_uy_y-1] - uy_o[i][n_uy_y-1]) / del_x;
                duy_dx_n = 0.0;
                ux[i][n_ux_y-1] = (Smx + del_y*del_z*2*mu/del_x*ux[i-1][n_ux_y-1] + del_y*del_z*2*mu/del_x*ux[i+1][n_ux_y-1] + del_x*del_z*-mu*duy_dx_s + del_x*del_z*mu/del_y*ux[i][n_ux_y-2] + del_x*del_z*mu*duy_dx_n + del_x*del_y*del_z*rho/del_t*ux_o[i][n_ux_y-1]) / (del_y*del_z*2*mu/del_x + del_y*del_z*2*mu/del_x + del_x*del_z*mu/del_y + del_x*del_z*mu/(0.5*del_y) + del_x*del_y*del_z*rho/del_t);
            }

            //Right most upper node ux
            rw = (ux_o[n_ux_x-1-1][n_ux_y-1] - ux_o[n_ux_x-1-2][n_ux_y-1]) / (ux_o[n_ux_x-1][n_ux_y-1] - ux_o[n_ux_x-2][n_ux_y-1]);
            ux_w = ux_o[n_ux_x-1-1][n_ux_y-1] + 0.5*psi(rw)*(ux_o[n_ux_x-1][n_ux_y-1] - ux_o[n_ux_x-1-1][n_ux_y-1]);
            re = 0.0;
            ux_e = ux_o[n_ux_x-1][n_ux_y-1];
            ux_s = 0.5*(ux_o[n_ux_x-1][n_ux_y-1] + ux_o[n_ux_x-1][n_ux_y-2]);
            uy_s = 0.5*(uy_o[n_ux_x-1][n_uy_y-1] + uy_o[n_ux_x-1+1][n_uy_y-1]);
            ux_n = 0.0;
            uy_n = 0.0;
            Smx = del_y*del_z*rho*ux_w*ux_w - del_y*del_z*rho*ux_e*ux_e + del_x*del_z*rho*uy_s*ux_s - del_x*del_z*rho*uy_n*ux_n + del_y*del_z*(p_o[n_ux_x-1][n_p_y-1] - p_o[n_ux_x-1+1][n_p_y-1]);
            duy_dx_s = (uy_o[n_ux_x-1+1][n_uy_y-1] - uy_o[n_ux_x-1][n_uy_y-1]) / del_x;
            duy_dx_n = 0.0;
            ux[n_ux_x-1][n_ux_y-1] = (Smx + del_y*del_z*2*mu/del_x*ux[n_ux_x-1-1][n_ux_y-1] + del_x*del_z*-mu*duy_dx_s + del_x*del_z*mu/del_y*ux[n_ux_x-1][n_ux_y-2] + del_x*del_z*mu*duy_dx_n + del_x*del_y*del_z*rho/del_t*ux_o[n_ux_x-1][n_ux_y-1]) / (del_y*del_z*2*mu/del_x + del_x*del_z*mu/del_y + del_x*del_z*mu/(0.5*del_y) + del_x*del_y*del_z*rho/del_t);

            //Middle inlet nodes ux
            for(int j = 1; j < n_ux_y - 1; ++j) {
                ux_w = 1.5 * ux_o[0][j] - 0.5 * ux_o[1][j];
                ux_e = ux_o[0][j];
                ux_s = 0.5*(ux_o[0][j] + ux_o[0][j-1]);
                uy_s = 0.5*(uy_o[0][j-1] + uy_o[1][j-1]);
                ux_n = 0.5*(ux_o[0][j] + ux_o[0][j+1]);
                uy_n = 0.5*(uy_o[0][j] + uy_o[1][j]);
                Smx = del_y*del_z*rho*ux_w*ux_w - del_y*del_z*rho*ux_e*ux_e + del_x*del_z*rho*uy_s*ux_s - del_x*del_z*rho*uy_n*ux_n + del_y*del_z*(p_o[0][j] - p_o[1][j]);
                duy_dx_s = (uy_o[1][j-1] - uy_o[0][j-1]) / del_x;
                duy_dx_n = (uy_o[1][j] - uy_o[0][j]) / del_x;
                ux[0][j] = (Smx + del_x*del_z*-mu*duy_dx_s + mu*del_x*del_z*duy_dx_n + del_x*del_z*mu/del_y*ux[0][j-1] + del_x*del_z*mu/del_y*ux[0][j+1] + del_x*del_y*del_z*rho/del_t*ux_o[0][j])/(del_x*del_z*mu/del_y + del_x*del_z*mu/del_y + del_x*del_y*del_z*rho/del_t);
            }

            //Lower left node
            ux_w = 1.5 * ux_o[0][0] - 0.5 * ux_o[1][0];
            ux_e = ux_o[0][0];
            ux_s = 0.0;
            uy_s = 0.0;
            ux_n = 0.5*(ux_o[0][0] + ux_o[0][1]);
            uy_n = 0.5*(uy_o[0][0] + uy_o[1][0]);
            Smx = del_y*del_z*rho*ux_w*ux_w - del_y*del_z*rho*ux_e*ux_e + del_x*del_z*rho*uy_s*ux_s - del_x*del_z*rho*uy_n*ux_n + del_y*del_z*(p_o[0][0] - p_o[1][0]);
            duy_dx_s = 0.0;
            duy_dx_n = (uy_o[1][0] - uy_o[0][0]) / del_x;
            ux[0][0] = (Smx + mu*del_x*del_z*duy_dx_n + del_x*del_z*mu/del_y*ux[0][1] + del_x*del_y*del_z*rho/del_t*ux_o[0][0])/(del_x*del_z*mu/(0.5*del_y) + del_x*del_z*mu/del_y + del_x*del_y*del_z*rho/del_t);

            //Lower left central node
            ux_w = ux_o[0][0] + 0.5*(ux_o[1][0] - ux_o[0][0]);
            re = (ux_o[1][0] - ux_o[0][0]) / (ux_o[2][0] - ux_o[1][0]);
            ux_e = ux_o[1][0] + 0.5*psi(re)*(ux_o[2][0] - ux_o[1][0]);
            ux_s = 0.0;
            uy_s = 0.0;
            ux_n = 0.5*(ux_o[1][0] + ux_o[1][1]);
            uy_n = 0.5*(uy_o[1][0] + uy_o[2][0]);
            Smx = del_y*del_z*rho*ux_w*ux_w - del_y*del_z*rho*ux_e*ux_e + del_x*del_z*rho*uy_s*ux_s - del_x*del_z*rho*uy_n*ux_n + del_y*del_z*(p_o[1][0] - p_o[2][0]);
            duy_dx_s = 0.0;
            duy_dx_n = (uy_o[2][0] - uy_o[1][0]) / del_x;
            ux[1][0] = (Smx + del_y*del_z*mu*2/del_x*ux[0][0] + del_y*del_z*mu*2/del_x*ux[2][0] + del_x*del_z*mu*duy_dx_n + del_x*del_z*mu/del_y*ux[1][1] + del_x*del_y*del_z*rho/del_t*ux_o[1][0])/(del_y*del_z*mu*2/del_x + del_y*del_z*mu*2/del_x + del_x*del_z*mu/(0.5*del_y) + del_x*del_z*mu/del_y + del_x*del_y*del_z*rho/del_t);

            //Lower central nodes
            for(int i = 2; i < n_ux_x - 1; ++i) {
                rw = (ux_o[i-1][0] - ux_o[i-2][0]) / (ux_o[i][0] - ux_o[i-1][0]);
                ux_w = ux_o[i-1][0] + 0.5*psi(rw)*(ux_o[i][0] - ux_o[i-1][0]);
                re = (ux_o[i][0] - ux_o[i-1][0]) / (ux_o[i+1][0] - ux_o[i][0]);
                ux_e = ux_o[i][0] + 0.5*psi(re)*(ux_o[i+1][0] - ux_o[i][0]);
                ux_s = 0.0;
                uy_s = 0.0;
                ux_n = 0.5*(ux_o[i][0] + ux_o[i][1]);
                uy_n = 0.5*(uy_o[i][0] + uy_o[i+1][0]);
                Smx = del_y*del_z*rho*ux_w*ux_w - del_y*del_z*rho*ux_e*ux_e + del_x*del_z*rho*uy_s*ux_s - del_x*del_z*rho*uy_n*ux_n + del_y*del_z*(p_o[i][0] - p_o[i+1][0]);
                duy_dx_s = 0.0;
                duy_dx_n = (uy_o[i+1][0] - uy_o[i][0]) / del_x;
                ux[i][0] = (Smx + del_y*del_z*mu*2/del_x*ux[i-1][0] + del_y*del_z*mu*2/del_x*ux[i+1][0] + del_x*del_z*mu*duy_dx_n + del_x*del_z*mu/del_y*ux[i][1] + del_x*del_y*del_z*rho/del_t*ux_o[i][0])/(del_y*del_z*mu*2/del_x + del_y*del_z*mu*2/del_x + del_x*del_z*mu/(0.5*del_y) + del_x*del_z*mu/del_y + del_x*del_y*del_z*rho/del_t);
            }

            //Right most lower node
            rw = (ux_o[n_ux_x-1-1][0] - ux_o[n_ux_x-1-2][0]) / (ux_o[n_ux_x-1][0] - ux_o[n_ux_x-1-1][0]);
            ux_w = ux_o[n_ux_x-1-1][0] + 0.5*psi(rw)*(ux_o[n_ux_x-1][0] - ux_o[n_ux_x-1-1][0]);
            re = 0.0;
            ux_e = ux_o[n_ux_x-1][0];
            ux_s = 0.0;
            uy_s = 0.0;
            ux_n = 0.5*(ux_o[n_ux_x-1][0] + ux_o[n_ux_x-1][1]);
            uy_n = 0.5*(uy_o[n_ux_x-1][0] + uy_o[n_ux_x-1+1][0]);
            Smx = del_y*del_z*rho*ux_w*ux_w - del_y*del_z*rho*ux_e*ux_e + del_x*del_z*rho*uy_s*ux_s - del_x*del_z*rho*uy_n*ux_n + del_y*del_z*(p_o[n_ux_x-1][0] - p_o[n_ux_x-1+1][0]);
            duy_dx_s = 0.0;
            duy_dx_n = (uy_o[n_ux_x-1+1][0] - uy_o[n_ux_x-1][0]) / del_x;
            ux[n_ux_x-1][0] = (Smx + del_y*del_z*mu*2/del_x*ux[n_ux_x-1-1][0] + del_x*del_z*mu*duy_dx_n + del_x*del_z*mu/del_y*ux[n_ux_x-1][1] + del_x*del_y*del_z*rho/del_t*ux_o[n_ux_x-1][0])/(del_y*del_z*mu*2/del_x + del_x*del_z*mu/(0.5*del_y) + del_x*del_z*mu/del_y + del_x*del_y*del_z*rho/del_t);

            //2nd central inlet nodes
            for(int j = 1; j < n_ux_y - 1; ++j) {
                ux_w = ux_o[0][j] + 0.5*(ux_o[1][j] - ux_o[0][j]);
                re = (ux_o[1][j] - ux_o[0][j]) / (ux_o[2][j] - ux_o[1][j]);
                ux_e = ux_o[1][j] + 0.5*psi(re)*(ux_o[2][j] - ux_o[1][j]);
                ux_s = 0.5*(ux_o[1][j] + ux_o[1][j-1]);
                uy_s = 0.5*(uy_o[1][j-1] + uy_o[2][j-1]);
                ux_n = 0.5*(ux_o[1][j] + ux_o[1][j+1]);
                uy_n = 0.5*(uy_o[1][j] + uy_o[2][j]);
                Smx = del_y*del_z*rho*ux_w*ux_w - del_y*del_z*rho*ux_e*ux_e + del_x*del_z*rho*uy_s*ux_s - del_x*del_z*rho*uy_n*ux_n + del_y*del_z*(p_o[1][j] - p_o[2][j]);
                duy_dx_s = (uy_o[2][j-1] - uy_o[1][j-1]) / del_x;
                duy_dx_n = (uy_o[2][j] - uy_o[1][j]) / del_x;
                ux[1][j] = (Smx + del_y*del_z*mu*2/del_x*ux[0][j] + del_y*del_z*mu*2/del_x*ux[2][j] + del_x*del_z*-mu*duy_dx_s + del_x*del_z*mu/del_y*ux[1][j-1] + del_x*del_z*mu*duy_dx_n + del_x*del_z*mu/del_y*ux[1][j+1] + del_x*del_y*del_z*rho/del_t*ux_o[1][j])/(del_y*del_z*mu*2/del_x + del_y*del_z*2*mu/del_x + del_x*del_z*mu/del_y + del_x*del_z*mu/del_y + del_x*del_y*del_z*rho/del_t);

            }

            //Central nodes
            for(int i = 2; i < n_ux_x - 1; ++i) {
                for(int j = 1; j < n_ux_y - 1; ++j) {
                	//Set velocities within object to zero
                    bool set_zero = ((x_coord[i][j] - x_loc)*(x_coord[i][j] - x_loc) + (y_coord[i][j] - y_loc)*(y_coord[i][j] - y_loc)) < (R*R);
                    if(set_zero) {
                        ux[i][j] = 0.0;
                        ux_o[i][j] = 0.0;
                    }
                    rw = (ux_o[i-1][j] - ux_o[i-2][j]) / (ux_o[i][j] - ux_o[i-1][j]);
                    ux_w = ux_o[i-1][j] + 0.5*psi(rw)*(ux_o[i][j] - ux_o[i-1][j]);
                    re = (ux_o[i][j] - ux_o[i-1][j]) / (ux_o[i+1][j] - ux_o[i][j]);
                    ux_e = ux_o[i][j] + 0.5*psi(re)*(ux_o[i+1][j] - ux_o[i][j]);
                    ux_s = 0.5*(ux_o[i][j] + ux_o[i][j-1]);
                    uy_s = 0.5*(uy_o[i][j-1] + uy_o[i+1][j-1]);
                    ux_n = 0.5*(ux_o[i][j] + ux_o[i][j+1]);
                    uy_n = 0.5*(uy_o[i][j] + uy_o[i+1][j]);
                    Smx = del_y*del_z*rho*ux_w*ux_w - del_y*del_z*rho*ux_e*ux_e + del_x*del_z*rho*uy_s*ux_s - del_x*del_z*rho*uy_n*ux_n + del_y*del_z*(p_o[i][j] - p_o[i+1][j]);
                    duy_dx_s = (uy_o[i+1][j-1] - uy_o[i][j-1]) / del_x;
                    duy_dx_n = (uy_o[i+1][j] - uy_o[i][j]) / del_x;
                    ux[i][j] = (Smx + del_y*del_z*mu*2/del_x*ux[i-1][j] + del_y*del_z*mu*2/del_x*ux[i+1][j] + del_x*del_z*-mu*duy_dx_s + del_x*del_z*mu/del_y*ux[i][j-1] + del_x*del_z*mu*duy_dx_n + del_x*del_z*mu/del_y*ux[i][j+1] + del_x*del_y*del_z*rho/del_t*ux_o[i][j])/(del_y*del_z*mu*2/del_x + del_y*del_z*2*mu/del_x + del_x*del_z*mu/del_y + del_x*del_z*mu/del_y + del_x*del_y*del_z*rho/del_t);
                    if(set_zero) {
                        ux[i][j] = 0.0;
                        ux_o[i][j] = 0.0;
                    }
                    //Calc min and max
                    if(ux[i][j] > max_ux) {
                    	max_ux = ux[i][j];
                    }
                    if(ux[i][j] < min_ux) {
                    	min_ux = ux[i][j];
                    }
                }
            }

            //Central outlet nodes
            for(int j = 1; j < n_ux_y - 1; ++j) {
                rw = (ux_o[n_ux_x-1-1][j] - ux_o[n_ux_x-1-2][j]) / (ux_o[n_ux_x-1][j] - ux_o[n_ux_x-1-1][j]);
                ux_w = ux_o[n_ux_x-1-1][j] + 0.5*psi(rw)*(ux_o[n_ux_x-1][j] - ux_o[n_ux_x-1-1][j]);
                re = 0.0;
                ux_e = ux_o[n_ux_x-1][j];
                ux_s = 0.5*(ux_o[n_ux_x-1][j] + ux_o[n_ux_x-1][j-1]);
                uy_s = 0.5*(uy_o[n_ux_x-1][j-1] + uy_o[n_ux_x-1+1][j-1]);
                ux_n = 0.5*(ux_o[n_ux_x-1][j] + ux_o[n_ux_x-1][j+1]);
                uy_n = 0.5*(uy_o[n_ux_x-1][j] + uy_o[n_ux_x-1+1][j]);
                Smx = del_y*del_z*rho*ux_w*ux_w - del_y*del_z*rho*ux_e*ux_e + del_x*del_z*rho*uy_s*ux_s - del_x*del_z*rho*uy_n*ux_n + del_y*del_z*(p_o[n_ux_x-1][j] - p_o[n_ux_x-1+1][j]);
                duy_dx_s = (uy_o[n_ux_x-1+1][j-1] - uy_o[n_ux_x-1][j-1]) / del_x;
                duy_dx_n = (uy_o[n_ux_x-1+1][j] - uy_o[n_ux_x-1][j]) / del_x;
                ux[n_ux_x-1][j] = (Smx + del_y*del_z*mu*2/del_x*ux[n_ux_x-1-1][j] + del_x*del_z*-mu*duy_dx_s + del_x*del_z*mu/del_y*ux[n_ux_x-1][j-1] + del_x*del_z*mu*duy_dx_n + del_x*del_z*mu/del_y*ux[n_ux_x-1][j+1] + del_x*del_y*del_z*rho/del_t*ux_o[n_ux_x-1][j])/(del_y*del_z*mu*2/del_x + del_x*del_z*mu/del_y + del_x*del_z*mu/del_y + del_x*del_y*del_z*rho/del_t);
            }

            /* Uy nodes */
            //Upper left central node
            uy_s = 0.5*(uy_o[1][n_uy_y-1] + uy_o[1][n_uy_y-2]);
            uy_n = 0.5*(uy_o[1][n_uy_y-1] + 0.0);
            ux_w = 0.5*(ux_o[0][n_ux_y-2] + ux_o[0][n_ux_y-1]);
            double uy_w = 0.5*(uy_o[0][n_uy_y-1] + uy_o[1][n_uy_y-1]);
            ux_e = 0.5*(ux_o[1][n_ux_y-2] + ux_o[1][n_ux_y-1]);
            double uy_e = 0.5*(uy_o[1][n_uy_y-1] + uy_o[2][n_uy_y-1]);
            double Smy = del_x*del_z*rho*uy_s*uy_s - del_x*del_z*rho*uy_n*uy_n + del_y*del_z*rho*ux_w*uy_w - del_y*del_z*rho*ux_e*uy_e + del_x*del_z*(p_o[1][n_p_y-2] - p_o[1][n_p_y-1]);
            double dux_dy_w = (ux_o[0][n_ux_y-1] - ux_o[0][n_ux_y-2]) / del_y;
            double dux_dy_e = (ux_o[1][n_ux_y-1] - ux_o[1][n_ux_y-2]) / del_y;
            uy[1][n_uy_y-1] = (Smy + del_x*del_z*mu*2/del_y*uy[1][n_uy_y-2] + del_y*del_z*-mu*dux_dy_w + del_y*del_z*mu*dux_dy_e + del_y*del_z*mu/del_x*uy[2][n_uy_y-1] + del_x*del_y*del_z*rho/del_t*uy_o[1][n_uy_y-1]) / (del_x*del_z*mu*2/del_y * 2 + del_y*del_z*mu/del_x * 2 + del_x*del_y*del_z*rho/del_t);

            //Upper central nodes
            for(int i = 2; i < n_uy_x - 1; ++i) {
                uy_s = 0.5*(uy_o[i][n_uy_y-1] + uy_o[i][n_uy_y-2]);
                uy_n = 0.5*(uy_o[i][n_uy_y-1] + 0.0);
                ux_w = 0.5*(ux_o[i-1][n_ux_y-2] + ux_o[i-1][n_ux_y-1]);
                uy_w = 0.5*(uy_o[i-1][n_uy_y-1] + uy_o[i][n_uy_y-1]);
                ux_e = 0.5*(ux_o[i][n_ux_y-2] + ux_o[i][n_ux_y-1]);
                uy_e = 0.5*(uy_o[i][n_uy_y-1] + uy_o[i+1][n_uy_y-1]);
                Smy = del_x*del_z*rho*uy_s*uy_s - del_x*del_z*rho*uy_n*uy_n + del_y*del_z*rho*ux_w*uy_w - del_y*del_z*rho*ux_e*uy_e + del_x*del_z*(p_o[i][n_p_y-2] - p_o[i][n_p_y-1]);
                dux_dy_w = (ux_o[i-1][n_ux_y-1] - ux_o[i-1][n_ux_y-2]) / del_y;
                dux_dy_e = (ux_o[i][n_ux_y-1] - ux_o[i][n_ux_y-2]) / del_y;
                uy[i][n_uy_y-1] = (Smy + del_x*del_z*mu*2/del_y*uy[i][n_uy_y-2] + del_y*del_z*-mu*dux_dy_w + del_y*del_z*mu*dux_dy_e + del_y*del_z*mu/del_x*uy[i-1][n_uy_y-1] + del_y*del_z*mu/del_x*uy[i+1][n_uy_y-1] + del_x*del_y*del_z*rho/del_t*uy_o[i][n_uy_y-1]) / (del_x*del_z*mu*2/del_y * 2 + del_y*del_z*mu/del_x * 2 + del_x*del_y*del_z*rho/del_t);

            }

            //Lower left central node
            uy_s = 0.5*(0.0 + uy_o[1][0]);
            uy_n = 0.5*(uy_o[1][0] + uy_o[1][1]);
            ux_w = 0.5*(ux_o[0][0] + ux_o[0][1]);
            uy_w = 0.5*(uy_o[0][0] + uy_o[1][0]);
            ux_e = 0.5*(ux_o[1][0] + ux_o[1][1]);
            uy_e = 0.5*(uy_o[1][0] + uy_o[2][0]);
            Smy = del_x*del_z*rho*uy_s*uy_s - del_x*del_z*rho*uy_n*uy_n + del_y*del_z*rho*ux_w*uy_w - del_y*del_z*rho*ux_e*uy_e + del_x*del_z*(p_o[1][0] - p_o[1][1]);
            dux_dy_w = (ux_o[0][1] - ux_o[0][0]) / del_y;
            dux_dy_e = (ux_o[1][1] - ux_o[1][0]) / del_y;
            uy[1][0] = (Smy + del_x*del_z*mu*2/del_y*uy[1][1] + del_y*del_z*-mu*dux_dy_w + del_y*del_z*mu*dux_dy_e + del_y*del_z*mu/del_x*uy[0][0] + del_y*del_z*mu/del_x*uy[2][0] + del_x*del_y*del_z*rho/del_t*uy_o[1][0]) / (del_x*del_z*mu*2/del_y * 2 + del_y*del_z*mu/del_x * 2 + del_x*del_y*del_z*rho/del_t);

            //Lower central nodes
            for(int i = 2; i < n_uy_x - 1; ++i) {
                uy_s = 0.5*(0.0 + uy_o[i][0]);
                uy_n = 0.5*(uy_o[i][0] + uy_o[i][1]);
                ux_w = 0.5*(ux_o[i-1][0] + ux_o[i-1][1]);
                uy_w = 0.5*(uy_o[i-1][0] + uy_o[i][0]);
                ux_e = 0.5*(ux_o[i][0] + ux_o[i][1]);
                uy_e = 0.5*(uy_o[i][0] + uy_o[i+1][0]);
                Smy = del_x*del_z*rho*uy_s*uy_s - del_x*del_z*rho*uy_n*uy_n + del_y*del_z*rho*ux_w*uy_w - del_y*del_z*rho*ux_e*uy_e + del_x*del_z*(p_o[i][0] - p_o[i][1]);
                dux_dy_w = (ux_o[i-1][1] - ux_o[i-1][0]) / del_y;
                dux_dy_e = (ux_o[i][1] - ux_o[i][0]) / del_y;
                uy[i][0] = (Smy + del_x*del_z*mu*2/del_y*uy[i][1] + del_y*del_z*-mu*dux_dy_w + del_y*del_z*mu*dux_dy_e + del_y*del_z*mu/del_x*uy[i-1][0] + del_y*del_z*mu/del_x*uy[i+1][0] + del_x*del_y*del_z*rho/del_t*uy_o[i][0]) / (del_x*del_z*mu*2/del_y * 2 + del_y*del_z*mu/del_x * 2 + del_x*del_y*del_z*rho/del_t);

            }

            //Central nodes
            for(int i = 1; i < n_uy_x - 1; ++i) {
                for(int j = 1; j < n_uy_y - 1; ++j) {
                	//Set velocities within object to zero
                    bool set_zero = ((x_coord_uy[i][j] - x_loc)*(x_coord_uy[i][j] - x_loc) + (y_coord_uy[i][j] - y_loc)*(y_coord_uy[i][j] - y_loc)) < (R*R);
                    if(set_zero) {
                        uy[i][j] = 0.0;
                        uy_o[i][j] = 0.0;
                    }
                    uy_s = 0.5*(uy_o[i][j-1] + uy_o[i][j]);
                    uy_n = 0.5*(uy_o[i][j] + uy_o[i][j+1]);
                    ux_w = 0.5*(ux_o[i-1][j] + ux_o[i-1][j+1]);
                    uy_w = 0.5*(uy_o[i-1][j] + uy_o[i][j]);
                    ux_e = 0.5*(ux_o[i][j] + ux_o[i][j+1]);
                    uy_e = 0.5*(uy_o[i][j] + uy_o[i+1][j]);
                    Smy = del_x*del_z*rho*uy_s*uy_s - del_x*del_z*rho*uy_n*uy_n + del_y*del_z*rho*ux_w*uy_w - del_y*del_z*rho*ux_e*uy_e + del_x*del_z*(p_o[i][j] - p_o[i][j+1]);
                    dux_dy_w = (ux_o[i-1][j+1] - ux_o[i-1][j]) / del_y;
                    dux_dy_e = (ux_o[i][j+1] - ux_o[i][j]) / del_y;
                    uy[i][j] = (Smy + del_x*del_z*mu*2/del_y*uy[i][j-1] + del_x*del_z*mu*2/del_y*uy[i][j+1] + del_y*del_z*-mu*dux_dy_w + del_y*del_z*mu*dux_dy_e + del_y*del_z*mu/del_x*uy[i-1][j] + del_y*del_z*mu/del_x*uy[i+1][j] + del_x*del_y*del_z*rho/del_t*uy_o[i][j]) / (del_x*del_z*mu*2/del_y * 2 + del_y*del_z*mu/del_x * 2 + del_x*del_y*del_z*rho/del_t);
                    if(set_zero) {
                        uy[i][j] = 0.0;
                    }
                    //Calc min and max
                    if(uy[i][j] > max_uy) {
                    	max_uy = uy[i][j];
                    }
                    if(uy[i][j] < min_uy) {
                    	min_uy = uy[i][j];
                    }
                }
            }

            //Outlet nodes uy
            for(int j = 0; j < n_uy_y; ++j) {
                uy[n_uy_x-1][j] = uy[n_uy_x-2][j];
            }

            iter_counter++;
        }

        //Update velocities for old timestep
        for(int i = 0; i < n_ux_x; ++i)
            for(int j = 0; j < n_ux_y; ++j) {
                ux_o[i][j] = ux[i][j];
            }

        for(int i = 0; i < n_uy_x; ++i)
            for(int j = 0; j < n_uy_y; ++j) {
                uy_o[i][j] = uy[i][j];
            }

        for(int i = 0; i < n_p_x; ++i)
            for(int j = 0; j < n_p_y; ++j) {
                p_o[i][j] = p[i][j];
            }

        //Calculate velocities at central scalar nodes
        for(int i = 1; i < n_p_x - 1; ++i) {
            for(int j = 1; j < n_p_y - 1; ++j) {
                ux_p[i][j] = 0.5*(ux[i-1][j] + ux[i][j]);
                uy_p[i][j] = 0.5*(uy[i][j-1] + uy[i][j]);
            }
        }

        /* Export data to file */
        //Export quiver data
        export_grid_data("grid_size_p.txt", n_p_x, n_p_y, nt);
        export_solver_data("data_ux_p_", ux_p, n_p_x, n_p_y, ts_counter);
        export_solver_data("data_uy_p_", uy_p, n_p_x, n_p_y, ts_counter);
        export_solver_data("data_x_p_", x_p, n_p_x, n_p_y, ts_counter);
        export_solver_data("data_y_p_", y_p, n_p_x, n_p_y, ts_counter);

        //Export Ux data
        export_grid_data("grid_size.txt", n_ux_x, n_ux_y, nt);
        export_solver_data("data_ux_", ux, n_ux_x, n_ux_y, ts_counter);
        export_solver_data("data_x_", x_coord, n_ux_x, n_ux_y, ts_counter);
        export_solver_data("data_y_", y_coord, n_ux_x, n_ux_y, ts_counter);
        export_max_and_min("max_and_min_ux.txt", max_ux, min_ux);

        //Export Uy data
        export_grid_data("grid_size_uy.txt", n_uy_x, n_uy_y, nt);
        export_solver_data("data_uy_", uy, n_uy_x, n_uy_y, ts_counter);
        export_solver_data("data_x_uy_", x_coord_uy, n_uy_x, n_uy_y, ts_counter);
        export_solver_data("data_y_uy_", y_coord_uy, n_uy_x, n_uy_y, ts_counter);
        export_max_and_min("max_and_min_uy.txt", max_uy, min_uy);

        ts_counter++;
    }

    /* Set results data */
    //Ux
    for(int i = 0; i < n_ux_x; ++i) {
        for(int j = 0; j < n_ux_y; ++j) {
            solver_data->ux[i][j] = ux[i][j];
            solver_data->ux_coord_x[i][j] = x_coord[i][j];
            solver_data->ux_coord_y[i][j] = y_coord[i][j];
        }
    }

    //Uy
    for(int i = 0; i < n_uy_x; ++i) {
        for(int j = 0; j < n_uy_y; ++j) {
            solver_data->uy[i][j] = uy[i][j];
            solver_data->uy_coord_x[i][j] = x_coord_uy[i][j];
            solver_data->uy_coord_y[i][j] = y_coord_uy[i][j];
        }
    }

    //P
    for(int i = 0; i < n_p_x; ++i) {
        for(int j = 0; j < n_p_y; ++j) {
            solver_data->p[i][j] = p[i][j];
            solver_data->p_coord_x[i][j] = x_p[i][j];
            solver_data->p_coord_y[i][j] = y_p[i][j];
        }
    }
}
