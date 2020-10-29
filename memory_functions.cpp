/*
 * memory_functions.cpp
 *
 *  Created on: Oct 29, 2020
 *      Author: d-w-h
 */

double** array2D(int nx, int ny) {

    double** f = new double*[nx];
    for(int i = 0; i < nx; ++i) {
        f[i] = new double[ny];
    }

    return f;
}
