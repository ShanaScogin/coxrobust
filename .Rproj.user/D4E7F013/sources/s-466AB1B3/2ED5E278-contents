/*
 * Copyright (C) 2006 Filip Borowicz
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, 
 * MA  02111-1307  USA
 */

/**
 * Calculates the value of the negative partial log-likelihood function
 * and the first and second derivative at parameter value beta
 *
 * @parma the beta vector
 * @parma the vector of times
 * @parma the vector of censoring indicators
 * @parma the covariate matrix z
 * @parma the number of individuals
 * @parma the number of covariates
 * @parma the value of the negative partial log-likelihood function
 * @parma the first derivative of the negative partial log-likelihood function
 * @parma the second derivative of the negative partial log-likelihood function
 *
 * @return void
 */

#include "coxrobust.h"

void ple(double *beta, double *time, int *status, double *covar, int *n_row,
         int *n_col, double *res, double *gradient, double *hessian) {

    int i, j, k;
    double zbeta, exp_zbeta, sum1, tmp;
    double *sum2, *tmp_mat;
    double **z, **sum3, **hess;
    
    sum2    = (double *)R_alloc(*n_col, sizeof(double));
    tmp_mat = (double *)R_alloc(*n_col * *n_col, sizeof(double));
    z       = dmatrix(covar, *n_row, *n_col);
    sum3    = dmatrix(tmp_mat, *n_col, *n_col);
    hess    = dmatrix(hessian, *n_col, *n_col);
    
    sum1 = 0;
    *res = 0;
    for (j=0; j<*n_col; j++) {
        sum2[j] = 0;
        gradient[j] = 0;
        for (k=0; k<=j; k++) {
            sum3[j][k] = 0;
            hess[j][k] = 0;
        }
    }
    
    for (i=*n_row-1; i>=0; i--) {
    
        zbeta = 0;
        for (j=0; j<*n_col; j++) {
            zbeta += beta[j] * z[j][i];
        }
        exp_zbeta = exp(zbeta);
    
        sum1 += exp_zbeta;
        for (j=0; j<*n_col; j++) {
            sum2[j] += z[j][i] * exp_zbeta;
            for (k=0; k<=j; k++) {
                sum3[j][k] += z[j][i] * z[k][i] * exp_zbeta;
            }
        }
    
        if ( status[i] != 0 ) {
    
            *res += log(sum1) - zbeta;
    
            for (j=0; j<*n_col; j++) {
    
                tmp = sum2[j]/sum1;
                gradient[j] += tmp - z[j][i];
    
                for (k=0; k<=j; k++) {
                    hess[j][k] += (sum3[j][k] - tmp*sum2[k])/sum1;
                }
            }
    
        }
    
    }
    
    for (j=1; j<*n_col; j++) {
        for (k=0; k<j; k++) {
            hess[k][j] = hess[j][k];
        }
    }

}
