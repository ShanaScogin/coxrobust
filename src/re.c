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
 * Calculates the value of the negative robust log-likelihood function
 * and the first and secod derivative at parameter value beta
 *
 * @parma the beta vector
 * @parma the vector of times
 * @parma the vector of censoring indicators
 * @parma the covariate matrix z
 * @parma the vector of exp(prev_beta'z) (prev_beta - beta obatined through previous step)
 * @parma the M value
 * @parma the number of individuals
 * @parma the number of covariates
 * @parma the type of weighting function
 * @parma the value of the negative robust log-likelihood function
 * @parma the first derivative of the negative robust log-likelihood function
 * @parma the second derivative of the negative robust log-likelihood function
 *
 * @return void
 */

#include "coxrobust.h"

void re(double *beta, double *time, int *status, double *covar,
        double *prev_exp_zbeta, double *M, int*n_row, int *n_col,
        int *a_type, double *res, double *gradient, double *hessian) {

    int i, ii, j, k;
    double sum1, a_ii, tmp;
    double *sum2, *zbeta, *tmp_mat;
    double *exp_zbeta;
    double **z, **sum3, **hess;
    WEIGHT_FUNCTION A;

    A = get_weight_function(*a_type);

    zbeta     = (double *)R_alloc(*n_row, sizeof(double));
    exp_zbeta = (double *)R_alloc(*n_row, sizeof(double));
    sum2      = (double *)R_alloc(*n_col, sizeof(double));
    tmp_mat   = (double *)R_alloc(*n_col * *n_col, sizeof(double));
    z 	      = dmatrix(covar, *n_row, *n_col);
    sum3 	  = dmatrix(tmp_mat, *n_col, *n_col);
    hess      = dmatrix(hessian, *n_col, *n_col);

    sum1 = 0;
    *res = 0;
    for  (j=0; j<*n_col; j++) {
        gradient[j] = 0;
    }

    for (i=*n_row-1; i>=0; i--) {
    
        zbeta[i] = 0;
        for (j=0; j<*n_col; j++) {
            zbeta[i] += beta[j] * z[j][i];
            sum2[j] = 0;
            for (k=0; k<=j; k++) {
                sum3[j][k] = 0;
            }
        }
    
        exp_zbeta[i] = exp(zbeta[i]);
    
        if ( status[i] != 0 ) {
    
            a_ii = A(time[i], prev_exp_zbeta[i], *M);
    
            sum1 = 0;
            for (ii=i; ii<*n_row; ii++) {
    
                tmp = A(time[i], prev_exp_zbeta[ii], *M)*exp_zbeta[ii];
                sum1 += tmp;
                for (j=0; j<*n_col; j++) {
                    sum2[j] += z[j][ii] * tmp;
                    for (k=0; k<=j; k++) {
                        sum3[j][k] += z[j][ii] * z[k][ii] * tmp;
                    }
                }
    
            }
    
            if ( sum1 == 0 ) {
                sum1 = 1.0;
            }
    
            *res += a_ii*(log(sum1) - zbeta[i]);
    
            for (j=0; j<*n_col; j++) {
                tmp = sum2[j]/sum1;
                gradient[j] += (tmp - z[j][i])*a_ii;
    
                for (k=0; k<=j; k++) {
                    hess[j][k] += ( (sum3[j][k] - tmp*sum2[k])/sum1 ) * a_ii;
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
