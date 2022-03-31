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

/*
 * @parma the exp(beta) vector
 * @parma the vector of times
 * @parma the vector of censoring indicators
 * @parma the covariate matrix z
 * @parma the vector of exp(prev_beta'z) (prev_beta - beta obatined through previous step)
 * @parma the M value
 * @parma the number of individuals
 * @parma the number of covariates
 * @parma the type of weighting function
 * @parma the result
 *
 * @return void
 */
#include "coxrobust.h"

void lin(double *exp_zbeta, double *time, int *status, double *covar,
         double *prev_exp_zbeta, double *M, int *n_row, int *n_col,
         int *a_type, double *result) {

    int i, j, k;
    double a_ii;
    double tmp;
    double sum1;
    double sum2;
    double *tmp1;
    double **z;
    double **res;
    WEIGHT_FUNCTION A;

    tmp1  = (double *)R_alloc(*n_row, sizeof(double));
    z     = dmatrix(covar, *n_row, *n_col);
    res   = dmatrix(result, *n_row, *n_col);

    A = get_weight_function(*a_type);

    for (j=0; j<*n_col; j++) {
        for (i=0; i<*n_row; i++) {
            
            if ( status[i] != 0 ) {

               a_ii = A(time[i], prev_exp_zbeta[i], *M);
               sum1 = 0;
               sum2 = 0;
               for (k=i; k<*n_row; k++) {
                   tmp = A(time[i], prev_exp_zbeta[k], *M) * exp_zbeta[k];
                   sum1 += tmp;
                   sum2 += z[j][k]*tmp;
               }

               if ( sum1 == 0 ) {
                  sum1 = 1.0;
               }

               res[j][i] = ( z[j][i] - sum2/sum1 )*a_ii;
               tmp1[i]   = ( (a_ii * exp_zbeta[i] ) / (sum1*sum1) ) * 
                           (sum2 - z[j][i]*sum1);

            } else {

               res[j][i] = 0;
               tmp1[i] = 0;

            }

        }

        for (i=0; i<*n_row; i++) {
            for (k=i; k<*n_row; k++) {
                res[j][i] -= A(time[i], prev_exp_zbeta[k], *M) * tmp1[k];
            }
        }
    
    }

}
