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
 * Robustly estimates the baseline cumulative hazard
 *
 * @parma exp(beta'z)
 * @parma the vector of times
 * @parma the vector of censoring indicators
 * @parma the exp(prev_beta'z) (prev_beta - beta obatined through previous step)
 * @parma the M value
 * @parma the number of individuals
 * @parma the type of weighting function
 * @parma the robustly estimated baseline cumulative hazard
 *
 * @return void
 */

#include "coxrobust.h"

void lambda(double *exp_zbeta, double *time, int *status,
            double *prev_exp_zbeta, double *M, int *n_row,
            int *a_type, double *lmb) {

    int i, k;
    double a_ii, sum;
    WEIGHT_FUNCTION A;
    
    A = get_weight_function(*a_type);
    
    for (i=0; i<*n_row; i++) {
    
        if ( status[i] != 0 ) {
    
            a_ii = A(time[i], prev_exp_zbeta[i], *M);
    
            if ( a_ii > 0 ) {
    
                sum = 0;
                for (k=i; k<*n_row; k++) {
                    sum += A(time[i], prev_exp_zbeta[k], *M) * exp_zbeta[k];
                }
    
                lmb[i] = (i == 0) ? a_ii/sum : lmb[i-1] + a_ii/sum;
    
            } else {
    
                lmb[i] = (i == 0) ? 0 : lmb[i-1];
    
            }
    
        } else {
    
            lmb[i] = (i == 0) ? 0 : lmb[i-1];
    
        }
    
    }

}
