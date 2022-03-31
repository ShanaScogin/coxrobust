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

#ifndef COXROBUST_H
#define COXROBUST_H

#include <math.h>
#include "R.h"
#include "R_ext/Memory.h"

#define IDENT 0
#define LINEAR 1
#define QUADRATIC 2
#define EXPONENTIAL 3

typedef double (*WEIGHT_FUNCTION)(double, double, double);

double **dmatrix(double *array, int n_col, int n_row);

WEIGHT_FUNCTION get_weight_function(int type);

double ident(double time, double exp_zbeta, double M);

double linear(double time, double exp_zbeta, double M);

double quadratic(double time, double exp_zbeta, double M);

double exponential(double time, double exp_zbeta, double M);

void ple(double *beta, double *time, int *status, double *covar, int *n_row,
         int *n_col, double *res, double *gradient, double *hessian);

void re(double *beta, double *time, int *status, double *covar,
        double *prev_exp_zbeta, double *M, int*n_row, int *n_col,
        int *a_type, double *res, double *gradient, double *hessian);

void lambda(double *exp_zbeta, double *time, int *status,
            double *prev_exp_zbeta, double *M, int *n_row, int *a_type,
            double *lmb);

void lin(double *exp_zbeta, double *time, int *status, double *covar,
         double *prev_exp_zbeta, double *M, int *n_row, int *n_col,
         int *a_type, double *result);

#endif /* COXROBUST_H */
