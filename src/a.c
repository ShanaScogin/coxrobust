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
 * Weighting function
 */
#include "coxrobust.h"

WEIGHT_FUNCTION get_weight_function(int type) {
    WEIGHT_FUNCTION f;
    
    switch (type) {
        case IDENT :
            f = ident;
            break;
        case LINEAR :
            f = linear;
            break;
        case QUADRATIC :
            f = quadratic;
            break;
        case EXPONENTIAL : 
            f = exponential;
            break;
        default : 
            f = ident;
    }

    return f;

}

double ident(double time, double exp_zbeta, double M) {
    return(1);
}

double linear(double time, double exp_zbeta, double M) {
    double a, tmp;
    tmp = time*exp_zbeta;

    if ( M > tmp ) {
        a = M - tmp;
    } else {
        a = 0;
    }
    return(a);
}

double quadratic(double time, double exp_zbeta, double M) {
    double a, tmp;
    tmp = M - time*exp_zbeta;

    if ( tmp > 0 ) {
        a = (tmp*tmp) / (M*M);
    } else {
        a = 0;
    }

    return(a);
}

double exponential(double time, double exp_zbeta, double M) {
    double a;
    a = exp( -fabs( (time*exp_zbeta) / M ) );
    return(a);
}
