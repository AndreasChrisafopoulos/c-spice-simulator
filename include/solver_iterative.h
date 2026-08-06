#ifndef SOLVE_ITERATIVE_H
#define SOLVE_ITERATIVE_H

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <math.h>        
#include <gsl/gsl_blas.h> 
#include <string.h>      
#include <stdio.h>


void buildPreconditioner(const gsl_matrix *A, gsl_vector *Minv);
int CG_solve(const gsl_matrix *A, const gsl_vector *b, gsl_vector *x,
             const gsl_vector *Minv, double itol, int maxIter);

int BiCG_solve(const gsl_matrix *A, const gsl_vector *b, gsl_vector *x,
               const gsl_vector *Minv, double itol, int maxIter);

int iterative_system_solver(gsl_matrix *A, gsl_vector *b, gsl_vector *x,
                            double itol, int maxIter, int spd_enable);

#endif
