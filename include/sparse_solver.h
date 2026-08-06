#ifndef SPARSE_SOLVER_H
#define SPARSE_SOLVER_H

#include <gsl/gsl_vector.h>
#include "csparse.h"
#include <stdlib.h>        // malloc, free
#include <stdio.h>         // printf
#include <string.h>        // memcpy
#include <math.h>          // fabs, sqrt
#include <gsl/gsl_blas.h>

void sparse_lu_solve(css *S, csn *N, gsl_vector *b, gsl_vector *x, int n);
int sparse_lu_factorize(cs *A_csc, css **S, csn **N);
int sparse_cholesky_factorize(cs *A_csc, css **S, csn **N);
void sparse_cholesky_solve(css *S, csn *N,
                           gsl_vector *b, gsl_vector *x, int n);

int solve_system_sparse(cs *A_csc, gsl_vector *b,
                        gsl_vector *x_out, int n, int *use_cholesky, css *S, csn *N);

void sparse_Ax(const cs *A, const double *x, double *y);                        
void sparse_ATx(const cs *A, const double *x, double *y);                       
void buildPreconditioner_sparse(const cs *A, gsl_vector *Minv);
int CG_sparse_solve(const cs *A, const gsl_vector *b, gsl_vector *x,
                    const gsl_vector *Minv, double itol, int maxIter);
int BiCG_sparse(const cs *A, const gsl_vector *b, gsl_vector *x,
                const gsl_vector *Minv, double itol, int maxIter);
int iterative_system_solver_sparse(cs *A_csc, gsl_vector *b, gsl_vector *x,
                                   double itol, int maxIter, int spd_enable);

#endif
