#ifndef SOLVE_DIRECT_H
#define SOLVE_DIRECT_H

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include <math.h>   // fabs


int solve_lu_gsl_inplace(gsl_matrix *A, gsl_vector *b, gsl_vector *x_out);
int solve_cholesky_gsl_inplace(gsl_matrix *A, const gsl_vector *b, gsl_vector *x_out);
int solve_system_gsl(gsl_matrix *A, gsl_vector *b, gsl_vector *x_out, int *use_cholesky);

int solve_lu_custom_inplace(gsl_matrix *A, gsl_vector *b, gsl_vector *x);
int solve_cholesky_custom_inplace(gsl_matrix *A, gsl_vector *b, gsl_vector *x);
int is_symmetric(const gsl_matrix *A, double tol);
int solve_system_custom(gsl_matrix *A, gsl_vector *b, gsl_vector *x_out, int *use_cholesky);

int forward_substitution(const gsl_matrix *A, gsl_vector *b, int use_cholesky);
int backward_substitution(const gsl_matrix *A, gsl_vector *y, gsl_vector *x, int use_cholesky);

#endif
