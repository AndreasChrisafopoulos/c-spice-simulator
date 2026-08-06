#include "solver_direct.h"
#include <unistd.h>
#include <stdio.h> 
#include <string.h>
#include <strings.h> 

/* In-place LU solve: matrix A is overwritten */
int solve_lu_gsl_inplace(gsl_matrix *A, gsl_vector *b, gsl_vector *x_out) 
{
    size_t n = A->size1;
    if (A->size1 != A->size2 || b->size != n || x_out->size != n) return -1;

    gsl_permutation *p = gsl_permutation_alloc(n);
    int sign;
    gsl_linalg_LU_decomp(A, p, &sign);    /* in-place LU */
    gsl_linalg_LU_solve(A, p, b, x_out);
    gsl_permutation_free(p);
    return 0;
}

/* In-place Cholesky solver for Ax = b (A must be symmetric positive definite). */
int solve_cholesky_gsl_inplace(gsl_matrix *A, const gsl_vector *b, gsl_vector *x_out) 
{
    size_t n = A->size1;
    if (A->size1 != A->size2 || b->size != n || x_out->size != n) return -1;

    /* Factorize A in-place: A <- L (lower triangular stored in lower part) */
    int status = gsl_linalg_cholesky_decomp(A);
    if (status != 0) {
        /* if not SPD, gsl returns e.g. GSL_EDOM; handle or propagate */
        return status;
    }

    /* Solve using the factorization (cholesky is treated as const here) */
    status = gsl_linalg_cholesky_solve(A, b, x_out);
    return status;
}

/* Solves Ax = b using GSL; tries Cholesky if enabled, otherwise falls back to LU. */
int solve_system_gsl(gsl_matrix *A, gsl_vector *b, gsl_vector *x_out, int *use_cholesky)
{
    gsl_error_handler_t *old_handler = gsl_set_error_handler_off();
    int status = 0;    

    if (*use_cholesky) {
        gsl_matrix *A_temp = gsl_matrix_alloc(A->size1, A->size2);
        gsl_matrix_memcpy(A_temp, A);
        status = solve_cholesky_gsl_inplace(A, b, x_out);
        if (status != 0) {
            *use_cholesky = 0;
            printf("Cholesky failed (matrix not SPD). Falling back to LU.\n");
            status = solve_lu_gsl_inplace(A_temp, b, x_out);
        }
    } else {
        status = solve_lu_gsl_inplace(A, b, x_out);
    }

    gsl_set_error_handler(old_handler);
    return status;
}

/* Performs forward substitution.
 * Solves L·y = b in-place, where L comes from either:
 * - Cholesky factorization (non-unit lower triangular), or
 * - LU factorization (unit lower triangular).
 */
int forward_substitution(const gsl_matrix *A, gsl_vector *b, int use_cholesky)
{
    size_t n = A->size1;

    double *Ad = A->data;
    size_t lda = A->tda;
    double *bd = b->data;

    if (use_cholesky)
    {
        // Cholesky: solve L * y = b
        // L is lower triangular, diag != 1
        for (size_t i = 0; i < n; i++)
        {
            double *Ai = Ad + i * lda;
            double sum = bd[i];

            // sum -= L[i][j] * y[j]
            for (size_t j = 0; j < i; j++)
                sum -= Ai[j] * bd[j];

            bd[i] = sum / Ai[i];
        }
    }
    else
    {
        // LU: solve L * y = b
        // L is unit lower triangular (diag = 1)
        for (size_t i = 0; i < n; i++)
        {
            double *Ai = Ad + i * lda;
            double sum = bd[i];

            for (size_t j = 0; j < i; j++)
                sum -= Ai[j] * bd[j];

            // Diagonal is 1 → no division
            bd[i] = sum;
        }
    }

    return 0;
}




/* Performs backward substitution.
 * Solves either:
 * - U·x = y for LU factorization, or
 * - Lᵀ·x = y for Cholesky factorization,
 * writing the result into x.
 */
int backward_substitution(const gsl_matrix *A, gsl_vector *y, gsl_vector *x, int use_cholesky)
{
    size_t n = A->size1;

    // copy y -> x
    memcpy(x->data, y->data, n * sizeof(double));

    double *Ad = A->data;
    size_t lda = A->tda;
    double *xd = x->data;

    if (!use_cholesky)
    {
        
        // LU case: solve U * x = y
        // U = upper triangular
        for (ssize_t i = n - 1; i >= 0; i--)
        {
            double *Ai = Ad + i * lda;  // row i
            double sum = xd[i];

            // sum -= U[i][j] * x[j]
            for (size_t j = i + 1; j < n; j++)
                sum -= Ai[j] * xd[j];

            double Uii = Ai[i];
            xd[i] = sum / Uii;
        }
    }
    else
    {
        
        // Cholesky case: solve Lᵀ * x = y
        // L is lower triangular
        // so U = Lᵀ is upper
        for (ssize_t i = n - 1; i >= 0; i--)
        {
            double sum = xd[i];

            // traverse column i of L (row-major gives L[j][i])
            for (size_t j = i + 1; j < n; j++)
            {
                double *Aj = Ad + j * lda;
                sum -= Aj[i] * xd[j];   // L[j][i]
            }

            double Lii = Ad[i * lda + i]; // diag element
            xd[i] = sum / Lii;
        }
    }

    return 0;
}


/* Custom LU solver (in-place), with partial pivoting. */
int solve_lu_custom_inplace(gsl_matrix *A, gsl_vector *b, gsl_vector *x)
{
    size_t n = A->size1;
    double *Ad = A->data;
    size_t tda = A->tda;
    #define Aij(i,j) Ad[(i)*(tda) + (j)]

    for (size_t k = 0; k + 1 < n; k++) {
        size_t pivot = k;
        double max_val = fabs(Aij(k,k));
        for (size_t i = k + 1; i < n; i++) {
            double val = fabs(Aij(i,k));
            if (val > max_val) { max_val = val; pivot = i; }
        }
        if (pivot != k) {
            // swap rows: swap memory blocks of length n (use tmp buffer)
            for (size_t j = 0; j < n; j++) {
                double tmp = Aij(k,j); Aij(k,j) = Aij(pivot,j); Aij(pivot,j) = tmp;
            }
            double ttmp = gsl_vector_get(b, k);
            gsl_vector_set(b, k, gsl_vector_get(b, pivot));
            gsl_vector_set(b, pivot, ttmp);
        }

        double Akk = Aij(k,k);
        if (fabs(Akk) < 1e-18) { return -3; }
        for (size_t i = k + 1; i < n; i++) {
            double factor = Aij(i,k) / Akk;
            Aij(i,k) = factor;
            double *Ai_ptr = &Ad[i*tda];
            double *Ak_ptr = &Ad[k*tda];
            for (size_t j = k + 1; j < n; j++) {
                Ai_ptr[j] = Ai_ptr[j] - factor * Ak_ptr[j];
            }
        }
    }
    #undef Aij

    // forward/backward using cblas_dtrsv as before (but ensure data contiguous)
    forward_substitution(A, b, 0);
    backward_substitution(A, b, x, 0);
    return 0;
}


/* Checks if matrix A is symmetric within tolerance tol. */
int is_symmetric(const gsl_matrix *A, double tol)
{
    size_t n = A->size1;
    const double *Ad = A->data;
    size_t tda = A->tda;

    for (size_t i = 0; i < n; i++) {
        const double *Ai = Ad + i * tda;
        for (size_t j = i + 1; j < n; j++) {
            const double *Aj = Ad + j * tda;
            double aij = Ai[j];
            double aji = Aj[i];
            if (fabs(aij - aji) > tol)
                return 0;
        }
    }
    return 1;
}

/* Custom in-place Cholesky (symmetric, positive definite). */
int solve_cholesky_custom_inplace(gsl_matrix *A, gsl_vector *b, gsl_vector *x)
{
    size_t n = A->size1;
    double *Ad = A->data;
    size_t tda = A->tda;

    // Check symmetry (fast pointer version)
    if (!is_symmetric(A, 1e-12)) {
        fprintf(stderr, "Matrix is not symmetric.\n");
        return -1;
    }

    // Cholesky factorization: A = L * Lᵀ
    for (size_t j = 0; j < n; j++) {

        double *Aj = Ad + j * tda;  // row j
        double sum = 0.0;

        // Compute diagonal: Ljj
        for (size_t k = 0; k < j; k++) {
            double Ljk = Aj[k];   // A[j][k]
            sum += Ljk * Ljk;
        }

        double Ajj = Aj[j] - sum;
        if (Ajj <= 0.0)
            return -2; // not PD

        double Ljj = sqrt(Ajj);
        Aj[j] = Ljj;

        // Compute below diagonal: A[i][j]
        for (size_t i = j + 1; i < n; i++) {
            double *Ai = Ad + i * tda;

            double sum2 = 0.0;
            for (size_t k = 0; k < j; k++) {
                sum2 += Ai[k] * Aj[k];   // Lik * Ljk
            }

            double Aij = Ai[j];
            Ai[j] = (Aij - sum2) / Ljj;
        }
    }

    // Solve L*y = b
    forward_substitution(A, b, 1);

    // Solve Lᵀ*x = y
    backward_substitution(A, b, x, 1);

    return 0;
}

/* Custom linear solver: Cholesky (if enabled) with fallback to custom LU. */
int solve_system_custom(gsl_matrix *A, gsl_vector *b, gsl_vector *x_out, int *use_cholesky)
{
    int status = 0;    

    if (*use_cholesky) {
        gsl_matrix *A_temp = gsl_matrix_calloc(A->size1, A->size2);
        gsl_matrix_memcpy(A_temp, A);
        status = solve_cholesky_custom_inplace(A, b, x_out);
        if (status != 0) {
            *use_cholesky = 0;
            printf("Cholesky failed (matrix not SPD). Falling back to LU.\n");
            status = solve_lu_custom_inplace(A_temp, b, x_out);
        }
    } else {
        status = solve_lu_custom_inplace(A, b, x_out);
    }

  
    return status;
}



