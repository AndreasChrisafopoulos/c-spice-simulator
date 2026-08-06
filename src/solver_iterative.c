#include "solver_iterative.h"

/* Computes the Jacobi preconditioner (Minv) from matrix A.
Minv[i] = 1 / A[i][i], or 1 if the diagonal is zero.
Used by iterative solvers (CG/BiCG) to improve convergence.*/
void buildPreconditioner(const gsl_matrix *A, gsl_vector *Minv)
{
    int n = A->size1;

    for (int i = 0; i < n; i++) {
        double aii = gsl_matrix_get(A, i, i);

        if (fabs(aii) > 1e-20)
            gsl_vector_set(Minv, i, 1.0 / aii);
        else
            gsl_vector_set(Minv, i, 1.0);  // fallback
    }
}



// CG with Jacobi Preconditioning (Stable)
int CG_solve(
    const gsl_matrix *A,
    const gsl_vector *b,
    gsl_vector *x,
    const gsl_vector *Minv,
    double itol,
    int maxIter
){
    int n = b->size;

    gsl_vector *r  = gsl_vector_alloc(n);
    gsl_vector *z  = gsl_vector_alloc(n);
    gsl_vector *p  = gsl_vector_alloc(n);
    gsl_vector *Ap = gsl_vector_alloc(n);

    /* r = b - A*x */
    gsl_vector_memcpy(r, b);
    gsl_blas_dgemv(CblasNoTrans, -1.0, A, x, 1.0, r);

    /* z = M^{-1} r  precond*/
    for(int i=0;i<n;i++)
        gsl_vector_set(z,i, gsl_vector_get(r,i)*gsl_vector_get(Minv,i));

    gsl_vector_memcpy(p, z);

    double rho;
    gsl_blas_ddot(r, z, &rho); // ger rho

    double bnorm = gsl_blas_dnrm2(b);
    if (bnorm == 0) bnorm = 1.0;

    for(int k=0;k<maxIter;k++) //main loop
    {
        /* Ap = A*p */
        gsl_blas_dgemv(CblasNoTrans, 1.0, A, p, 0.0, Ap);

        double pAp;
        gsl_blas_ddot(p, Ap, &pAp);

        if (fabs(pAp) < 1e-30) {
            printf("CG breakdown p^T A p = 0\n");
            break;
        }

        double alpha = rho / pAp;

        /* x = x + alpha p */
        gsl_blas_daxpy(alpha, p, x); //renew x

        /* r = r - alpha Ap */
        gsl_blas_daxpy(-alpha, Ap, r); //renew r

        /* check convergence */
        double rnorm = gsl_blas_dnrm2(r);
        if (rnorm / bnorm < itol) {
            gsl_vector_free(r);
            gsl_vector_free(z);
            gsl_vector_free(p);
            gsl_vector_free(Ap);
            return k;
        }

        /* z = Minv r */
        for(int i=0;i<n;i++)
            gsl_vector_set(z,i, gsl_vector_get(r,i)*gsl_vector_get(Minv,i));

        double rho_new;
        gsl_blas_ddot(r, z, &rho_new);

        double beta = rho_new / rho;

        /* p = z + beta p */
        gsl_vector_scale(p, beta);
        gsl_vector_add(p, z);

        rho = rho_new;
    }

    gsl_vector_free(r);
    gsl_vector_free(z);
    gsl_vector_free(p);
    gsl_vector_free(Ap);

    return -1;
}

/* Bi-CG with Jacobi Preconditioning (stable) */
int BiCG_solve(
    const gsl_matrix *A,
    const gsl_vector *b,
    gsl_vector *x,
    const gsl_vector *Minv,
    double itol,
    int maxIter
){
    int n = b->size;

    gsl_vector *r    = gsl_vector_alloc(n);
    gsl_vector *rt   = gsl_vector_alloc(n);
    gsl_vector *z    = gsl_vector_alloc(n);
    gsl_vector *zt   = gsl_vector_alloc(n);
    gsl_vector *p    = gsl_vector_alloc(n);
    gsl_vector *pt   = gsl_vector_alloc(n);
    gsl_vector *q    = gsl_vector_alloc(n);
    gsl_vector *qt   = gsl_vector_alloc(n);

    /* r = b - A x */
    gsl_vector_memcpy(r, b);
    gsl_blas_dgemv(CblasNoTrans, -1.0, A, x, 1.0, r);

    /* r̃ = r */
    gsl_vector_memcpy(rt, r);

    double bnorm = gsl_blas_dnrm2(b);
    if (bnorm == 0) bnorm = 1.0;

    double rho_old = 0.0;

    for (int iter = 1; iter <= maxIter; iter++)
    {
        /* Solve M z = r (diagonal preconditioner) */
        for (int i=0;i<n;i++)
            gsl_vector_set(z, i, gsl_vector_get(r,i) * gsl_vector_get(Minv,i));

        /* Solve Mᵀ z̃ = r̃ (same for diagonal M) */
        for (int i=0;i<n;i++)
            gsl_vector_set(zt, i, gsl_vector_get(rt,i) * gsl_vector_get(Minv,i));

        /* rho = r̃ᵀ z */
        double rho;
        gsl_blas_ddot(rt, z, &rho);
        if (fabs(rho) < 1e-30) goto fail;

        if (iter == 1)
        {
            gsl_vector_memcpy(p, z);
            gsl_vector_memcpy(pt, zt);
        }
        else
        {
            double beta = rho / rho_old;

            gsl_vector_scale(p, beta);
            gsl_vector_add(p, z);

            gsl_vector_scale(pt, beta);
            gsl_vector_add(pt, zt);
        }

        rho_old = rho;

        /* q = A p */
        gsl_blas_dgemv(CblasNoTrans, 1.0, A, p, 0.0, q);

        /* q̃ = Aᵀ p̃ */
        gsl_blas_dgemv(CblasTrans, 1.0, A, pt, 0.0, qt);

        /* omega = p̃ᵀ q */
        double omega;
        gsl_blas_ddot(pt, q, &omega);
        if (fabs(omega) < 1e-30) goto fail;

        double alpha = rho / omega;

        /* x = x + α p */
        gsl_blas_daxpy(alpha, p, x);

        /* r = r − α q */
        gsl_blas_daxpy(-alpha, q, r);

        /* r̃ = r̃ − α q̃ */
        gsl_blas_daxpy(-alpha, qt, rt);

        /* convergence test */
        double rnorm = gsl_blas_dnrm2(r);
        if (rnorm / bnorm < itol)
        {
            gsl_vector_free(r);
            gsl_vector_free(rt);
            gsl_vector_free(z);
            gsl_vector_free(zt);
            gsl_vector_free(p);
            gsl_vector_free(pt);
            gsl_vector_free(q);
            gsl_vector_free(qt);
            return iter;
        }
    }

fail:
    gsl_vector_free(r);
    gsl_vector_free(rt);
    gsl_vector_free(z);
    gsl_vector_free(zt);
    gsl_vector_free(p);
    gsl_vector_free(pt);
    gsl_vector_free(q);
    gsl_vector_free(qt);
    return -1;
}

/* Iterative linear system solver (CG for SPD, BiCG otherwise) with Jacobi preconditioning. */
int iterative_system_solver(gsl_matrix *A, gsl_vector *b, gsl_vector *x, double itol, int maxIter, int spd_enable) {
    int n = A->size1;
    int iters;
    gsl_vector *Minv = gsl_vector_calloc(n);
    buildPreconditioner(A, Minv);

    if (spd_enable) 
        iters = CG_solve(A, b, x, Minv, itol, maxIter);
    else
        iters = BiCG_solve(A, b, x, Minv, itol, maxIter);

    return iters;
}