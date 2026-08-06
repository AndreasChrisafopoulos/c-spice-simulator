/* Sparse linear solvers (LU / Cholesky / CG / BiCG) for MNA systems using CSparse. */
#include "sparse_solver.h"

int sparse_lu_factorize(cs *A_csc, css **S, csn **N)
{

    *S = cs_sqr(2, A_csc, 0);
    *N = cs_lu(A_csc, *S, 1);

    if (!(*S) || !(*N))
        return -1;

    return 0;
}


void sparse_lu_solve(css *S, csn *N, gsl_vector *b, gsl_vector *x, int n)
{
    int i;
    double *bb = malloc(n * sizeof(double));
    double *y  = malloc(n * sizeof(double));

    // copy b -> bb 
    for (i = 0; i < n; i++)
        bb[i] = gsl_vector_get(b, i);

    // Pb 
    cs_ipvec(N->pinv, bb, y, n);

    // Ly = Pb 
    cs_lsolve(N->L, y);

    // Ux = y 
    cs_usolve(N->U, y);

    // x = Q*y 
    cs_ipvec(S->q, y, bb, n);

    // copy sol
    for (i = 0; i < n; i++)
        gsl_vector_set(x, i, bb[i]);

    free(bb);
    free(y);
}


int sparse_cholesky_factorize(cs *A_csc, css **S, csn **N)
{
   // cs *A_csc = cs_compress(A_triplet);
    //cs_dupl(A_csc);

    *S = cs_schol(1, A_csc);
    *N = cs_chol(A_csc, *S);

    //cs_spfree(A_csc);

    if (!(*S) || !(*N))
        return -1;   // FAILED

    return 0;        // OK
}


void sparse_cholesky_solve(css *S, csn *N, gsl_vector *b, gsl_vector *x, int n)
{
    int i;
    double *bb = malloc(n * sizeof(double));
    double *y  = malloc(n * sizeof(double));

    // copy b -> bb
    for (i = 0; i < n; i++)
        bb[i] = gsl_vector_get(b, i);

    // Pb
    cs_ipvec(S->pinv, bb, y, n);

    // Ly = Pb
    cs_lsolve(N->L, y);

    // Lᵀx = y
    cs_ltsolve(N->L, y);

    // x = Pᵀ y
    cs_pvec(S->pinv, y, bb, n);

    // copy solution
    for (i = 0; i < n; i++)
        gsl_vector_set(x, i, bb[i]);

    free(bb);
    free(y);
}

/* Sparse system solver: Cholesky if SPD, otherwise LU fallback */
int solve_system_sparse(cs *A_csc,
                        gsl_vector *b,
                        gsl_vector *x_out,
                        int n,
                        int *use_cholesky, css *S, csn *N)
{
    //css *S = NULL;
    //csn *N = NULL;
    int status=0;

    if (*use_cholesky) {
        //status = sparse_cholesky_factorize(A_csc, &S, &N);

        if (status == 0) {
            sparse_cholesky_solve(S, N, b, x_out, n);
        } else {
            printf("Sparse Cholesky failed. Falling back to LU.\n");
            *use_cholesky = 0;

            status = sparse_lu_factorize(A_csc, &S, &N);
            if (status != 0) return -1;

            sparse_lu_solve(S, N, b, x_out, n);
        }
    } else {
        //status = sparse_lu_factorize(A_csc, &S, &N);
        if (status != 0) return -1;

        sparse_lu_solve(S, N, b, x_out, n);
    }

   // cs_sfree(S);
    //cs_nfree(N);

    return 0;
}

/* y = A * x  (sparse CSC) */
void sparse_Ax(const cs *A, const double *x, double *y)
{
    int i, j, p;
    for (i = 0; i < A->m; i++) y[i] = 0.0;

    for (j = 0; j < A->n; j++)
        for (p = A->p[j]; p < A->p[j+1]; p++)
            y[A->i[p]] += A->x[p] * x[j];
}

/* y = Aᵀ * x  (sparse CSC) */
void sparse_ATx(const cs *A, const double *x, double *y)
{
    int  j, p;
    for (j = 0; j < A->n; j++) y[j] = 0.0;

    for (j = 0; j < A->n; j++)
        for (p = A->p[j]; p < A->p[j+1]; p++)
            y[j] += A->x[p] * x[A->i[p]];
}


/* Jacobi preconditioner for sparse CSC matrix */
void buildPreconditioner_sparse(const cs *A, gsl_vector *Minv)
{
    int n = A->n;
    int j, p;

    for (j = 0; j < n; j++)
        gsl_vector_set(Minv, j, 1.0);  // default

    for (j = 0; j < n; j++)
        for (p = A->p[j]; p < A->p[j+1]; p++)
            if (A->i[p] == j && fabs(A->x[p]) > 1e-20)
                gsl_vector_set(Minv, j, 1.0 / A->x[p]);
}

int CG_sparse_solve(
    const cs *A,              // sparse CSC matrix
    const gsl_vector *b,
    gsl_vector *x,
    const gsl_vector *Minv,
    double itol,
    int maxIter
){
    int n = b->size;

    double *r  = calloc(n, sizeof(double));
    double *z  = calloc(n, sizeof(double));
    double *p  = calloc(n, sizeof(double));
    double *Ap = calloc(n, sizeof(double));

    int i;
    double rho, rho_new, alpha, beta;

    // r = b - A x 
    sparse_Ax(A, x->data, r);
    for (i = 0; i < n; i++)
        r[i] = gsl_vector_get(b, i) - r[i];

    // z = M^{-1} r 
    for (i = 0; i < n; i++)
        z[i] = r[i] * gsl_vector_get(Minv, i);

    memcpy(p, z, n * sizeof(double));

    rho = 0.0;
    for (i = 0; i < n; i++)
        rho += r[i] * z[i];

    double bnorm = gsl_blas_dnrm2(b);
    if (bnorm == 0.0) bnorm = 1.0;

    for (int k = 0; k < maxIter; k++)
    {
        // Ap = A p 
        sparse_Ax(A, p, Ap);

        double pAp = 0.0;
        for (i = 0; i < n; i++)
            pAp += p[i] * Ap[i];

        if (fabs(pAp) < 1e-30) break;

        alpha = rho / pAp;

        for (i = 0; i < n; i++) {
            x->data[i] += alpha * p[i];
            r[i]       -= alpha * Ap[i];
        }

        double rnorm = 0.0;
        for (i = 0; i < n; i++)
            rnorm += r[i] * r[i];
        rnorm = sqrt(rnorm);

        if (rnorm / bnorm < itol) {
            free(r); free(z); free(p); free(Ap);
            return k;
        }

        for (i = 0; i < n; i++)
            z[i] = r[i] * gsl_vector_get(Minv, i);

        rho_new = 0.0;
        for (i = 0; i < n; i++)
            rho_new += r[i] * z[i];

        beta = rho_new / rho;

        for (i = 0; i < n; i++)
            p[i] = z[i] + beta * p[i];

        rho = rho_new;
    }

    free(r); free(z); free(p); free(Ap);
    return -1;
}

int BiCG_sparse(
    const cs *A,                  // sparse CSC matrix
    const gsl_vector *b,
    gsl_vector *x,
    const gsl_vector *Minv,        // Jacobi preconditioner (1/diag)
    double itol,
    int maxIter
){
    int n = b->size;

    double *r  = calloc(n, sizeof(double));
    double *rt = calloc(n, sizeof(double));
    double *z  = calloc(n, sizeof(double));
    double *zt = calloc(n, sizeof(double));
    double *p  = calloc(n, sizeof(double));
    double *pt = calloc(n, sizeof(double));
    double *q  = calloc(n, sizeof(double));
    double *qt = calloc(n, sizeof(double));

    // r = b - A x 
    sparse_Ax(A, x->data, r);
    for (int i = 0; i < n; i++)
        r[i] = gsl_vector_get(b, i) - r[i];

    // r~ = r 
    memcpy(rt, r, n * sizeof(double));

    double bnorm = gsl_blas_dnrm2(b);
    if (bnorm == 0.0) bnorm = 1.0;

    double rho_old = 0.0;

    for (int iter = 1; iter <= maxIter; iter++)
    {
        // z = M⁻¹ r , zt = M⁻¹ rt 
        for (int i = 0; i < n; i++) {
            z[i]  = gsl_vector_get(Minv, i) * r[i];
            zt[i] = gsl_vector_get(Minv, i) * rt[i];
        }

        // rho = rtᵀ z 
        double rho = 0.0;
        for (int i = 0; i < n; i++)
            rho += rt[i] * z[i];

        if (fabs(rho) < 1e-30) break;

        if (iter == 1) {
            memcpy(p,  z,  n * sizeof(double));
            memcpy(pt, zt, n * sizeof(double));
        } else {
            double beta = rho / rho_old;
            for (int i = 0; i < n; i++) {
                p[i]  = z[i]  + beta * p[i];
                pt[i] = zt[i] + beta * pt[i];
            }
        }

        rho_old = rho;

        // q = A p , qt = Aᵀ pt 
        sparse_Ax(A,  p,  q);
        sparse_ATx(A, pt, qt);

        double omega = 0.0;
        for (int i = 0; i < n; i++)
            omega += pt[i] * q[i];

        if (fabs(omega) < 1e-30) break;

        double alpha = rho / omega;

        // update x, r, rt 
        for (int i = 0; i < n; i++) {
            x->data[i] += alpha * p[i];
            r[i]       -= alpha * q[i];
            rt[i]      -= alpha * qt[i];
        }

        // convergence check 
        double rnorm = 0.0;
        for (int i = 0; i < n; i++)
            rnorm += r[i] * r[i];
        rnorm = sqrt(rnorm);

        if (rnorm / bnorm < itol) {
            free(r); free(rt); free(z); free(zt);
            free(p); free(pt); free(q); free(qt);
            return iter;
        }
    }

    free(r); free(rt); free(z); free(zt);
    free(p); free(pt); free(q); free(qt);
    return -1;
}

/* Iterative sparse solver: CG for SPD matrices, BiCG otherwise */
int iterative_system_solver_sparse(
    cs *A_csc,
    gsl_vector *b,
    gsl_vector *x,
    double itol,
    int maxIter,
    int spd_enable
){
    int n = A_csc->n;
    int iters;

    gsl_vector *Minv = gsl_vector_calloc(n);
    buildPreconditioner_sparse(A_csc, Minv);

    if (spd_enable)
    {
        iters = CG_sparse_solve(A_csc, b, x, Minv, itol, maxIter);

        if (iters < 0) {
            printf("Sparse CG failed → fallback to BiCG.\n");
            iters = BiCG_sparse(A_csc, b, x, Minv, itol, maxIter);
        }
    }
    else
    {
        iters = BiCG_sparse(A_csc, b, x, Minv, itol, maxIter);
    }

    gsl_vector_free(Minv);
    return iters;
}







