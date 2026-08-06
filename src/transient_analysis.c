#include "eval.h"
#include "solver_direct.h"
#include "solver_iterative.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <stdio.h>
#include <sys/stat.h>
#include "stamp.h"
#include "sparse_solver.h"

#define V(v,i) gsl_vector_get(v,i)

/**
 * Stamps transient (time-dependent) sources into the RHS vector b.
 *
 * Updates the right-hand side vector for the current time t by evaluating
 * all transient sources and adding their contributions:
 *
 *  - Voltage sources are stamped directly into the augmented rows
 *    corresponding to MNA voltage source equations.
 *  - Current sources are injected into the appropriate node equations
 *    via sweep stamping.
 *
 * Notes:
 *  - The matrix structure is NOT modified.
 *  - The vector b is assumed to be pre-initialized with DC contributions.
 *  - Only time-dependent increments are applied.
 */
void stamp_transient_sources(
    gsl_vector *b,
    TranSrc *tran_sources, int tran_src_number,
    Element *m2, Element *elements,
    int node_count,
    double t
){
    double v_stamp, i_stamp;

    for ( int i = 0 ; i < tran_src_number ; i++) {

        int idx = tran_sources[i].index;

        if (tran_sources[i].type) {                    // Voltage source
            TwoPortsSource *v = (TwoPortsSource *) m2[idx].data;
            v_stamp = eval_source(v, t) ;
            
            int row = node_count - 1 + idx ; 
            gsl_vector_set(b, row, v_stamp);

        } else {                                       // Current source
            TwoPortsSource *I = (TwoPortsSource *) elements[idx].data;
            i_stamp = eval_source(I, t)- I->value;
            sweep_stamp(b, I->ports[0], I->ports[1], i_stamp, 0, 0);
        }
    }
}


/**
 * Performs transient analysis of the circuit using Modified Nodal Analysis.
 *
 * Computes the time-domain response of the circuit by solving a sequence of
 * linear systems derived from the MNA formulation and the selected time
 * integration method.
 *
 * Workflow:
 *  - Solves the DC operating point to obtain the initial state.
 *  - Builds the transient system matrix A = Gdc + αC.
 *  - Advances the solution in time using backward Euler or trapezoidal rule.
 *  - At each time step:
 *      * Stamps time-dependent sources into the RHS vector.
 *      * Forms the transient RHS using the previous solution.
 *      * Solves the linear system using direct or iterative methods.
 *      * Stores selected node voltages to output files.
 *
 * Notes:
 *  - The matrices Gdc and C are modified in-place during setup.
 *  - The RHS vector b_dc is treated as the DC contribution baseline.
 *  - Matrix A is rebuilt at each time step for solver compatibility.
 *  - Output files and gnuplot scripts are generated automatically.
 */
int run_transient_analysis(
    gsl_matrix *Gdc, gsl_matrix *C, 
    gsl_vector *b_dc,
    Element *m2, Element *elements,
    Node *nodes, int node_count,
    TRAN_Analysis *tran, TranSrc *tran_sources, int tran_src_number,
    int use_iterative, int *use_cholesky, int use_custom, double itol,
    char outdir[256]
){
    int n = b_dc->size;
    double scale = (tran->method ? 2/tran->time_step : 1/tran->time_step);

    gsl_matrix *Acopy = gsl_matrix_alloc(n,n);
    gsl_vector *bcopy = gsl_vector_alloc(n);
    gsl_matrix *G2 = gsl_matrix_alloc(n,n);
    gsl_vector *x_prev = gsl_vector_alloc(n);
    gsl_vector *x_curr = gsl_vector_alloc(n);
    gsl_matrix *A = gsl_matrix_alloc(n,n);
    gsl_vector *b = gsl_vector_alloc(n);
    gsl_vector *e_curr = gsl_vector_alloc(n);
    gsl_vector *e_prev = gsl_vector_alloc(n);

   
   // CREATE TEXT OUTPUT + GNUPLOT FILES
    FILE *fp[tran->plot_count];
    char gpfile[tran->plot_count][256];

    for(int p=0;p<tran->plot_count;p++){
        char name[256];
        sprintf(name,"%s%d/tran_V%s.txt", outdir, use_iterative,
                nodes[tran->plot_nodes[p]].name);
        fp[p] = fopen(name,"w");

        sprintf(gpfile[p],"%s%d/tran_V%s.plt", outdir, use_iterative,
                nodes[tran->plot_nodes[p]].name);
        FILE *g = fopen(gpfile[p],"w");
        fprintf(g,
            "set title 'Transient Node %s'\n"
            "set xlabel 'Time (s)'\n"
            "set ylabel 'Voltage (V)'\n"
            "set grid\n"
            "set terminal png\n"
            "set output '%s%d/tran_V%s.png'\n"
            "plot '%s' using 1:2 with lines title 'V(%s)'\n",
            nodes[tran->plot_nodes[p]].name,
            outdir, use_iterative, nodes[tran->plot_nodes[p]].name,
            name, nodes[tran->plot_nodes[p]].name
        );
        fclose(g);
    }


    gsl_matrix_memcpy(Acopy, Gdc);
    gsl_vector_memcpy(bcopy, b_dc);
    solve_system_gsl(Acopy,bcopy,x_prev,use_cholesky);

    gsl_matrix_scale(C, scale);
    if (tran->method){
        gsl_matrix_memcpy(G2, Gdc);
        gsl_matrix_scale(G2, -1);
        gsl_matrix_add(G2, C);
        gsl_vector_memcpy(e_prev, b_dc);
    }
    
    gsl_matrix_add(Gdc, C);
    


    for (double t=0 ; t < (tran->final_time +tran->time_step); t = t + tran->time_step) {
        gsl_matrix_memcpy(A, Gdc);
        gsl_vector_memcpy(b,b_dc);
        if (tran->method) {
            // e_curr = e(t_k)
            gsl_vector_memcpy(e_curr, b_dc);
            stamp_transient_sources(e_curr, tran_sources, tran_src_number, m2, elements, node_count, t);

            // b = e_curr + e_prev
            gsl_vector_memcpy(b, e_curr);
            gsl_vector_add(b, e_prev);

            // b += G2 * x_prev
            gsl_blas_dgemv(CblasNoTrans, 1.0, G2, x_prev, 1.0, b);

            // e_prev = e_curr
            gsl_vector_memcpy(e_prev, e_curr);

            /*
            gsl_vector_memcpy(e_curr,b_dc);
            stamp_transient_sources(e_curr, tran_sources, tran_src_number, m2, elements, node_count, t);
            gsl_vector_add(e_prev, e_curr);
            gsl_vector_memcpy(b, e_prev);
            gsl_vector_memcpy(e_prev, e_curr);
            gsl_blas_dgemv(CblasNoTrans, 1, G2, x_prev, 1.0, b);*/
        } else {
            stamp_transient_sources(b, tran_sources, tran_src_number, m2, elements, node_count, t);
            gsl_blas_dgemv(CblasNoTrans, 1, C, x_prev, 1.0, b);
        }
        if (use_iterative)
            iterative_system_solver(A, b, x_curr, itol, 2*n, *use_cholesky);
        else if (use_custom)
            solve_system_custom(A, b, x_curr, use_cholesky);
        else
            solve_system_gsl(A, b, x_curr, use_cholesky);
        gsl_vector_memcpy(x_prev, x_curr);

        for(int p=0;p<tran->plot_count;p++){
            int id=tran->plot_nodes[p]-1;
            fprintf(fp[p],"%lf %lf\n",t,V(x_curr,id));
        }
    }

    // CLOSE AND GNUPLOT RENDER 
    for(int p=0;p<tran->plot_count;p++){
        fclose(fp[p]);
        char cmd[256];
        sprintf(cmd,"gnuplot %s", gpfile[p]);
        int ret=system(cmd);//PNG
        (void)ret;
    }


    gsl_matrix_free(Acopy);
    gsl_vector_free(bcopy);
    gsl_matrix_free(G2);
    gsl_vector_free(x_prev);
    gsl_vector_free(x_curr);
    gsl_matrix_free(A);
    gsl_vector_free(b);
    gsl_vector_free(e_curr);
    gsl_vector_free(e_prev);

    return 0;
}


/**
 * Sparse-matrix implementation of run_transient_analysis().
 *
 * Uses CSparse matrices and vectors while preserving identical
 * numerical formulation, time-stepping scheme and solver logic.
 *
 * See run_transient_analysis() for full algorithm description.
 */
int run_transient_analysis_sparse(
    cs *Gdc,              // ORIGINAL G matrix (DO NOT MODIFY)
    cs *C,                // ORIGINAL C matrix (DO NOT MODIFY)
    gsl_vector *b_dc,
    Element *m2, Element *elements,
    Node *nodes, int node_count,
    TRAN_Analysis *tran,
    TranSrc *tran_sources, int tran_src_number,
    int use_iterative,
    int *use_cholesky,
    double itol,
    char outdir[256]
){
    int n = b_dc->size;
    double alpha = (tran->method ? 2.0/tran->time_step
                                 : 1.0/tran->time_step);
    int status0 = 0, status1 = 0;
    
    //   Allocate vectors
    gsl_vector *x_prev = gsl_vector_calloc(n);
    gsl_vector *x_curr = gsl_vector_calloc(n);
    gsl_vector *b      = gsl_vector_calloc(n);
    gsl_vector *tmp    = gsl_vector_calloc(n);
    gsl_vector *e_prev = gsl_vector_calloc(n);
    gsl_vector *e_curr = gsl_vector_calloc(n);

   
    cs *C_scaled = cs_add(C, C, alpha, 0.0);   // αC
    cs_dupl(C_scaled);

    // G2 = αC - Gdc  (uses ORIGINAL Gdc) 
    cs *G2 = NULL;
    if (tran->method) {
        G2 = cs_add(C_scaled, Gdc, 1.0, -1.0);
        cs_dupl(G2);
        gsl_vector_memcpy(e_prev, b_dc);
        //stamp_transient_sources(e_prev, tran_sources, tran_src_number,m2, elements, node_count, 0);
    }

    // A = Gdc + αC  
    cs *A = cs_add(Gdc, C_scaled, 1.0, 1.0);
    cs_dupl(A);

    css *Sdc = NULL;
    csn *Ndc = NULL;
    css *S = NULL;
    csn *N = NULL;

    if(*use_cholesky) {
        status0 = sparse_cholesky_factorize(Gdc, &Sdc, &Ndc);
        status1 = sparse_cholesky_factorize(A, &S, &N);
        *use_cholesky = (!status0 && !status1);
        if(status0 == -1 || status1 == -1) {
            printf("Sparse Cholesky failed. Falling back to LU.\n");
            sparse_lu_factorize(Gdc, &Sdc, &Ndc);
            sparse_lu_factorize(A, &S, &N);
        }
    } else {
        sparse_lu_factorize(Gdc, &Sdc, &Ndc);
        sparse_lu_factorize(A, &S, &N);
    }

    //   DC operating point
    solve_system_sparse(Gdc, b_dc, x_prev, n, use_cholesky, Sdc, Ndc);


    //   Output files
    FILE *fp[tran->plot_count];
    char gpfile[tran->plot_count][256];

    for (int p = 0; p < tran->plot_count; p++) {
        char name[256];
        sprintf(name, "%s%d/sparse_tran_V%s.txt",
                outdir, use_iterative,
                nodes[tran->plot_nodes[p]].name);
        fp[p] = fopen(name, "w");

        sprintf(gpfile[p], "%s%d/sparse_tran_V%s.plt",
                outdir, use_iterative,
                nodes[tran->plot_nodes[p]].name);
        FILE *g = fopen(gpfile[p], "w");
        fprintf(g,
            "set title 'Transient Node %s'\n"
            "set xlabel 'Time (s)'\n"
            "set ylabel 'Voltage (V)'\n"
            "set grid\n"
            "set terminal png\n"
            "set output '%s%d/sparse_tran_V%s.png'\n"
            "plot '%s' using 1:2 with lines title 'V(%s)'\n",
            nodes[tran->plot_nodes[p]].name,
            outdir, use_iterative,
            nodes[tran->plot_nodes[p]].name,
            name,
            nodes[tran->plot_nodes[p]].name
        );
        fclose(g);
    }

    
    //   Time stepping
    for (double t = 0;
         t < tran->final_time + tran->time_step;
         t += tran->time_step)
    {
        gsl_vector_memcpy(b, b_dc);

        if (tran->method) {
           
            gsl_vector_memcpy(e_curr, b_dc);
            stamp_transient_sources(e_curr, tran_sources, tran_src_number, m2, elements, node_count, t);

            gsl_vector_add(e_prev, e_curr);
            gsl_vector_memcpy(b, e_prev);
            gsl_vector_memcpy(e_prev, e_curr);
           
            sparse_Ax(G2, x_prev->data, tmp->data);
            gsl_vector_add(b, tmp);
            
        }
        else {
            stamp_transient_sources(b, tran_sources, tran_src_number, m2, elements, node_count, t);
            sparse_Ax(C_scaled, x_prev->data, tmp->data);
            gsl_vector_add(b, tmp);
        }

        if (use_iterative) {
            gsl_vector_set_zero(x_curr);
            iterative_system_solver_sparse(A, b, x_curr, itol, n, *use_cholesky);
        }else
            solve_system_sparse(A, b, x_curr, n, use_cholesky, S, N);

        gsl_vector_memcpy(x_prev, x_curr);

        for (int p = 0; p < tran->plot_count; p++) {
            int id = tran->plot_nodes[p] - 1;
            fprintf(fp[p], "%e %e\n", t, V(x_curr, id));
        }
    }

   
    //   Cleanup
    for (int p = 0; p < tran->plot_count; p++) {
        fclose(fp[p]);
        char cmd[256];
        sprintf(cmd, "gnuplot %s", gpfile[p]);
        int ret=system(cmd);
        (void)ret;
    }

    cs_spfree(C_scaled);
    cs_spfree(A);
    if (G2) cs_spfree(G2);

    gsl_vector_free(x_prev);
    gsl_vector_free(x_curr);
    gsl_vector_free(b);
    gsl_vector_free(tmp);
    gsl_vector_free(e_prev);
    gsl_vector_free(e_curr);
    cs_sfree(Sdc);
    cs_nfree(Ndc);
    cs_sfree(S);
    cs_nfree(N);


    return 0;
}
