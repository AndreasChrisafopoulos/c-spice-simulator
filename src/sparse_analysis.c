#include "sparse_analysis.h"
#include "sparse_solver.h"
#include "stamp.h"
#include "utils.h"
#include <stdio.h>

/* DC operating point analysis using sparse matrices (direct or iterative solver). */
int run_dc_analysis_sparse(
    int *use_cholesky,
    int use_custom,
    Element *m2_elements,
    Node *nodes,
    int node_count,
    int m2_count,
    cs *A_csc,          // sparse matrix in triplet form
    gsl_vector *b,
    int use_iterative,
    double itol,
    char output_dir[256]
){
    int size = A_csc->n;
    int status = 0;
    css *S = NULL;
    csn *N = NULL;
    int method=0;

    gsl_vector *x_out = gsl_vector_calloc(size);
    
    // SOLVER SELECTION 

    if (use_iterative) {
        // iterative sparse CG / BiCG 
        status = iterative_system_solver_sparse(
            A_csc, b, x_out, itol, size, *use_cholesky
        );
    }
    else {
        // direct sparse LU / Cholesky 
        if(*use_cholesky){
            method = sparse_cholesky_factorize(A_csc, &S, &N);
            *use_cholesky = !method;
            if(method == -1) {
                printf("Sparse Cholesky failed. Falling back to LU.\n");
                sparse_lu_factorize(A_csc, &S, &N);
            }
        }else
            sparse_lu_factorize(A_csc, &S, &N);

        status = solve_system_sparse(
            A_csc, b, x_out, size, use_cholesky, S, N
        );
    }

    //  OUTPUT 
    if (!status || use_iterative) {
        char path[512];
        sprintf(path, "%s%d/dc_op.txt", output_dir, use_iterative);

        FILE *fp = fopen(path, "w");
        if (!fp) {
            perror("dc_op.txt");
            gsl_vector_free(x_out);
            return -2;
        }

        double *x = x_out->data;

        // 1) Node voltages 
        for (int i = 1; i < node_count; i++) {
            double v = x[i - 1];
            str_to_lower(nodes[i].name);
            fprintf(fp, "%s = %.5e\n", nodes[i].name, v);
        }

        // 2) Branch currents 
        int base = node_count - 1;
        for (int j = 0; j < m2_count; j++) {
            double I = x[base + j];
            TwoPortsSource *tp = (TwoPortsSource *)m2_elements[j].data;
            str_to_lower(tp->name);
            fprintf(fp, "%s#branch = %.5e\n", tp->name, I);
        }

        fclose(fp);
    }


    cs_sfree(S);
    cs_nfree(N);

    gsl_vector_free(x_out);
    return status;
}



/* Executes all DC sweep analyses (sparse MNA, direct or iterative). */
void run_all_sweeps_sparse(
    Element *elements,
    Element *m2_elements,
    cs *A_csc,
    gsl_vector *b_original,
    ElementID **elementIDs,
    //int use_custom,
    int *use_cholesky,
    int node_count,
    DCSweep *dc_sweeps,
    int dc_sweep_count,
    int element_count,
    Node *nodes,
    int use_iterative,
    double itol,
    char output_dir[256]
){
    css *S = NULL;
    csn *N = NULL;
    int status=0;

    if(*use_cholesky){
        status = sparse_cholesky_factorize(A_csc, &S, &N);
        *use_cholesky = !status;
        if(status == -1) {
            printf("Sparse Cholesky failed. Falling back to LU.\n");
            sparse_lu_factorize(A_csc, &S, &N);
        }
    }
    else
        sparse_lu_factorize(A_csc, &S, &N);

    
    for (int i = 0; i < dc_sweep_count; i++) {

        int b_index = node_count +
            find_elementID(dc_sweeps[i].source_name,
                           element_count, elementIDs);

        char outfile[128];
        sprintf(outfile, "%s%d/dc_sweep_%s",
                output_dir, use_iterative,
                dc_sweeps[i].source_name);

        char first = toupper(dc_sweeps[i].source_name[0]);

        run_dc_sweep_analysis_sparse(
            elements,
            m2_elements,
            dc_sweeps[i].start,
            dc_sweeps[i].end,
            dc_sweeps[i].step,
            A_csc,
            b_original,
            //use_custom,
            use_cholesky,
            node_count,
            outfile,
            b_index,
            dc_sweeps[i].plot_count,
            dc_sweeps[i].plot_nodes,
            (first == 'V'),
            nodes,
            use_iterative,
            itol,
            S,
            N
        );
    }

    
    cs_sfree(S);
    cs_nfree(N);

}

/* Performs a DC sweep using sparse solvers and writes node voltages for plotting. */
int run_dc_sweep_analysis_sparse(
    Element *elements,
    Element *m2_elements,
    double start,
    double end,
    double step,
    cs *A_csc,
    gsl_vector *b_original,
    int *use_cholesky,
    int node_count,
    const char *outfile,
    int b_index,
    int plot_count,
    int plot_nodes[MAX_SWEEP_PLOTS],
    int sweep_type,
    Node *nodes,
    int use_iterative,
    double itol, css *S, csn *N
){
    int n = A_csc->n;
    int iters =0;
    gsl_vector *b_work = gsl_vector_calloc(n);
    gsl_vector *x_out  = gsl_vector_calloc(n);

    FILE *fps[MAX_SWEEP_PLOTS];
    for (int p = 0; p < plot_count; p++) {
        char fname[256];
        sprintf(fname, "%s_v%s.txt", outfile, nodes[plot_nodes[p]].name);
        fps[p] = fopen(fname, "w");
    }

    // compress  
    //cs *A_csc = cs_compress(A_triplet);
    //cs_dupl(A_csc);

    const double eps = 1e-12;

    for (double val = start; val <= end + eps; val += step)
    {
        gsl_vector_memcpy(b_work, b_original);

        // apply sweep stamp 
        if (sweep_type) {
            TwoPortsSource *tp =
                (TwoPortsSource *)m2_elements[b_index-node_count].data;

            sweep_stamp(
                b_work,
                tp->ports[0],
                tp->ports[1],
                val,
                1,
                b_index-1
            );
        } else {
            TwoPortsSource *tp =
                (TwoPortsSource *)elements[b_index-node_count].data;

            sweep_stamp(
                b_work,
                tp->ports[0],
                tp->ports[1],
                val - tp->value,
                0,
                b_index-1
            );
        }

        int status=0;
        if (use_iterative) {
                iters = iterative_system_solver_sparse(
                A_csc, b_work, x_out, itol, n, *use_cholesky
            );
        } else {
            status = solve_system_sparse(
                A_csc, b_work, x_out, n, use_cholesky, S, N
            );
        }

        if (status) continue;

        for (int p = 0; p < plot_count; p++) {
            int node_id = plot_nodes[p];
            double v = gsl_vector_get(x_out, node_id - 1);
            fprintf(fps[p], "%e  %e\n", val, v);
        }
    }

    for (int p = 0; p < plot_count; p++)
        fclose(fps[p]);

    gsl_vector_free(b_work);
    gsl_vector_free(x_out);
    return iters;
}

