#include "analysis.h"
#include "solver_direct.h"
#include "solver_iterative.h"
#include "stamp.h"
#include "utils.h"
#include <stdio.h>


/**
 * Solves the DC operating point of the circuit.
 * Selects the appropriate solver based on user options and
 * writes node voltages and branch currents to an output file.
 */
int run_dc_analysis(
    int *use_cholesky, 
    int use_custom, 
    Element *m2_elements, 
    Node *nodes, 
    int node_count, 
    int m2_count, 
    gsl_matrix *A, 
    gsl_vector *b, 
    int use_iterative, 
    double itol, 
    char output_dir[256]
) 
{
    int size = A->size1;
    int status=0;

    // Allocate solution vector
    gsl_vector *x_out = gsl_vector_calloc(size);

    // Select solver based on user options
    if (use_iterative) 
        status = iterative_system_solver(A, b, x_out, itol, size, *use_cholesky);
    else if (use_custom) 
        status = solve_system_custom(A, b, x_out, use_cholesky); //status: -3 lu failed, table singular 
    else 
        status = solve_system_gsl(A, b, x_out, use_cholesky);

    // Write results if solution is available
    if (!status || use_iterative) {
        char path[512];
        sprintf(path, "%s%d/dc_op.txt", output_dir, use_iterative);
        FILE *fp = fopen(path, "w");
        if (!fp) {
            perror("dc_op.txt");
            return -2;
        }

        double *x = x_out->data;

        // 1) Voltages  
        for (int i = 1; i < node_count; i++) {
            double v = x[i - 1];
            str_to_lower(nodes[i].name);
            if (m2_count==0) {
                fprintf(fp, "%s = %.6e\n",  nodes[i].name, v);
            } else 
                fprintf(fp, "%s = %.6e\n", nodes[i].name, v);
        }

        // 2) Branch currents 
        int base = node_count - 1;
        for (int j = 0; j < m2_count; j++) {
            double I = x[base + j];
            TwoPortsSource *tp = (TwoPortsSource *)m2_elements[j].data;
            str_to_lower(tp->name);
            fprintf(fp, "%s#branch = %.6e\n", tp->name, I);
        }
        fclose(fp);
    }
    
    gsl_vector_free(x_out);

    return status;
}



/*
 * Executes all DC sweep analyses defined in the netlist.
 *
 * For each .DC directive, the corresponding source is located,
 * an output prefix is generated, and a DC sweep analysis is
 * performed for all requested plot nodes.
 *
 * Supports both voltage and current source sweeps and delegates
 * the actual sweep computation to run_dc_sweep_analysis().
 */
void run_all_sweeps(
    Element *elements, 
    Element *m2_elements,
    gsl_matrix *A_original, 
    gsl_vector *b_original, 
    ElementID **elementIDs,
    int use_custom, 
    int *use_cholesky,
    int node_count, 
    DCSweep *dc_sweeps, 
    int dc_sweep_count, 
    int element_count, 
    Node *nodes,
    int use_iterative, 
    double itol, 
    char output_dir[256]
)
{
    for (int i = 0; i < dc_sweep_count; i++) {
        int b_index = node_count + find_elementID(dc_sweeps[i].source_name, element_count, elementIDs );

        if (dc_sweeps[i].plot_count == 0) {
            printf("  No .PLOT entries\n");
        } else {
            char outfile[64];
            sprintf(outfile, "%s%d/dc_sweep_%s", output_dir, use_iterative, dc_sweeps[i].source_name);

            char first = toupper(dc_sweeps[i].source_name[0]);
            if(first =='V')
                run_dc_sweep_analysis(elements, m2_elements, dc_sweeps[i].start, dc_sweeps[i].end, dc_sweeps[i].step, A_original, b_original, use_custom, use_cholesky, node_count, outfile, b_index, dc_sweeps[i].plot_count, dc_sweeps[i].plot_nodes, 1, nodes, use_iterative, itol);
            else 
                run_dc_sweep_analysis(elements, m2_elements, dc_sweeps[i].start, dc_sweeps[i].end, dc_sweeps[i].step, A_original, b_original, use_custom, use_cholesky, node_count, outfile, b_index, dc_sweeps[i].plot_count, dc_sweeps[i].plot_nodes, 0, nodes, use_iterative, itol);
        }
    }

}


/*
 * Performs a DC sweep analysis on the circuit.
 * 
 * The function repeatedly solves the DC operating point while sweeping
 * a source value from start to end with the given step. For each sweep
 * point, the system matrix and RHS vector are restored, updated with
 * the new source value, and solved using the selected solver.
 *
 * Node voltages specified by .PLOT directives are written to output files,
 * and corresponding gnuplot scripts are generated to visualize the results.
 */
int run_dc_sweep_analysis(
    Element *elements, Element *m2_elements,
    double start, double end, double step,
    gsl_matrix *A_original,
    gsl_vector *b_original,
    int use_custom,
    int *use_cholesky,
    int node_count,
    const char *outfile,
    int b_index,
    int plot_count,
    int plot_nodes[MAX_SWEEP_PLOTS],
    int sweep_type, Node *nodes,
    int use_iterative,
    double itol
)
{
    size_t n = A_original->size1;

    gsl_matrix *A_work = gsl_matrix_calloc(n, n);
    gsl_vector *b_work = gsl_vector_calloc(n);
    gsl_vector *x_out  = gsl_vector_calloc(n);


    FILE *gps[MAX_SWEEP_PLOTS];
    FILE *fps[MAX_SWEEP_PLOTS];
    for (int p = 0; p < plot_count; p++) {
        char fname[256];
        sprintf(fname, "%s_v%s.txt", outfile, nodes[plot_nodes[p]].name);
        fps[p] = fopen(fname, "w");

        if (!fps[p]) {
            perror(fname);
            return -1;
        }

        char gscript[256];
        sprintf(gscript, "%s_v%s.plt", outfile, nodes[plot_nodes[p]].name);
        gps[p] = fopen(gscript, "w");

        fprintf(gps[p],
            "set title 'DC Sweep: %s node %s'\n"
            "set xlabel 'Source value'\n"
            "set ylabel 'Voltage (V)'\n"
            "set terminal png\n"
            "set output '%s_v(%s).png'\n"
            "plot '%s' using 1:2 with lines title '%s'\n",
            outfile, nodes[plot_nodes[p]].name,
            outfile, nodes[plot_nodes[p]].name,
            fname,
            nodes[plot_nodes[p]].name
        );
        fclose(gps[p]);
    }


    for (double val = start; val <= end+step ; val += step)
    {
        // restore A and B
        gsl_matrix_memcpy(A_work, A_original);
        gsl_vector_memcpy(b_work, b_original);

        // apply new sweep stamp
        if(sweep_type){
            TwoPortsSource *tp = (TwoPortsSource *)m2_elements[b_index-node_count].data;
           // printf("\nname= %s", tp->name);

            sweep_stamp(
                b_work,
                tp->ports[0],
                tp->ports[1],
                (val),
                1,
                (b_index-1)   // voltage sources
            );
        //  printf("b_[] = %lf\n", gsl_vector_get(b_work, b_index-1));
        } else if (sweep_type==0){
            TwoPortsSource *tp = (TwoPortsSource *)elements[b_index-node_count].data;
            sweep_stamp(
                b_work,
                tp->ports[0],
                tp->ports[1],
                val-(tp->value),
                0,
                (b_index-1)   
            );
        }
        // solve
        int status=0, maxIter = A_original->size1;
        if (use_iterative)
            iterative_system_solver(A_work, b_work, x_out, itol, maxIter, *use_cholesky );
        else if (use_custom)
            status = solve_system_custom(A_work, b_work, x_out, use_cholesky);
        else
            status = solve_system_gsl(A_work, b_work, x_out, use_cholesky);

        if (status) {
            fprintf(fps[0], "%e ERROR\n", val);
            continue;
        }

        // Write node voltage

        if (plot_count == 0) {
            printf("  No .PLOT entries\n");
        } else {
            for (int p = 0; p < plot_count; p++) {
                int node_id = plot_nodes[p];
                double v_node = gsl_vector_get(x_out, node_id-1 );
                fprintf(fps[p], "%e  %e\n", val, v_node);
            }
        }
    }
    for (int p = 0; p < plot_count; p++) {
        fclose(fps[p]);

        char cmd[256];
        sprintf(cmd, "gnuplot %s_v%s.plt",
            outfile,
            nodes[plot_nodes[p]].name
        );

        int ret = system(cmd);
        (void)ret;

    }

    // fclose(fp);
    gsl_matrix_free(A_work);
    gsl_vector_free(b_work);
    gsl_vector_free(x_out);

    return 0;
}

