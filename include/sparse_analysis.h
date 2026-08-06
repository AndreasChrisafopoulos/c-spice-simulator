#ifndef SPARSE_ANALYSIS_H
#define SPARSE_ANALYSIS_H

#include "structs.h"
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include "csparse.h"

int run_dc_analysis_sparse(int *use_cholesky, int use_custom,
    Element *m2_elements, Node *nodes, int node_count, int m2_count,
    cs *A_csc, gsl_vector *b, int use_iterative,
    double itol, char output_dir[256]);

void run_all_sweeps_sparse(Element *elements, Element *m2_elements,
    cs *A_csc, gsl_vector *b_original, ElementID **elementIDs,
    int *use_cholesky, int node_count, DCSweep *dc_sweeps,
    int dc_sweep_count, int element_count, Node *nodes,
    int use_iterative, double itol, char output_dir[256]);

int run_dc_sweep_analysis_sparse(Element *elements, Element *m2_elements,
    double start, double end, double step, cs *A_csc,
    gsl_vector *b_original, int *use_cholesky, int node_count,
    const char *outfile, int b_index, int plot_count,
    int plot_nodes[MAX_SWEEP_PLOTS], int sweep_type,
    Node *nodes, int use_iterative, double itol, css *S, csn *N);

#endif
