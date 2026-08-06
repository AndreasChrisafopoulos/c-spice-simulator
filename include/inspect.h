#ifndef INSPECT_H
#define INSPECT_H

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include "structs.h"  
// printing
void print_first_pass_results(
    int element_count,
    int node_count,
    int m2_count,
    int use_cholesky,
    int use_custom,
    int use_iterative,
    double itol,
    int use_sparse,
    int do_trans,
    TRAN_Analysis tran,
    int tran_src_number  
);
void print_dc_sweeps(DCSweep *sweeps, int sweep_count, Node *nodes);
void print_transient_spec(TwoPortsSource *s);
void print_transient_info(TRAN_Analysis *tran, Node *nodes);
void print_elements(Element *elements, Element *m2_elements, int m2_size, int size);
void print_all_nodes(NodeHashtable **table, int table_size);
void print_elementIDs(ElementID **table, int table_size);
void print_DCTables(gsl_matrix *A, gsl_vector *b, int DCtable_size);
void print_solution(const gsl_vector *x);
void print_nodes_by_id(Node *nodes, int node_count);
void print_tran_sources(TranSrc *ts, int count);

#endif