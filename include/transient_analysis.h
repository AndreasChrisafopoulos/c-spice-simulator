#ifndef TRANSIENT_ANALYSIS_H
#define TRANSIENT_ANALYSIS_H

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include "structs.h"
#include "csparse.h"


int run_transient_analysis(
    gsl_matrix *Gdc, gsl_matrix *C, 
    gsl_vector *b_dc,
    Element *m2, Element *elements,
    Node *nodes, int node_count,
    TRAN_Analysis *tran, TranSrc *tran_sources, int tran_src_number,
    int use_iterative, int *use_cholesky, int use_custom, double itol,
    char outdir[256]
);
void stamp_transient_sources(
    gsl_vector *b,
    TranSrc *tran_sources, int tran_src_number,
    Element *m2, Element *elements,
    int node_count,
    double t
);

int run_transient_analysis_sparse(
    cs *Gdc,              
    cs *C,                
    gsl_vector *b_dc,
    Element *m2, Element *elements,
    Node *nodes, int node_count,
    TRAN_Analysis *tran,
    TranSrc *tran_sources, int tran_src_number,
    int use_iterative,
    int *use_cholesky,
    double itol,
    char outdir[256]
);

#endif
