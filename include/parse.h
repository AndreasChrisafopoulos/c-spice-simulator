#ifndef PARSE_H
#define PARSE_H

#include "structs.h"        
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include "csparse.h"        


void first_pass(const char *filename, int *element_count, int *node_count,
                int *VsrcAndIndNumber, int *use_cholesky, int *use_ccustom,
                int *use_iterative, double *itol, int *use_sparse,   int *has_trans,
                TRAN_Analysis *tran, int *tran_src_number  );

elemenT getType(char p);

int insert_node(NodeHashtable **table, const char *name, int table_size, int *node_counter);

void insert_elementID(int index, const char *name, int table_size, ElementID **table);
int parse_transient_spec(const char *line, TwoPortsSource *src);
void insert_element(Element *elements, Element *m2_elements, int *m2_index, int *index,
                    const char *line, NodeHashtable **nodes_hashtable,
                    int nodeTable_size, int elemntTable_size, int *node_counter,
                    ElementID **elementIDs, gsl_matrix *A, gsl_vector *b,
                    int *m2_index_forStamp, Node *nodes,
                    int use_sparse, cs *A_triplet, gsl_matrix *C,
                    cs *C_triplet, int do_trans, int *m2_index_forStampC, TranSrc *tran_sources, int *tran_src_cntr);

void parse_netlist(const char *filename, Element *elements, Element *m2_elements,
                   NodeHashtable **nodes_hashtable, int nodeTable_count,
                   int elementTable_size, int *node_counter,
                   ElementID **elementIDs, gsl_matrix *A, gsl_vector *b,
                   int *m2_index_forStamp, Node *nodes, DCSweep *dc_sweeps,
                   int *dc_sweep_count, int use_sparse, cs *A_triplet, gsl_matrix *C,
                   cs *C_triplet, int do_trans, TRAN_Analysis *tran, TranSrc *tran_sources);

#endif 
