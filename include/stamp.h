#ifndef STAMP_H
#define STAMP_H

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include "csparse.h"
#include <stdio.h>
#include <stdlib.h>
#include "structs.h"
void stamp(char type, int n1, int n2, gsl_matrix *A, gsl_vector *b, int *m2_index, double value, int node_count, int use_sparse, cs *A_triplet);

int sweep_stamp(gsl_vector *b, int n1, int n2, double value, int type, int index);

void stamp_C(char type, int n1, int n2,
             gsl_matrix *C,
             int *m2_index,
             double value,
             int node_count,
             int use_sparse,
             cs *C_triplet);

#endif
