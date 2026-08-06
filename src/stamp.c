#include "stamp.h"

#define V(v,i) gsl_vector_get(v,i)
#define M(m,i,j) gsl_matrix_get(m,i,j)

/* MNA stamping for R, I, V, L elements (supports dense & sparse matrices) */
void stamp(char type, int n1, int n2,
           gsl_matrix *A,
           gsl_vector *b,
           int *m2_index,
           double value,
           int node_count,
           int use_sparse,
           cs *A_triplet)
{
    if (type == 'R') {

        double g = 1.0 / value;

        if (n1 != 0) {
            if (use_sparse)
                cs_entry(A_triplet, n1-1, n1-1, g);
            else
                gsl_matrix_set(A, n1-1, n1-1,
                    gsl_matrix_get(A, n1-1, n1-1) + g);
        }

        if (n2 != 0) {
            if (use_sparse)
                cs_entry(A_triplet, n2-1, n2-1, g);
            else
                gsl_matrix_set(A, n2-1, n2-1,
                    gsl_matrix_get(A, n2-1, n2-1) + g);
        }

        if (n1 != 0 && n2 != 0) {
            if (use_sparse) {
                cs_entry(A_triplet, n1-1, n2-1, -g);
                cs_entry(A_triplet, n2-1, n1-1, -g);
            } else {
                gsl_matrix_set(A, n1-1, n2-1,
                    gsl_matrix_get(A, n1-1, n2-1) - g);
                gsl_matrix_set(A, n2-1, n1-1,
                    gsl_matrix_get(A, n2-1, n1-1) - g);
            }
        }

    } else if (type == 'I') {

        if (n1 != 0)
            gsl_vector_set(b, n1-1,
                gsl_vector_get(b, n1-1) - value);
        if (n2 != 0)
            gsl_vector_set(b, n2-1,
                gsl_vector_get(b, n2-1) + value);

    } else if (type == 'V' || type == 'L') {

        int row = (node_count - 1) + *m2_index;

        if (n1 != 0) {
            if (use_sparse) {
                cs_entry(A_triplet, n1-1, row, 1.0);
                cs_entry(A_triplet, row, n1-1, 1.0);
            } else {
                gsl_matrix_set(A, n1-1, row,
                    gsl_matrix_get(A, n1-1, row) + 1.0);
                gsl_matrix_set(A, row, n1-1,
                    gsl_matrix_get(A, row, n1-1) + 1.0);
            }
        }

        if (n2 != 0) {
            if (use_sparse) {
                cs_entry(A_triplet, n2-1, row, -1.0);
                cs_entry(A_triplet, row, n2-1, -1.0);
            } else {
                gsl_matrix_set(A, n2-1, row,
                    gsl_matrix_get(A, n2-1, row) - 1.0);
                gsl_matrix_set(A, row, n2-1,
                    gsl_matrix_get(A, row, n2-1) - 1.0);
            }
        }

        if (type == 'V')
            gsl_vector_set(b, row, value);

        *m2_index = *m2_index + 1;
    }
}

/* Applies DC sweep contribution to RHS vector b (current or voltage source) */
int sweep_stamp(gsl_vector *b, int n1, int n2, double value, int type, int index) {

    if (type==0){
        if (n1 != 0)
            gsl_vector_set(b, n1-1,
                gsl_vector_get(b, n1-1) - value);
        if (n2 != 0)
            gsl_vector_set(b, n2-1,
                gsl_vector_get(b, n2-1) + value);
    }else if(type){
        gsl_vector_set(b, index, value);
    }
return 0;
}


/* C-matrix stamping for transient analysis (C and L elements) */
void stamp_C(char type, int n1, int n2,
             gsl_matrix *C,
             int *m2_index,
             double value,
             int node_count,
             int use_sparse,
             cs *C_triplet)
{
    //  CAPACITOR 
    if (type == 'C') {
        double c = value;

        if (n1 != 0) {
            if (use_sparse)
                cs_entry(C_triplet, n1-1, n1-1, c);
            else
                gsl_matrix_set(C, n1-1, n1-1,
                    gsl_matrix_get(C, n1-1, n1-1) + c);
        }

        if (n2 != 0) {
            if (use_sparse)
                cs_entry(C_triplet, n2-1, n2-1, c);
            else
                gsl_matrix_set(C, n2-1, n2-1,
                    gsl_matrix_get(C, n2-1, n2-1) + c);
        }

        if (n1 != 0 && n2 != 0) {
            if (use_sparse) {
                cs_entry(C_triplet, n1-1, n2-1, -c);
                cs_entry(C_triplet, n2-1, n1-1, -c);
            } else {
                gsl_matrix_set(C, n1-1, n2-1,
                    gsl_matrix_get(C, n1-1, n2-1) - c);
                gsl_matrix_set(C, n2-1, n1-1,
                    gsl_matrix_get(C, n2-1, n1-1) - c);
            }
        }
    }

    //  INDUCTOR (dynamic part) 
    else if (type == 'L') {
        int row = (node_count - 1) + *m2_index-1;
        double L = value;
        if (use_sparse)
                cs_entry(C_triplet, row, row, -L);
        else {
        //  C stamp (dynamic part) 
        gsl_matrix_set(C, row, row,
            gsl_matrix_get(C, row, row) - L);
        }
        //(*m2_index)++;
    }

}

