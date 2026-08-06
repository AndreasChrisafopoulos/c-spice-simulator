#define _POSIX_C_SOURCE 199309L
#include "functions.h"

#include <sys/stat.h>
#include <sys/types.h>


#include <time.h>



double get_time_sec() {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}


void print_cs_csc(const cs *A, const char *name)
{
    printf("\n=== %s (CSC format) ===\n", name);
    printf("m=%d n=%d nnz=%d\n", A->m, A->n, A->p[A->n]);

    for (int j = 0; j < A->n; j++) {
        for (int p = A->p[j]; p < A->p[j+1]; p++) {
            printf("(%d,%d) = %.6e\n", A->i[p], j, A->x[p]);
        }
    }
}



int main(int argc, char *argv[]) {
    if (argc < 3) {
        printf("Usage: %s <part_number> <filename>\n", argv[0]);
        return 1;
    }

    int part = atoi(argv[1]);     // 3
    char *fname = argv[2];        //  "part3_LU_10K_nodes.cir"

    char prefix[64];
    snprintf(prefix, sizeof(prefix), "Part%d_Netlists/", part);

    char file[256];
    snprintf(file, sizeof(file), "%s%s", prefix, fname);

    printf("Final file path: %s\n", file);


    mkdir("OUT", 0755);

    char output_dir[512];
    create_output_directory_for_netlist(fname, "OUT", output_dir);




    NodeHashtable **nodes_hashtable;              //hashtable of all nodes
    Node *nodes;
    Element *elements;         // table of all elemets except from V and L
    Element *m2_elements;   // V and L
    ElementID **elementIDs;    //hash table, stores all element names and their index in elements table 
    int nodeTable_size = 0, node_ids = 0;
    int element_count = 0;
    int m2_count=0, node_count=0;
    int DCtable_size=0;
    int m2_index=0;
    int use_cholesky=0, use_custom=0;
    DCSweep dc_sweeps[MAX_SWEEPS];
    int dc_sweep_count=0;
    int use_iterative=0;
    double itol=0.001;
    int use_sparse=0;
    int do_trans=0;
    TRAN_Analysis tran;
    int tran_src_number;


    first_pass(file, &element_count, &node_count, &m2_count, &use_cholesky, &use_custom, &use_iterative, &itol, &use_sparse, &do_trans, &tran, &tran_src_number);

 
    print_first_pass_results(element_count, node_count, m2_count, use_cholesky, use_custom, use_iterative, itol, use_sparse, do_trans, tran, tran_src_number);
    elements = (Element *)malloc(sizeof(Element) * (element_count - m2_count));
    m2_elements = (Element *)malloc(sizeof(Element) * m2_count);
    nodes_hashtable = (NodeHashtable **)calloc(2 * node_count, sizeof(NodeHashtable *));
    elementIDs = (ElementID **)calloc(element_count, sizeof(ElementID*));
    nodes = (Node *)calloc(node_count, sizeof(Node));
    nodes[0].id=0;
    strcpy(nodes[0].name,"ground");
    nodeTable_size = 2 * node_count;


    TranSrc *tran_sources = malloc(tran_src_number * sizeof(TranSrc));

    DCtable_size = m2_count + node_count - 1;
    gsl_vector *b = gsl_vector_calloc(DCtable_size);
    gsl_matrix *A = NULL;
    cs *A_triplet = NULL;

    gsl_matrix *C = NULL;
    cs *C_triplet = NULL;

    gsl_matrix *G = NULL;
    cs *G_triplet = NULL;

    if (use_sparse) {
        int nzmax = 10 * DCtable_size;   
        A_triplet = cs_spalloc(DCtable_size, DCtable_size, nzmax, 1, 1);
        C_triplet = cs_spalloc(DCtable_size, DCtable_size, nzmax, 1, 1);
        G_triplet = cs_spalloc(DCtable_size, DCtable_size, nzmax, 1, 1);
    } else {
        A = gsl_matrix_calloc(DCtable_size, DCtable_size);
        C = gsl_matrix_calloc(DCtable_size, DCtable_size);
        G = gsl_matrix_calloc(DCtable_size, DCtable_size);
        
        
    }


    parse_netlist(file, elements, m2_elements, nodes_hashtable, nodeTable_size, element_count, &node_ids, elementIDs, A, b, &m2_index, nodes, dc_sweeps, &dc_sweep_count, use_sparse, A_triplet, C, C_triplet, do_trans, &tran, tran_sources);
    double t0 = get_time_sec();
    
    print_transient_info(&tran, nodes);

    if(use_sparse){
        gsl_vector *b_cpy = gsl_vector_calloc(DCtable_size);
        gsl_vector_memcpy(b_cpy, b);
        int maxIter=DCtable_size;
        cs *A_csc = cs_compress(A_triplet);
        cs_dupl(A_csc);
      
        run_dc_analysis_sparse(&use_cholesky, use_custom, m2_elements, nodes, node_count, m2_count, A_csc, b, use_iterative, itol, output_dir);
        run_all_sweeps_sparse(elements, m2_elements, A_csc, b, elementIDs, &use_cholesky, node_count, dc_sweeps, dc_sweep_count, element_count, nodes, use_iterative, itol, output_dir);
          if (do_trans) {
            cs *C_csc = cs_compress(C_triplet);
            cs_dupl(C_csc);
            run_transient_analysis_sparse(A_csc, C_csc, b_cpy, m2_elements, elements, nodes, node_count, &tran, tran_sources, tran_src_number, use_iterative, &use_cholesky, itol, output_dir);
            cs_spfree(C_csc);
        }
         cs_spfree(A_csc);
         gsl_vector_free(b_cpy);
    } else{
        gsl_matrix *A_cpy = gsl_matrix_calloc(DCtable_size, DCtable_size);
        gsl_vector *b_cpy = gsl_vector_calloc(DCtable_size);
        gsl_matrix_memcpy(A_cpy, A);
        gsl_vector_memcpy(b_cpy, b);
        if (do_trans) {
            //gsl_matrix *A_cpy2 = gsl_matrix_calloc(DCtable_size, DCtable_size);
            gsl_vector *b_cpy2 = gsl_vector_calloc(DCtable_size);
            gsl_matrix_memcpy(G, A);
            gsl_vector_memcpy(b_cpy2, b);
            run_transient_analysis(G, C, b_cpy2, m2_elements, elements, nodes, node_count, &tran, tran_sources, tran_src_number, use_iterative, &use_cholesky, use_custom, itol, output_dir);
            gsl_vector_free(b_cpy2);
        }
        run_dc_analysis(&use_cholesky, use_custom, m2_elements, nodes, node_count, m2_count, A, b, use_iterative, itol, output_dir);
        run_all_sweeps(elements, m2_elements, A_cpy, b_cpy, elementIDs, use_custom, &use_cholesky, node_count, dc_sweeps, dc_sweep_count, element_count, nodes, use_iterative, itol, output_dir);   
        gsl_vector_free(b_cpy);
        gsl_matrix_free(A_cpy);    
    }
    double t1 = get_time_sec();
    printf("Time: %.6f sec\n", t1 - t0);


    if (use_sparse) {
        cs_spfree(A_triplet);
        cs_spfree(C_triplet);
        cs_spfree(G_triplet);
    } else {
        gsl_matrix_free(A);
        gsl_matrix_free(C);
        gsl_matrix_free(G);
    }

    gsl_vector_free(b);
    
    free(tran_sources);

    cleanup(elements, m2_elements, m2_count, element_count, nodes_hashtable, elementIDs, nodeTable_size);
    free(nodes);

    return 0;
}

/* ================= DEBUG HELPERS ================= */

/* Print DC matrices */
// print_DCTables(A, b, DCtable_size);

/* Print nodes */
// print_nodes_by_id(nodes, node_count);

/* Print elements */
// print_elements(elements, m2_elements, m2_count, element_count - m2_count);

/* Print element IDs */
// print_elementIDs(elementIDs, element_count);

/* Print transient info */
// print_transient_info(&tran, nodes);

/* Print sparse matrix (CSC) */
/*
cs *C_csc = cs_compress(C_triplet);
cs_dupl(C_csc);
print_cs_csc(C_csc, "C sparse");
cs_spfree(C_csc);
*/

/* Solve & print DC solution manually */
/*
gsl_vector *x_out = gsl_vector_calloc(DCtable_size);
solve_system_gsl(A, b, x_out, &use_cholesky);
print_solution(x_out);
gsl_vector_free(x_out);
*/

