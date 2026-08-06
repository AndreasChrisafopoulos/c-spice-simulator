#ifndef UTILS_H
#define UTILS_H

#include "structs.h"
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

unsigned int hash_string(const char *s, int table_size);
void str_to_lower(char *s);
int find_elementID(const char *name, int table_size, ElementID **table);

int create_output_directory_for_netlist(
    const char *netlist_filename,
    const char *base_dir,
    char *out_dir
);


// cleanup
void free_nodes(NodeHashtable **table, int table_size);
void free_elements(Element *elements, int element_count);
void free_elementIDs(ElementID **table, int table_size);
void cleanup(Element *elements, Element *m2_elements, int m2_count,
             int element_count, NodeHashtable **nodes,
             ElementID **elementIDs, int node_table_size);

#endif
