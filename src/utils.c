#include "utils.h"
#include <sys/stat.h>
#include <sys/types.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <strings.h>

/**
 * Finds the index of an element by name using the ElementID hash table.
 *
 * Performs a case-insensitive lookup in the hash bucket corresponding
 * to the given name and returns the stored index of the element.
 */
int find_elementID(const char *name, int table_size, ElementID **table)
{
    unsigned int idx = hash_string(name, table_size);
    ElementID *curr = table[idx];

    while (curr) {
        if (strcasecmp(curr->name, name) == 0) {
            return curr->index;
        }
        curr = curr->next;
    }

    return -1;   // not found
}

unsigned int hash_string(const char *s, int table_size) {
    unsigned long hash = 5381;
    int c;
    while ((c = (unsigned char)*s++))
        hash = ((hash << 5) + hash) + c;
    return (unsigned int)(hash % table_size);
}

void str_to_lower(char *s) {
    for (; *s; ++s)
        *s = tolower((unsigned char)*s);
}

/**
 *  Creates two output directories for the analysis results.
 *
 * Given a netlist filename (e.g., "circuit.cir"), this function removes the
 * extension, appends "_outputfiles", and creates two folders:
 *
 *      <base>0   and   <base>1
 *
 * This is done so that later the program can store:
 *  - non-iterative analysis results in <base>0
 *  - iterative analysis results in <base>1
 *
 * The function stores the base folder name (without 0/1) in out_dir
 * and works on both Linux/Unix and Windows.
 *
 * return 0 on success, -1 on error.
 */
int create_output_directory_for_netlist(
    const char *netlist_filename,
    const char *base_dir,
    char *out_dir
)
{
    char base_name[256];
    strcpy(base_name, netlist_filename);

    char *dot = strrchr(base_name, '.');
    if (dot) *dot = '\0';

    // OUT/<name>_outputfiles
    snprintf(out_dir, 512, "%s/%s_outputfiles", base_dir, base_name);

    char out_dir0[512];
    char out_dir1[512];

    snprintf(out_dir0, sizeof(out_dir0), "%s0", out_dir);
    snprintf(out_dir1, sizeof(out_dir1), "%s1", out_dir);

    struct stat st = {0};

    if (stat(out_dir0, &st) == -1) {
        #ifdef _WIN32
                if (_mkdir(out_dir0) != 0) {
        #else
                if (mkdir(out_dir0, 0755) != 0) {
        #endif
                    perror(out_dir0);
                    return -1;
                }
    }

    memset(&st, 0, sizeof(st));

    if (stat(out_dir1, &st) == -1) {
        #ifdef _WIN32
                if (_mkdir(out_dir1) != 0) {
        #else
                if (mkdir(out_dir1, 0755) != 0) {
        #endif
                    perror(out_dir1);
                    return -1;
                }
    }

    return 0;
}


// Memory cleanup
void free_nodes(NodeHashtable **table, int table_size) {
    for (int i = 0; i < table_size; i++) {
        NodeHashtable *curr = table[i];
        while (curr) {
            NodeHashtable *next = curr->next;
            free(curr);
            curr = next;
        }
    }
    free(table);
}

void free_elements(Element *elements, int element_count) {
    for (int i = 0; i < element_count; i++) {
        free(elements[i].data);
    }
    free(elements);
}

void free_elementIDs(ElementID **table, int table_size) {
    for (int i = 0; i < table_size; i++) {
        ElementID *curr = table[i];
        while (curr) {
            ElementID *next = curr->next;
            free(curr);
            curr = next;
        }
    }
    free(table);
}

void free_DCtables(double **A, double *b, int DCtable_size) {
    if (A != NULL) {
        for (int i = 0; i < DCtable_size; i++) {
            free(A[i]);   
        }
        free(A);  
    }

    if (b != NULL) {
        free(b);  
    }
}


void cleanup(Element *elements, Element *m2_elements, int m2_count, int element_count, NodeHashtable **nodes, ElementID **elementIDs, int node_table_size) {
    free_elements(elements, (element_count-m2_count));
    free_elements(m2_elements, m2_count);
    free_nodes(nodes, node_table_size);
    free_elementIDs(elementIDs, element_count);
}
