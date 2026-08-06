#include "structs.h"
#include "inspect.h"
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>


void print_dc_sweeps(DCSweep *sweeps, int sweep_count, Node *nodes)
{
    printf("\n===== DC SWEEP PARSE RESULT =====\n");

    for (int i = 0; i < sweep_count; i++) {
        printf("\nSWEEP #%d\n", i);
        printf("  Source       : %s\n", sweeps[i].source_name);
        printf("  Start        : %g\n", sweeps[i].start);
        printf("  End          : %g\n", sweeps[i].end);
        printf("  Step         : %g\n", sweeps[i].step);

        if (sweeps[i].plot_count == 0) {
            printf("  No .PLOT entries\n");
        } else {
            printf("  Will plot    :\n");

            for (int p = 0; p < sweeps[i].plot_count; p++) {
                int node_id = sweeps[i].plot_nodes[p];

                
                printf("       V(%s)  [node_id=%d]\n", 
                   nodes[node_id].name, node_id);
            }
        }
    }

    printf("\n=================================\n");
}

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
) {
    printf("\n===== FIRST PASS RESULTS =====\n");

    printf("Total elements:      %d\n", element_count);
    printf("Total nodes:         %d\n", node_count);
    printf("Voltage sources:     %d\n", m2_count);

    printf("\n--- OPTIONS ---\n");
    printf("SPARSE:              %s\n", use_sparse ? "YES" : "NO");
    printf("SPD (Cholesky):      %s\n", use_cholesky ? "YES" : "NO");
    printf("CUSTOM:              %s\n", use_custom ? "YES" : "NO");
    printf("ITERATIVE:           %s\n", use_iterative ? "YES" : "NO");
    printf("ITOL:                %g\n", itol);
    printf("TRANS:               %s  time step: %g final time: %g\n", do_trans ? "YES" : "NO", tran.time_step, tran.final_time);
    printf("TRANS SOURCES:       %d\n", tran_src_number);
    printf("==============================\n\n");
}

void print_transient_spec(TwoPortsSource *s)
{
    if (s->tr_type == TR_NONE) {
        printf("      [DC source]\n");
        return;
    }

    switch(s->tr_type)
    {
        case TR_EXP:
            printf("      EXP (i1=%.6g i2=%.6g td1=%.6g tc1=%.6g td2=%.6g tc2=%.6g)\n",
                s->tr.exp.i1, s->tr.exp.i2, s->tr.exp.td1, s->tr.exp.tc1,
                s->tr.exp.td2, s->tr.exp.tc2);
            break;

        case TR_SIN:
            printf("      SIN (i1=%.6g ia=%.6g freq=%.6g td=%.6g df=%.6g phase=%.6g)\n",
                s->tr.sin.i1, s->tr.sin.ia, s->tr.sin.freq,
                s->tr.sin.td, s->tr.sin.df, s->tr.sin.phase);
            break;

        case TR_PULSE:
            printf("      PULSE (i1=%.6g i2=%.6g td=%.6g tr=%.6g tf=%.6g pw=%.6g per=%.6g)\n",
                s->tr.pulse.i1, s->tr.pulse.i2, s->tr.pulse.td,
                s->tr.pulse.tr, s->tr.pulse.tf, s->tr.pulse.pw, s->tr.pulse.per);
            break;

        case TR_PWL:
            printf("      PWL ");
            for(int i = 0; i < s->tr.pwl.points; i++)
                printf("(%.6g %.6g) ", s->tr.pwl.t[i], s->tr.pwl.i[i]);
            printf("\n");
            break;

        default:
            printf("      [Unknown transient type]\n");
    }
}




void print_transient_info(TRAN_Analysis *tran, Node *nodes) {

    if(tran->plot_count == 0){
        printf("\n[TRAN] No transient plot nodes found.\n");
        return;
    }

    printf("\n================ TRAN ANALYSIS SETTINGS ================\n");
    printf("Simulation time step : %le\n", tran->time_step);
    printf("Final time           : %le\n", tran->final_time);
    printf("Nodes selected for plotting (%d):\n", tran->plot_count);
    printf("Method = %d, 1=tr, 0=be\n", tran->method);

    for(int i = 0; i < tran->plot_count; i++){
        int node_id = tran->plot_nodes[i];

        printf("  • V(%s)   (node ID: %d)\n", nodes[node_id].name, node_id);
    }

    printf("========================================================\n\n");
}


void print_elements(Element *elements, Element *m2_elements, int m2_size, int size){
    for (int i = 0; i < size; i++){
        switch (elements[i].type) {
            case RES:
            case CAP:
            case IND: {
                TwoPortsElement *tp = (TwoPortsElement *)elements[i].data;
                printf("    - %s (value: %.6g)\n", tp->name, tp->value);
                printf("connected with : %d %d nodes\n\n", (tp->ports[0]), tp->ports[1]);
                break;
            }
            case VSRC:
            case ISRC: {
                TwoPortsSource *tp = (TwoPortsSource *)elements[i].data;
                printf("    - %s (value: %.6g)\n", tp->name, tp->value);
                printf("connected with : %d %d nodes\n\n", tp->ports[0], tp->ports[1]);
                print_transient_spec(tp);
                break;
            }
            case DIODE: {
                Diode *d = (Diode *)elements[i].data;
                printf("    - Diode %s (model: %s, area: %.6g)\n", d->name, d->model_name, d->area);
                printf("connected with : %d %d nodes\n\n", d->ports[0], d->ports[1]);
                break;
            }
            case MOS: {
                Mos *m = (Mos *)elements[i].data;
                printf("    - MOS %s (model: %s, L=%.6g, W=%.6g)\n", m->name, m->model_name, m->l, m->w);
                printf("connected with : %d %d %d %d nodes\n\n", m->ports[0], m->ports[1], m->ports[2], m->ports[3]);
                break;
                }
            case BJT: {
                Bjt *q = (Bjt *)elements[i].data;
                printf("    - BJT %s (model: %s, area:%.6g )\n", q->name, q->model_name, q->area);
                break;
            }
            default:
                printf("    - Unknown element type %d\n", elements[i].type);
        }
    }

    for (int i = 0; i < m2_size; i++){
        switch (m2_elements[i].type) {
            case IND: {
                TwoPortsElement *tp = (TwoPortsElement *)m2_elements[i].data;
                printf("    - %s (value: %.6g)\n", tp->name, tp->value);
                printf("connected with : %d %d nodes\n\n", (tp->ports[0]), tp->ports[1]);
                break;
            }
            case VSRC: {
                TwoPortsSource *tp = (TwoPortsSource *)m2_elements[i].data;
                printf("    - %s (value: %.6g)\n", tp->name, tp->value);
                printf("connected with : %d %d nodes\n\n", tp->ports[0], tp->ports[1]);
                break;
            }
            default:
                printf("    - Unknown element type %d\n", m2_elements[i].type);
        }
    }
}

void print_all_nodes(NodeHashtable **table, int table_size) {
    printf("\n Node Table (Hash Table) \n");

    for (int i = 0; i < table_size; i++) {
        NodeHashtable *curr = table[i];
        if (!curr)
            continue; // empty bucket

        printf("\nBucket %d:\n", i);

        while (curr) {
            printf("Node: %-10s (ID: %d)\n", curr->name, curr->id);     
            curr = curr->next;
        }
    }

    printf("\n\n");
}

void print_elementIDs(ElementID **table, int table_size) {
    for (int i = 0; i < table_size; i++) {
        ElementID *curr = table[i];
        if (curr != NULL) {
            printf("Bucket %d:\n", i);
            while (curr != NULL) {
                printf("  Name: %s, Index: %d\n", curr->name, curr->index);
                curr = curr->next;
            }
        }
    }
}

void print_DCTables(gsl_matrix *A, gsl_vector *b, int DCtable_size) {
    for (int i = 0; i < DCtable_size; i++) {
        for (int j = 0; j < DCtable_size; j++) {
            printf("%.6g ", gsl_matrix_get(A, i, j));
        }
        printf("   b= %.6g\n", gsl_vector_get(b, i));
    }
}

void print_solution(const gsl_vector *x) {
    printf("\n--- Solution vector x ---\n");
    for (size_t i = 0; i < x->size; i++) {
        printf("x[%zu] = %.6g\n", i, gsl_vector_get(x, i));
    }
}

void print_nodes_by_id(Node *nodes, int node_count)
{
    printf("\nNodes by id:\n");
    for(int i=0; i<node_count; i++) {
        printf("Node id: %d, name: %s\n", nodes[i].id, nodes[i].name);
    }
}

void print_tran_sources(TranSrc *ts, int count){
    printf("\n====== TRANSIENT SOURCES ======\n");
    printf("Total: %d\n", count);

    for(int i=0;i<count;i++){
        printf("[%d] index=%d  type=%d\n", 
               i, ts[i].index, ts[i].type);
    }

    printf("===============================\n\n");
}