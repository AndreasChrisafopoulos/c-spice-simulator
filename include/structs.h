#ifndef STRUCTS_H
#define STRUCTS_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <strings.h>  
#include "uthash.h"

#define MAX_LINE_LEN 256
#define MAX_NAME_LEN 32
#define MAX_SWEEPS 32
#define MAX_SWEEP_PLOTS 32

enum element {RES, CAP, IND, VSRC, ISRC, DIODE, MOS, BJT, ZERO};
typedef enum element elemenT;

typedef enum {TR_NONE, TR_EXP, TR_SIN, TR_PULSE, TR_PWL } TransientType;

typedef struct {
    double i1, i2, td1, tc1, td2, tc2;
} TransientEXP;

typedef struct {
    double i1, ia, freq, td, df, phase;
} TransientSIN;

typedef struct {
    double i1, i2, td, tr, tf, pw, per;
} TransientPULSE;

typedef struct {
    double t[20];
    double i[20];
    int points;
} TransientPWL;

typedef struct {
    elemenT type;
    void* data;
} Element;

typedef struct Node {
    char name[MAX_NAME_LEN];
    int id;
    struct Node *next;
} NodeHashtable;

typedef struct {
    char name[MAX_NAME_LEN];
    int ports[2];
    double value;
} TwoPortsElement;

typedef struct {
    char name[MAX_NAME_LEN];
    int ports[2];
    double value;
    TransientType tr_type;
    union {
        TransientEXP exp;
        TransientSIN sin;
        TransientPULSE pulse;
        TransientPWL pwl;
    } tr;

} TwoPortsSource;

typedef struct {
    char name[MAX_NAME_LEN];
    int ports[2];
    char model_name[MAX_NAME_LEN];
    double area;
} Diode;

typedef struct {
    char name[MAX_NAME_LEN];
    int ports[4];
    char model_name[MAX_NAME_LEN];
    double l, w;
} Mos;

typedef struct {
    char name[MAX_NAME_LEN];
    int ports[3];
    char model_name[MAX_NAME_LEN];
    double area;
} Bjt;

typedef struct ElementID {
    char name[MAX_NAME_LEN];
    int index;
    struct ElementID *next;
} ElementID;

typedef struct NodeSortbyID {
    char name[MAX_NAME_LEN];
    int id;
} Node;

typedef struct {
    char name[64];          // node name (key)
    UT_hash_handle hh;      // uthash handle
} NodeEntry;

typedef struct {
    char source_name[MAX_NAME_LEN];
    double start;
    double end;
    double step;
    int plot_nodes[MAX_SWEEP_PLOTS];   // node indices
    int plot_count;
} DCSweep;

typedef struct {
    double time_step;
    double final_time;
    int plot_nodes[MAX_SWEEP_PLOTS];
    int plot_count;
    int method;
} TRAN_Analysis;

typedef struct {
    int index;     // index  ( elements[] or m2[])
    int type;      // VSRC/ISRC/SIN/PWL/...
} TranSrc;






#endif