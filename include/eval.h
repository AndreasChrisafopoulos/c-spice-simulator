#ifndef EVAL_H
#define EVAL_H

#include "structs.h"

// === Individual waveform evaluators ===
double eval_EXP(double t, TransientEXP *e);
double eval_SIN(double t, TransientSIN *s);
double eval_PULSE(double t, TransientPULSE *p);
double eval_PWL(double t, TransientPWL *p);

// === Main dispatcher ===
double eval_source(TwoPortsSource *src, double t);

#endif
