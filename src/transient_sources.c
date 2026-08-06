#define _USE_MATH_DEFINES   // (Windows only)
#include <math.h>
#include "structs.h"
#include "eval.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

double eval_EXP(double t, TransientEXP *e){
    if (t < e->td1) 
        return e->i1;

    else if (t < e->td2) 
        return e->i1 + (e->i2 - e->i1)*(1 - exp(-(t-e->td1)/e->tc1));

    return e->i1 + (e->i2 - e->i1)*(exp(-(t-e->td2)/e->tc2) - exp(-(t-e->td1)/e->tc1));
}


double eval_SIN(double t, TransientSIN *s){
    double tp = t - s->td;
    if (t < s->td) 
        return s->i1 + s->ia * sin( s->phase*M_PI/180.0);

    return s->i1 + s->ia * sin(2*M_PI*s->freq*tp + s->phase*M_PI/180.0) * exp(-s->df*tp);
}

double eval_PULSE(double t, TransientPULSE *p){
    double T = p->per;
    double tt = fmod(t - p->td, T);

    if (t < p->td) return p->i1;
    if (tt < p->tr) return p->i1 + (p->i2-p->i1)*(tt/p->tr);
    if (tt < p->tr + p->pw) return p->i2;
    if (tt < p->tr + p->pw + p->tf) return p->i2 - (p->i2-p->i1)*((tt-p->tr-p->pw)/p->tf);
    return p->i1;
}

double eval_PWL(double t, TransientPWL *p){
    if (p->points == 1)
        return p->i[0];

    if (t <= p->t[0])
        return p->i[0];

    for (int k = 0; k < p->points - 1; k++) {
        if (t < p->t[k+1])
            return p->i[k];
    }

    return p->i[p->points - 1];
}

/*
double eval_PWL(double t, TransientPWL *p){
    if (p->points == 1) return p->i[0];
    if (t <= p->t[0]) return p->i[0];
    for(int k=0;k<p->points-1;k++){
        if (t >= p->t[k] && t <= p->t[k+1]){
            double ratio = (t-p->t[k])/(p->t[k+1]-p->t[k]);
            return p->i[k] + ratio*(p->i[k+1]-p->i[k]);
        }
    }
    return p->i[p->points-1];
}*/

double eval_source(TwoPortsSource *src, double t){
    switch(src->tr_type){
        case TR_EXP:   return eval_EXP(t,&src->tr.exp);
        case TR_SIN:   return eval_SIN(t,&src->tr.sin);
        case TR_PULSE: return eval_PULSE(t,&src->tr.pulse);
        case TR_PWL:   return eval_PWL(t,&src->tr.pwl);
        default:       return src->value; // DC only
    }
}
