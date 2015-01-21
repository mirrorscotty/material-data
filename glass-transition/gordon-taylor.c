#include "glass-transition.h"
#include <stdlib.h>

gordontaylor* GTSemolina()
{
    gordontaylor *gt;
    gt = (gordontaylor*) calloc(sizeof(gordontaylor), 1);

    gt->Tg1 = 435;
    gt->Tg2 = 138;
    gt->kGT = 3.4;

    return gt;
}

void DestroyGT(gordontaylor *gt)
{
    free(gt);
    return;
}

double GordonTaylor(gordontaylor *gt, double Xdb)
{
    double w1, w2, Tg;
    w2 = Xdb/(1+Xdb);
    w1 = 1-w2;

    Tg = (w1*gt->Tg1 + gt->kGT*w2*gt->Tg2)/(w1+gt->kGT*w2);
    return Tg;
}

double GordonTaylorInv(gordontaylor *gt, double T)
{
    double Xwb;
    Xwb = (gt->Tg1-T)/(gt->Tg1-gt->Tg2*gt->kGT+(gt->kGT-1)*T);
    return Xwb/(1-Xwb);
}

