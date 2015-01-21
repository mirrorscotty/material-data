/**
 * @file glass-transition.h
 */

#ifndef GLASS_TRANSITION_H
#define GLASS_TRANSITION_H

/**
 * Values for the Gordon-Taylor equation
 */
typedef struct {
    double Tg1;
    double Tg2;
    double kGT;
} gordontaylor;

gordontaylor* GTSemolina();
void DestroyGT(gordontaylor*);
double GordonTaylor(gordontaylor*, double);
double GordonTaylorInv(gordontaylor*, double);

#endif

