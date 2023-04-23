#ifndef AUXILIARYFUNCTION_H_INCLUDED
#define AUXILIARYFUNCTION_H_INCLUDED
#include <cstddef>

/* ******************************************************************** */
/* Auxiliary defined function Macros                                          */
/* ******************************************************************** */

#define strSIGN(x) (((x) >= (0)) ? "+":"-")
#define SIGN(x) (((x) >= (0)) ?1:-1)
#define ABS(x) ((SIGN(x))*(x))

/* ******************************************************************** */
/* Some constants                                                       */
/* ******************************************************************** */
#define myPI     3.14159265359
#define invPI    0.31830988618
#define SQRT2    1.41421356237
#define invSQRT2 0.70710678119


/* ******************************************************************** */

/* Some constants defined for use in code                               */
/* ******************************************************************** */
#define ACCURACY 1.0e-10
#define SET_PRECISION 15 //for printing
/* ******************************************************************** */
/* Useful functions: GIVEN                                              */
/* ******************************************************************** */
double Power(double x, long n); //Calculates x^n
/* Absolute Value Functions                                             */
long double ldAbs(long double x);
double dAbs(double x);
float fAbs(float x);
int iAbs(int x);
short sAbs(short x);
long lAbs(long x);
long long llAbs(long long x);
/* Swap functions                                                       */
void ldSwap(long double *x, long double *y);
void dSwap(double *x, double *y);
void fSwap(float *x, float *y);
void llSwap(long long *x, long long *y);
void lSwap(long *x, long *y);
void iSwap(int *x, int *y);
void Swap(void *x, void *y, size_t n); //assumes x and y of same size and type
/* ******************************************************************** */
/* Transcendental functions: TO BE IMPLEMENTED IN LAB003                */
double Exp(double);
double Sine(double);
double Cosine(double);
double Tangent(double);
double CoTangent(double);
double Cosec(double);
double Secant(double);
double Cosh(double);
double Sinh(double);
double Tanh(double);
double Cosech(double);
double Sech(double);
double Coth(double);
#endif // AUXILIARYFUNCTION_H_INCLUDED
