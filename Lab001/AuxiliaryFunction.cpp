#include "AuxiliaryFunction.h"
#include <cstring>
/* ******************************************************************** */
/* Useful functions                                                     */
/* ******************************************************************** */
double Power(double x, long n) {
    if (n < 0) {x = 1.0/x; n = -n;}
    if (n == 0) return 1.0;
    else {
        double y = 1.0;
        while(n) {
            if ((n%2) == 1) y *= x;
            x *= x;
            n >>= 1;
        }
        return y;
    }
}
/* AUXILIARY FUNCTIONS                                                                                */
long double ldAbs(long double x) { return x* (SIGN(x)); }
double dAbs(double x) {return (SIGN(x)) * x;}
float fAbs(float x) {return (SIGN(x)) * x;}
int iAbs(int x) {return (SIGN(x)) * x;}
short sAbs(short x) {return (SIGN(x)) * x;}
long lAbs(long x) {return (SIGN(x)) * x;}
long long llAbs(long long x) {return (SIGN(x)) * x;}
void ldSwap(long double *x, long double *y) { long double t = *x; *x = *y; *y = t;}
void dSwap(double *x, double *y) { double t = *x; *x = *y; *y = t;}
void fSwap(float *x, float *y) { float t = *x; *x = *y; *y = t;}
void llSwap(long long *x, long long *y) { long long t = *x; *x = *y; *y = t;}
void lSwap(long *x, long *y) { long t = *x; *x = *y; *y = t;}
void iSwap(int *x, int *y) { int t = *x; *x = *y; *y = t;}
void Swap(void *x, void *y, size_t n){
    void *t = ::operator new(n); //ideally check for allocation error
    memmove(t,x,n);
    memmove(x,y,n);
    memmove(y,t,n);
::operator delete(t);
}
