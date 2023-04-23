#ifndef POLYNOMIAL_H_INCLUDED
#define POLYNOMIAL_H_INCLUDED
#include <iostream>
#include <vector>
class Polynomial{
std::vector<double> C; //Stores the coefficients
    void Resize(long);
public:
    //Constructors
Polynomial(); //creates a zero polynomial
Polynomial(double Data[], unsigned int Size); //creates a polynomial from the array
Polynomial(long Size); // creates a zero polynomial of size given
Polynomial(const Polynomial&);
Polynomial(double,long position); 
//creates a polynomial with the double coefficient at exponent position
/* Auxiliary functions */
    void Compact(); //removes highest leading zero coefficients
    long Degree(); //returns the degree of the polynomial
    bool isZero() {return Degree() == 0 && C[0] == 0;} // true if polynomial is zero else false
/* Arithmetic functions */
    friend Polynomial& operator += (Polynomial&,Polynomial );
    friend Polynomial operator + (Polynomial, Polynomial);
    friend Polynomial& operator -= (Polynomial&,Polynomial );
    friend Polynomial operator - (Polynomial, Polynomial);
    friend Polynomial& operator *= (Polynomial&,double );
    friend Polynomial operator * (Polynomial, double);
    friend Polynomial operator * (double,Polynomial);
    friend Polynomial& operator /= (Polynomial&,double );
    friend Polynomial operator / (Polynomial, double);
    friend Polynomial& operator *= (Polynomial&,Polynomial );
    friend Polynomial operator * (Polynomial, Polynomial);
    friend Polynomial& operator /= (Polynomial&,Polynomial );
    friend Polynomial operator / (Polynomial, Polynomial);
    friend Polynomial& operator %= (Polynomial&,Polynomial );
    friend Polynomial operator % (Polynomial, Polynomial);
/* Comparison functions
A polynomial with the larger degree is large
If degrees are equal, then the polynomial with the higher coefficient is larger */
    friend bool operator == (Polynomial&, Polynomial& ) ;
    friend bool operator != (Polynomial&, Polynomial&);
    friend bool operator >  (Polynomial&, Polynomial&);
    friend bool operator <  (Polynomial&, Polynomial&);
    friend bool operator >= (Polynomial&, Polynomial&);
    friend bool operator <= (Polynomial&, Polynomial&);
 /* Assignment operator */
    Polynomial& operator = (const Polynomial& x);
/* Horner's Rules */
    double Eval(double);
    double operator () (double x) { return Eval(x); }
    Polynomial Deflate(double& x);
        //deflate polynomial about double. The x carries back the original polynomial value.
        //The deflated polynomial is returned, original remains the same
    Polynomial Derivative();
    Polynomial Integrate(double fx0 = 0.0); 
// The constant of integration is the double specified
/* I/O Functions */
    void PolyPrint();
    friend std::istream& operator>> (std::istream& , Polynomial& P);
    friend std::ostream& operator<< (std::ostream& , const Polynomial&  P);
        // reads the polynomial as a sequence of (coefficient) numbers spaced by spaces and
        // terminated by newline.
        // The zeroth coefficient term is first, every degree term coefficient is specified.
        // My implementation assumes that if the sign of the coefficient is specified, there
        // is no gap between the sign and the number. CAN YOU DO BETTER
};
#endif // POLYNOMIAL_H_INCLUDED
