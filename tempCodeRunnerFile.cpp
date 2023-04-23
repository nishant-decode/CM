#include <cstddef>
#include <cstring>
#include <windows.h>
#include <exception>
#include <string>
#include <iostream>
#include <vector>
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <iomanip>
#include <iostream>

using namespace std;
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
// AUXILIARYFUNCTION_H_INCLUDED
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

class myTimer{
    LARGE_INTEGER F; //Frequency
    LARGE_INTEGER sT; // Start Time
    LARGE_INTEGER eT; // End Time
    double interval; // End Time - Start Time in seconds
public:
myTimer() { QueryPerformanceFrequency(&F); }
    void StartTimer(){ QueryPerformanceCounter(&sT); }
    void EndTimer(){ QueryPerformanceCounter(&eT); }
    double GetInterval() {
        return (double) (eT.QuadPart - sT.QuadPart) / F.QuadPart;
    }
};

class myException: public std::exception {
protected:
    /* Error message      */
std::string msg;
public:
    /* Constructor (C strings). */
    explicit myException(const char* message) : msg(message) {}
    /* Constructor (C++ STL strings). */
    explicit myException(const std::string& message) : msg(message) {}
    /* Destructor. Virtual to allow for subclassing if required. */
    virtual ~myException() noexcept {}
    /* Returns a pointer to the (constant) error description. */
    virtual const char* what() const noexcept { return msg.c_str(); }
};

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
        // is no gap between the sign and the number. 
};

/* *********************************************************** */
/* PRIVATE METHODS                                             */
/* *********************************************************** */
void Polynomial::Resize(long Size){ C.resize(Size); }

/* *********************************************************** */
/* CONSTRUCTORS                                                */
/* *********************************************************** */

Polynomial::Polynomial(){ C.push_back(0.0);}

Polynomial::Polynomial(double Data[], unsigned int Size){
    for (unsigned int i = 0; i < Size; i++) C.push_back(Data[i]);
Compact();
}

Polynomial::Polynomial(const Polynomial& A){
    C = A.C;
}

Polynomial::Polynomial(long Size){
    C.resize(Size);
}

Polynomial::Polynomial(double a,long position){
    C.resize(position + 1);
    C.at(position) = a;
}

/* *********************************************************** */
/* Some Auxiliary Functions                                    */
/* *********************************************************** */
void Polynomial::Compact(){ //removes leading zero coefficients
    for (long i = C.size() - 1; i>=0; i--) {
        if (C[i] == 0) C.pop_back();
        else break;
    }
    if (C.empty()) C.push_back(0.0);
}

long Polynomial::Degree(){ return C.size() - 1; }
/* *********************************************************** */
/* I/O overloading                                             */
/* *********************************************************** */
std::ostream& operator<< (std::ostream& os, const Polynomial&  P){
    if (P.C.size() == 0) return os << "0" << std::endl;
    long i = P.C.size() - 1;
    while ( i> 0 && P.C[i] == 0) i--; //remove highest 0 terms
    if (i == 0) return os << P.C[i];
    else {
        if (ABS(P.C[i]) == 1 && SIGN(P.C[i]) == 1) os << "x^(" << i << ") ";
        else if (ABS(P.C[i]) == 1) os << "- x^(" << i << ") ";
        //else if (SIGN(P.C[i]) == -1) os << string(strSIGN(P.C[i])) << " " << ABS(P.C[i]) << " x^(" << i << ") ";
        else os <<std::setprecision(SET_PRECISION) << P.C[i] << " x^(" << i << ") ";
    }
    i--;
    for (; i >= 1;i--) {
        if (ABS(P.C[i]) == 0) continue;
        if (ABS(P.C[i]) == 1) os <<std::string(strSIGN(P.C[i])) << " x^(" << i << ") ";
        else os <<std::string(strSIGN(P.C[i])) << " "
<<std::setprecision(SET_PRECISION) << ABS(P.C[i]) << " x^(" << i << ") ";
    }
    if (P.C[i] == 0) return os;
    else return os <<std::string(strSIGN(P.C[i])) << " "
<<std::setprecision(SET_PRECISION) << ABS(P.C[0]);
}

void Polynomial::PolyPrint(){
    if (C.size() == 0) {std::cout << "0" << std::endl; return;}
    long i = C.size() - 1;
    while ( i> 0 && C[i] == 0) i--; //remove highest 0 terms
    if (i == 0) {std::cout << C[i] << std::endl; return;}
    else {
        if (ABS(C[i]) == 1 && SIGN(C[i]) == 1) std::cout << "x^(" << i << ") ";
        else if (ABS(C[i]) == 1) std::cout << "- x^(" << i << ") ";
        //else if (SIGN(C[i]) == -1) cout << string(strSIGN(C[i])) << " " << ABS(C[i]) << " x^(" << i << ") ";
        else std::cout << C[i] << " x^(" << i << ") ";
    }
    i--;
    for (; i >= 1;i--) {
        if (ABS(C[i]) == 0) continue;
        if (ABS(C[i]) == 1) std::cout << std::string(strSIGN(C[i])) << " x^(" << i << ") ";
        else std::cout << std::string(strSIGN(C[i])) << " " << ABS(C[i]) << " x^(" << i << ") ";
    }
    if (C[i] == 0) std::cout << std::endl;
    else std::cout << std::string(strSIGN(C[i])) << " " << ABS(C[0]) << std::endl;
}

std::istream& operator>> (std::istream& Input, Polynomial &P){
    P.C.clear();
std::string aLine;
    //read a line
    if (getline(Input, aLine)){
        /* split line into words called myNumber
           Remember the signs should be together with the number if present
           simplistic
        */
std::stringstream tmpstream(aLine); //requires #include <sstream>
        for (std::string myNumber; tmpstream >> myNumber; ) P.C.push_back(std::stod(myNumber));
    }
    return Input;
}

/* ************************************************************************** */
/* Assignment operator                                                        */
/* ************************************************************************** */
Polynomial&Polynomial::operator = (const Polynomial& x) {
        C = x.C;
        return *this;
}  
/* *********************************************************** */
/* Arithmetic operators: TO BE IMPLEMENTED                     */
/* *********************************************************** */
Polynomial &operator+=(Polynomial &a, Polynomial b)
{
    // TO BE IMPLEMENTED BY STUDENT
    long maxSize = std::max(a.C.size(), b.C.size());
    a.C.resize(maxSize, 0.0);
    for (long i = 0; i < b.C.size(); i++)
        a.C[i] += b.C[i];
    a.Compact();
    return a;   
}

Polynomial operator+(Polynomial a, Polynomial b)
{
    // TO BE IMPLEMENTED BY STUDENT
    a += b;
    return a;
}

Polynomial &operator-=(Polynomial &a, Polynomial b)
{
    // TO BE IMPLEMENTED BY STUDENT
    long maxSize = std::max(a.C.size(), b.C.size());
    a.C.resize(maxSize, 0.0);
    for (long i = 0; i < b.C.size(); i++)
        a.C[i] -= b.C[i];
    a.Compact();
    return a;
}

Polynomial operator-(Polynomial a, Polynomial b)
{
    // TO BE IMPLEMENTED BY STUDENT
    a -= b;
    return a;
}

Polynomial &operator*=(Polynomial &a, double b)
{
    // TO BE IMPLEMENTED BY STUDENT
    for (long i = 0; i < a.C.size(); i++)
        a.C[i] *= b;
    a.Compact();
    return a;
}

Polynomial operator*(Polynomial a, double b)
{
    // TO BE IMPLEMENTED BY STUDENT
    a *= b;
    return a;
}

Polynomial operator*(double b, Polynomial a)
{
    // TO BE IMPLEMENTED BY STUDENT
    a *= b;
    return a;
}

Polynomial &operator/=(Polynomial &a, double b)
{
    // TO BE IMPLEMENTED BY STUDENT
    if (b == 0)
    {
        throw std::runtime_error("Divide by zero error: Polynomial division by zero");
    }
    for (long i = 0; i < a.C.size(); i++)
        a.C[i] /= b;
    a.Compact();
    return a;
}

Polynomial operator/(Polynomial a, double b)
{
    // TO BE IMPLEMENTED BY STUDENT
    a /= b;
    return a;
}

Polynomial &operator*=(Polynomial &a, Polynomial b)
{
    // TO BE IMPLEMENTED BY STUDENT
    long maxSize = a.C.size() + b.C.size() - 1;
    std::vector<double> result(maxSize, 0.0);
    for (long i = 0; i < a.C.size(); i++)
    {
        for (long j = 0; j < b.C.size(); j++)
        {
            result[i + j] += a.C[i] * b.C[j];
        }
    }
    a.C = result;
    a.Compact();
    return a;
}

Polynomial operator*(Polynomial a, Polynomial b)
{
    // TO BE IMPLEMENTED BY STUDENT
    a *= b;
    return a;
}

Polynomial &operator/=(Polynomial &a, Polynomial b)
{
    // TO BE IMPLEMENTED BY STUDENT
    // Check for division by zero
    if (b.isZero())
        throw std::invalid_argument("Divisor cannot be zero");

    // Divide the polynomials
    while (a.Degree() >= b.Degree() && !a.isZero())
    {
        double ratio = a.C[a.Degree()] / b.C[b.Degree()];
        long degree_diff = a.Degree() - b.Degree();

        // Subtract the quotient times the divisor from the dividend
        for (unsigned int i = 0; i <= b.Degree(); i++)
        {
            a.C[degree_diff + i] -= ratio * b.C[i];
        }

        a.Compact();
    }

    return a;
}

Polynomial operator/(Polynomial a, Polynomial b)
{
    // TO BE IMPLEMENTED BY STUDENT
    a /= b;
    return a;
}

Polynomial &operator%=(Polynomial &a, Polynomial b)
{
    // TO BE IMPLEMENTED BY STUDENT
    long aDegree = a.Degree();
    long bDegree = b.Degree();
    while (aDegree >= bDegree) {
        double coeff = a.C[aDegree] / b.C[bDegree];
        long degDiff = aDegree - bDegree;
        Polynomial term(coeff, degDiff);
        a -= (b * term);
        aDegree = a.Degree();
    }
    a.Compact();
    return a;

}

Polynomial operator%(Polynomial a, Polynomial b)
{
    // TO BE IMPLEMENTED BY STUDENT
    a %= b;
    return a;
}
/* ***************************************************************** */
/* Evaluation: TO BE IMPLEMENTED                                     */
/* ***************************************************************** */

double Polynomial::Eval(double x)
{
    // TO BE IMPLEMENTED BY STUDENT
    double result = 0.0;
    int n = Degree();

    for (int i = n; i >= 0; i--)
    {
        result = result * x + C[i];
    }

    return result;
}

// Polynomial Deflate(double& x)
// {
    
// }

Polynomial Polynomial::Derivative()
{
    // TO BE IMPLEMENTED BY STUDENT
    long n = Degree();
    if (n == 0)
        return Polynomial(0);
    Polynomial result(n);
    for (long i = 0; i < n; i++)
        result.C[i] = (i + 1) * C[i + 1];
    return result;
}

Polynomial Polynomial::Integrate(double x)
{
    // TO BE IMPLEMENTED BY STUDENT
    vector<double> resultCoeffs(C.size() + 1); // Store the coefficients of the integrated polynomial
    resultCoeffs[0] = x; // Set the constant of integration
    
    for (size_t i = 0; i < C.size(); i++)
    {
        resultCoeffs[i + 1] = C[i] / (i + 1); // Integrate each coefficient
    }
    
    Polynomial result;
    result.C = resultCoeffs;
    return result;
}

/* ******************************************************************* */
/* Relational operators : TO BE IMPLEMENTED                            */
/* ******************************************************************* */
bool operator == (Polynomial &a, Polynomial &b){ 
    // TO BE IMPLEMENTED BY STUDENT
    if (a.Degree() != b.Degree()) {
        return false;
    }
    for (long i = a.Degree(); i >= 0; i--) {
        if (a.C[i] != b.C[i]) {
            return false;
        }
    }
    return true;
}

bool operator != (Polynomial& a, Polynomial& b) { return !(a == b); }

bool operator > (Polynomial& a, Polynomial& b) { 
    // TO BE IMPLEMENTED BY STUDENT
if (a.Degree() > b.Degree()) {
        return true;
    }
    else if (a.Degree() < b.Degree()) {
        return false;
    }
    else {
        for (long i = a.Degree(); i >= 0; i--) {
            if (a.C[i] > b.C[i]) {
                return true;
            }
            else if (a.C[i] < b.C[i]) {
                return false;
            }
        }
    }
    return false;
}
bool operator < (Polynomial& a, Polynomial& b) { return !(a == b || a > b);};
bool operator >= (Polynomial& a, Polynomial& b) { return a == b || a > b;}
bool operator <= (Polynomial& a, Polynomial& b) { return !(a > b);}
int main()
{
    myTimer T;
T.StartTimer();
    cout << "CONSTRUCTORS-----POLYNOMIAL CONSTRUCTION\n";
{ // A zero polynomial
        Polynomial A;
        cout << "1. A = " << A << " A.Degree = " << A.Degree() << endl;
        double a[] = {0,0,0};
        Polynomial B(a,3);
        cout << "2. B = " << B << " B.Degree = " << B.Degree() << endl;
        // enter a sequence of zeros followed by ENTER
        //Polynomial C;
        //cin >> C;
        //cout << "3. C = " << C << " C.Degree = " << C.Degree() << endl;
        Polynomial D(10);
        cout << "4. D = " << D << " D.Degree = " << D.Degree() << endl;
        double b[] = {1,2,0,4};
        Polynomial E(b,4);
        cout << "5. E = " << E << " E.Degree = " << E.Degree() << endl;
        Polynomial F(20,7);
        cout << "6. F = " << F << " F.Degree = " << F.Degree() << endl;
    }
    cout << "-----POLYNOMIAL ADDITION------\n";
    {
        double a[] = {1,-1,3,-2};
        Polynomial A(a,4);
        double b[] = {2,1,-3,5};
        Polynomial B(b,4);
        cout << "7. A= " << A << endl
<< "   B = " << B << endl
<< "   A+B = " << (A+B) << endl;
        cout << "8. A+= B :" << (A+=B) << endl
<< "   A = " << A << endl;
    }
    cout << "-----POLYNOMIAL SUBTRACTION------\n";
    {
        double a[] = {1,-1,3,-2};
        Polynomial A(a,4);
        double b[] = {2,1,-3,5};
        Polynomial B(b,4);
        cout << "9. A= " << A << endl
<< "   B = " << B << endl
<< "   A-B = " << (A-B) << endl;
        cout << "10.A+= B :" << (A-=B) << endl
<< "   A = " << A << endl;
    }
        cout << "-----POLYNOMIAL MULTIPLICATION------\n";
    {
        double a[] = {1,1,1};
        Polynomial A(a,3);
        double b[] = {1,1};
        Polynomial B(b,2);
        cout << "11.A= " << A << endl
<< "   B = " << B << endl
<< "   A*B = " << (A*B) << endl;
        cout << "12.A*= B :" << (A*=B) << endl << "   A = " << A << endl;
    }
        cout << "-----POLYNOMIAL DIVISON/REMAINDER------\n";
    {
        double a[] = {2,0,3,4,5};
        Polynomial A(a,5);
        double b[] = {1,1,-2};
        Polynomial B(b,3);
        cout << "13.A= " << A << endl
<< "   B = " << B << endl
<< "   A/B = " << (A/B) << endl
<< "   A%B = " << (A%B) << endl;
        cout << "14.(A/B)*A + A%B = :" << (A/B)*B + A%B << endl;
    }
        cout << "-----POLYNOMIAL DIVISON BY ZERO------\n";
    {
        double a[] = {2,0,3,4,5};
        Polynomial A(a,5);
        double b[] = {0,0,-0};
        Polynomial B(b,3);
try{
            Polynomial C = B / A;
            cout << "15. C = " << C << endl;
            Polynomial D = A / B;
            cout << "16. D = " << D << endl;
        } catch (myException &e) {
            cout <<e.what() << endl;
        }
    }
    cout << "-----POLYNOMIAL REMAINDER BY ZERO------\n";
    {
        double a[] = {2,0,3,4,5};
        Polynomial A(a,5);
        double b[] = {0,0,-0};
        Polynomial B(b,3);
try{
            Polynomial C = B % A;
            cout << "15. C = " << C << endl;
            Polynomial D = A % B;
            cout << "16. D = " << D << endl;
        } catch (myException &e) {
            cout <<e.what() << endl;
        }
    }
        cout << "-----POLYNOMIAL EVALUATION ETC------\n";
    {
        double a[] = {1,2,1};
        Polynomial A(a,3);
        cout << "17.A = " << A << endl
<< "   A(1) = " << A(1) << " A(-1) = " << A(-1) << endl;
        Polynomial B;
        double x = 1;
        cout << "18.A = " << A << endl;
        B = A.Deflate(x);
        cout << "   Deflate A at -1: " << B << "; A(1) =" << x << endl
<< "   Deflated Poly  : " << B << endl;
        B = A;
        cout << "19. B = " << B << endl;
        B = B.Integrate(5);
        cout << "20. Integrate(B) = " << B << endl;
        B = B.Derivative();
        cout << "21. Diff(B) = " << B << endl;
    }
T.EndTimer();
    cout << "Time Taken = " <<T.GetInterval() << "(Seconds)" << endl;
    return 0;
}
/* ****************************************************************************** */

