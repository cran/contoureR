#include <Rcpp.h>
#include "structs.h"

//Namespaces
using namespace Rcpp;
using namespace std;

#ifndef __functions__
#define __functions__

//Function declarations
Vec3 crossProd(Vec3 a, Vec3 b, Vec3 c);

void getVecsByRefZeroBased(const IntegerMatrix& dm, const NumericMatrix& xyz, const int& ix, Vec3 *out);

bool isEqual(double a, double b);
bool isIn(NumericVector& x, double& v);
bool isIn(vector<double>& x,double& v);

//Operators
IntegerMatrix operator+=(const IntegerMatrix& a, const int& b);
IntegerMatrix operator-=(const IntegerMatrix& a, const int& b);

#endif // __functions__
