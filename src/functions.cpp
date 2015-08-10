#include <Rcpp.h>
#include <limits>
#include <math.h>
#include <algorithm>
#include "constants.h"          //Global Constants
#include "functions.h"
#include "structs.h"

//Namespaces
using namespace Rcpp;
using namespace std;

bool isEqual(double a, double b){ 
  return abs(a-b) <= D_TOL; 
}
bool isIn(NumericVector& x, double& v){
  return std::find(x.begin(), x.end(), v) != x.end();
}
bool isIn(vector<double>& x, double& v){
  return std::find(x.begin(), x.end(), v) != x.end();
}

void getVecsByRefZeroBased(const IntegerMatrix& dm, const NumericMatrix& xyz, const int& ix, Vec3 * out){ 
  for(int j = 0, k = 0; j < 3; j++){
    k = dm(ix,j);
    out[j] = Vec3(xyz(k,0),xyz(k,1),xyz(k,2)); 
  }
}

Vec3 crossProd(Vec3 a, Vec3 b, Vec3 c){
  Vec3 result, B = b - a, C = c - a;
  result.x = B.y*C.z - B.z*C.y;
  result.y = B.z*C.x - B.x*C.z;
  result.z = B.x*C.y - B.x*C.y;
  return result;
}

IntegerMatrix operator+=(const IntegerMatrix& a, const int& b){
    IntegerMatrix r = a;
    for(int i = 0; i < a.ncol(); i++){ 
      for(int j = 0; j < a.nrow(); j++){ 
        r(j,i) += b; 
      } 
    }
    return r;
}

IntegerMatrix operator-=(const IntegerMatrix& a, const int& b){
  return a += (-b);
}











