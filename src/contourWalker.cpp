#include <Rcpp.h> 
#include <math.h>
#include <time.h> 
#include <map>
#include <vector>
#include <limits>

#include "constants.h"          //Global Constants
#include "functions.h"          //Functions
#include "structs.h"
#include "mtrand.h"

//Namespaces
using namespace Rcpp;
using namespace std;

//Forward Delcarations
NumericVector           checkAndUniqueLevels(NumericVector& input);
IntegerMatrix           checkAndSortDM(IntegerMatrix& dm, NumericMatrix& xyz);
NumericMatrix           pertubate(IntegerMatrix& dm, NumericMatrix& xyz, NumericVector& levels, double& upToPercent);

inline Vec3 interpolateZVal(Vec3& a, Vec3& b, double& z){ Vec3 c = b - a; c *= (z-a.z)/c.z; c += a; c.z = z; return c; }

// [[Rcpp::export]]
NumericMatrix contourWalker(IntegerMatrix& dm, NumericMatrix& xyz, NumericVector& levels, double maximumPertubation=1e-5){
  //Create Local Copies & Variables
  IntegerMatrix dmNew     = IntegerMatrix(dm);  
  NumericMatrix xyzNew    = NumericMatrix(xyz); 
  NumericVector levelsNew = NumericVector(levels); 
  
  //Run Checks
  dmNew         = checkAndSortDM(dmNew,xyzNew);
  levelsNew     = checkAndUniqueLevels(levelsNew);
  xyzNew        = pertubate(dmNew,xyzNew,levelsNew,maximumPertubation);
  
  //Variables
  vector<Vec3> contour;
  vector<ContourData> contourData;
  unsigned int i, levelID, groupID, delID, pathID;
  
  //Build Nodes
  vector<Node>  nodes;
  for(i = 0; i < xyzNew.nrow(); i++)
    nodes.push_back(Node(i,Vec3(xyzNew(i,0), xyzNew(i,1), xyzNew(i,2))));
  
  //Build Deleaunay Triangles
  vector<Del>   dels;
  for(i = 0; i < dmNew.nrow(); i++)
    dels.push_back(Del(nodes[dmNew(i,0)],nodes[dmNew(i,1)],nodes[dmNew(i,2)]));
  
  //Establish the Network
  for(i = 0; i < dels.size() - 1; i++){
    for(int j = i+1; !dels[i].isFull() && j < dels.size(); j++){
      if(!dels[j].isFull())
        (&(dels[i]))->makeNeighbour( &(dels[j]) ); //TRY
    }
  }
  
  //Draw the Contours and pass the result back to the contour container by reference
  for(levelID = 0; levelID < levelsNew.size(); levelID++){
      for(delID = 0, groupID = 0; delID < dels.size(); delID++){
        (&(dels[delID]))->drawContour(levelsNew[levelID],contour,NULL);
        for(pathID = 0; pathID < contour.size(); pathID++)
          contourData.push_back(ContourData(levelID,groupID,pathID,contour[pathID].x,contour[pathID].y,contour[pathID].z));
        groupID += (contour.size() > 0) ? 1 : 0;
        contour.clear();
      }
      for(delID = 0; delID < dels.size(); delID++){
        dels[delID].reset();
      }
  }
  
  //Pre-allocate the result array
  NumericMatrix result(contourData.size(),6);
  
  //Put the contourData vector into numeric matrix and return the result 
  for(pathID = 0; pathID < contourData.size(); pathID++)
    result.row(pathID) = (&contourData[pathID])->toNumericVector();
  
  //Done
  return result;
}


NumericVector checkAndUniqueLevels(NumericVector& input){
  if(input.size() == 0)
    throw std::out_of_range("The Levels Vector is Empty, please specify at least one level to contour.");
  vector<double> tmp = vector<double>();
  for(int i = 0; i < input.size(); i++)
    if(!isIn(tmp,input(i))){ tmp.push_back(input(i)); }
  NumericVector ret(tmp.begin(),tmp.end());
  return ret;
}

IntegerMatrix checkAndSortDM(IntegerMatrix& dm, NumericMatrix& xyz){
  if(dm.ncol() != 3 || xyz.ncol() != 3)
    throw std::out_of_range("Expecting 3 columns in both 'dm' and 'xyz' matrixes.");
  
  //Check 1-based indexing of the dm matrix (coming from R)
  if(min(dm) == 1)
    dm += -1;
  
  //Check values in dm, will not raise out of bounds exception in xyz
  if(max(dm) >= xyz.nrow())
    throw std::out_of_range("Values in 'dm' would result in out of bounds errors in xyz.");
  
  IntegerMatrix dmSorted(dm.nrow(),dm.ncol());
  Centroid* centroids = new Centroid[dm.nrow()]();
  Vec3* vecs = new Vec3[3]();
  
  //Assemble the centroids array, from the vectors
  for(int i = 0; i < dm.nrow(); i++){
     getVecsByRefZeroBased(dm,xyz,i,vecs);
     centroids[i] = Centroid(i); //To Sort
      for(int j = 0; j < 3; j++)
        centroids[i] += vecs[j]/3.0;
  }
  
  //Execute the sort using the comparitor which is part of the centroids struct
  std::sort(centroids,centroids + dm.nrow());
  
  //Restore the data
  for(int i = 0; i < dm.nrow(); i++)
    dmSorted.row(i) = dm.row(centroids[i].index);
  dm = dmSorted;
  
  //Cleanup the dynamic arrays
  delete [] centroids;
  delete [] vecs;
  
  //Done, Return
  return dm;
}


inline double nudge(double& input, double& upToPcnt){
  MTRand_closed drand;
  double  random1 = drand() - 1.0, 
          random2 = drand() - 1.0;
  return input*(1.0 + random1*abs(upToPcnt)/100.0) + 2*random2*D_TOL;
}

NumericMatrix pertubate(IntegerMatrix& dm, NumericMatrix& xyz, NumericVector& levels, double& upToPercent){
  int pb, i, limit = dm.nrow(), cnt=0;
  double pcnt = abs(upToPercent);
  while(cnt <= limit && abs(upToPercent - D_TOL) > 0){ 
    for(i = 0, pb = 0; i < dm.nrow(); i++){
      int    ix[3] = {dm(i,0),dm(i,1),dm(i,2)};
      double z[3]  = {xyz(ix[0],2),xyz(ix[1],2),xyz(ix[2],2)};
      if(isEqual(z[0],z[1])  || isEqual(z[0],z[2])  || isEqual(z[1],z[2]) || 
         isIn(levels,z[0])   || isIn(levels,z[1])   || isIn(levels,z[2])){
         xyz(ix[0],2) = nudge(z[0],pcnt);
         xyz(ix[1],2) = nudge(z[1],pcnt);
         xyz(ix[2],2) = nudge(z[2],pcnt);
         pb++; 
      }
    }; cnt++;
    if(pb == 0) break;
  }
  return xyz;
}




