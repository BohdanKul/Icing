#ifndef HYPERCUBE_H
#define HYPERCUBE_H

#include <iostream>
#include <iomanip>
#include <map>

using namespace std;

// direction map - association between a relative direction and a number
const map<char, int> DMap = {
                                {'r', 0},
                                {'d', 1},
                                {'l', 2},
                                {'u', 3},
                                {'n', 4},
                                {'p', 5}
                            };

// the inverse map - from a number to a relative direction
const map<int, char> IDMap = {
                                {0, 'r'},
                                {1, 'd'},
                                {2, 'l'},
                                {3, 'u'},
                                {4, 'n'},
                                {5, 'p'}
                             };

// opposite direction map - associates with a each direction its opposite
const map<char, char> ODMap = {
                                {'l', 'r'},
                                {'r', 'l'},
                                {'u', 'd'},
                                {'d', 'u'},
                                {'n', 'p'},
                                {'p', 'n'}
                              };

// the same as the previous map but in terms of numbers
const map<int, int>  RDMap  = {
                                {0, DMap.at(ODMap.at(IDMap.at(0)))},
                                {1, DMap.at(ODMap.at(IDMap.at(1)))},
                                {2, DMap.at(ODMap.at(IDMap.at(2)))},
                                {3, DMap.at(ODMap.at(IDMap.at(3)))},
                                {4, DMap.at(ODMap.at(IDMap.at(4)))},
                                {5, DMap.at(ODMap.at(IDMap.at(5)))}
                              };


/****************************************************************************************************
 * A class defining a lattice of a hypercube. Currently designed for 2 and 3 dimensional cases.   
 ***************************************************************************************************/
class HyperCube{
    private:
        int X; 
        int Y;
        int Z;
        int N;              // total number of spins 
        int S;              // number of spins on the boundary
        int dim;            // dimensionality of the cube
        char up; char down; // directions along which OBC are
        void GetFullCrd(int sindex, int& x, int& y, int& z); // convert a 1-d index to a 3-d coordinate
    
    public:
        HyperCube(int _X, int _Y, int _Z);
        char GetUp(){   return up;};
        char GetDown(){ return down;};
        int  GetSize(){ return N;};
        int  GetDim(){  return dim;};
        int  GetBoundarySize(){ return S;};

        vector<int> GetSpinNghbs(int sindex);  // get all the neighbours of a spin
        void GetBoundaryPair(int sindex, int& bspin, int& tspin); // get two boundary spins 
};

/****************************************************************************************************
 * A 3-d lattice class
 ***************************************************************************************************/
HyperCube::HyperCube(int _X, int _Y, int _Z){
    X = _X;
    Y = _Y;
    Z = _Z;
    N = X*Y*Z;
    if (Z !=1 ){ dim = 3; up = 'n'; down = 'p'; S = X*Y;}
    else       { dim = 2; up = 'u'; down = 'd'; S = X;}
}

/****************************************************************************************************
 * Get all the neighbours of a spin
 ***************************************************************************************************/
vector<int> HyperCube::GetSpinNghbs(int sindex){
    vector<int> nghbs;
    nghbs.resize(dim*2);
    
    int x; int y; int z;
    GetFullCrd(sindex, x, y, z);
    
    int siteB;
    // right neighbour
    siteB = z*X*Y + y*X + (x+1)%X; 
    nghbs[DMap.at('r')] = siteB;
    
    // left neighbour
    if (x==0) siteB = z*X*Y + y*X + X-1;
    else      siteB = z*X*Y + y*X + x-1; 
    nghbs[DMap.at('l')] = siteB;

    // top neighbour
    siteB = z*X*Y + ((y+1)%Y)*X + x; 
    nghbs[DMap.at('u')] = siteB;

    // bottom neighbour
    if (y==0) siteB = z*X*Y + (Y-1)*X + x;
    else      siteB = z*X*Y + (y-1)*X + x;
    nghbs[DMap.at('d')] = siteB;
 
    if (dim>2){
        // neighbour on the next depth layer
        siteB = ((z+1)%Z)*X*Y + y*X + x;
        nghbs[DMap.at('n')] = siteB;

        // neighbour on the previous depth layer
        if (z==0) siteB = (Z-1)*X*Y + y*X + x;
        else      siteB = (z-1)*X*Y + y*X + x;
        nghbs[DMap.at('p')] = siteB;
    }
    return nghbs;
}

/****************************************************************************************************
 * Get two boundary spins at a flat index sindex 
 ***************************************************************************************************/
void HyperCube::GetBoundaryPair(int sindex, int& bspin, int& tspin){
     bspin = sindex;
     tspin = N-S + sindex;
}

/****************************************************************************************************
 * Convert a flat index to a 3-dimensional coordinate
 ***************************************************************************************************/
void HyperCube::GetFullCrd(int sindex, int& x, int& y, int& z){ 
     z = (int) sindex / (X*Y);
     
     sindex -= z*(X*Y);
     y = (int) sindex / X;
     
     sindex -= y*X;
     x = sindex;
};

#endif
