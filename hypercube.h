#ifndef HYPERCUBE_H
#define HYPERCUBE_H

#include <iostream>
#include "unitcell.h"

using namespace std;

class HyperCube{
    private:
        vector<UnitCell> lattice;
        vector<pair<int, int>> boundary;  // a vector of spin pairs lying on the boundary in the z-direction
        int X; int Y; int Z; int N;
        int r;
        int dim;
        char up; char down;

        void GetFullCrd(int index, int& x, int& y, int& z){ // convert a 1-d index to a 3-d coordinate
            z = (int) index / (X*Y);
            
            index -= z*(X*Y);
            y = (int) index / X;
            
            index -= y*X;
            x = index;
        };
    public:
        HyperCube(int _r, int _X, int _Y, int _Z);
        char GetUp(){   return up;};
        char GetDown(){ return down;};
        int  GetSize(){ return N;};
        int  GetDim(){  return dim;};

        UnitCell& GetUnitCell(int sindex){ return lattice[sindex];}; // get a unit cell by its flat index
        vector<pair<int,int>>& GetBoundary(){ return boundary;};     // get the boundary

};

/****************************************************************************************************
 * A 3-d lattice class
 ***************************************************************************************************/
HyperCube::HyperCube(int _r, int _X, int _Y, int _Z = 1){
    X = _X;
    Y = _Y;
    Z = _Z;
    if (Z !=1 ){ dim = 3; up = 'n'; down = 'p';}
    else       { dim = 2; up = 'u'; down = 'd';}

    N = X*Y*Z;
    r = _r;

    int x; int y; int z; int siteB;
    for (auto sindex=0; sindex!=N; sindex++){
        GetFullCrd(sindex, x, y, z);

        // Construct a 3-dimensional unit cell
        UnitCell UC(dim);

        // add the right neighbour
        siteB = z*X*Y + y*X + (x+1)%X; 
        UC.SetNghb('r', Crd(r,siteB));

        // add the left neighbour
        if (x==0) siteB = z*X*Y + y*X + X-1;
        else      siteB = z*X*Y + y*X + x-1; 
        UC.SetNghb('l', Crd(r,siteB));

        // add the top neighbour
        siteB = z*X*Y + ((y+1)%Y)*X + x; 
        UC.SetNghb('u', Crd(r,siteB));

        // add the bottom neighbour
        if (y==0) siteB = z*X*Y + (Y-1)*X + x;
        else      siteB = z*X*Y + (y-1)*X + x;
        UC.SetNghb('d', Crd(r,siteB));
        
        if (dim == 3){
            // add a neighbour on the next depth layer
            siteB = ((z+1)%Z)*X*Y + y*X + x;
            UC.SetNghb('n', Crd(r,siteB));

            // add a neighbour on the previous depth layer
            if (z==0) siteB = (Z-1)*X*Y + y*X + x;
            else      siteB = (z-1)*X*Y + y*X + x;
            UC.SetNghb('p', Crd(r,siteB));
        }
        // store it in the vector
        lattice.push_back(UC);
    }

    int siteA;
    if (dim == 3){
        for (auto y=0; y!=Y; y++){
            for (auto x=0; x!=X; x++){
                siteA =     0*X*Y + y*X + x;
                siteB = (Z-1)*X*Y + y*X + x;
                boundary.push_back(make_pair(siteA, siteB));
            }
        }
    }
    else{
        for (auto x=0; x!=X; x++){
            siteA = x;
            siteB = X*(Y-1)+x;
            boundary.push_back(make_pair(siteA, siteB));
        }
    }
}
#endif
