#ifndef LATTICE_H
#define LATTICE_H

#include <iostream>
#include <iomanip>

#include "unitcell.h"

using namespace std;


/****************************************************************************************************
 * A rectangular lattice of bonds 
 ***************************************************************************************************/
class Rectangle{
    
    private:
        vector<UnitCell> lattice;         // a lattice of connections
        vector<pair<int, int>> boundary;  // a vector of spin pairs lying on the y-boundary
        int  width; int  height; int N;   // geometric parameters
        int  r;                           // replica index the lattice belongs to

        bool PBC; 

        void GetFullCrd(int index, int& x, int& y){ // convert a 1-d index to a 2-d coordinate
            y = (int) index / width;
            x =       index % width;
        };

     public:
        Rectangle(int _r, int _width, int _height, bool _PBC);  // constructs a rectangular lattice


        UnitCell& GetUnitCell(int sindex){                      // get a unit cell by its flat index
            return lattice[sindex];
        };
        
        vector<pair<int,int>>& GetBoundary(){                   // get spins pairs lying on the y-boundary
            return boundary;
        };

        void print();
};


/****************************************************************************************************
 * A 2-d lattice class
 ***************************************************************************************************/
Rectangle::Rectangle(int _r, int _width, int _height, bool _PBC){
    
    // store lattice dimensions
    width  = _width;
    height = _height;
    N      = width * height;
    r      = _r;
    PBC    = _PBC; 

    // build a nearest neighbours 2-d lattice 
    int x; int y; int siteB;
    for (auto sindex=0; sindex!=N; sindex++){
        GetFullCrd(sindex, x, y);
        
        // add a new spin unit cell 
        UnitCell UC;
        
        // add the right neighbour
        siteB = y*width+(x+1)%width; 
        UC.SetNghb('r', Crd(r,siteB));

        // add the left neighbour
        if (x==0) siteB = y*width + width-1;
        else      siteB = y*width +     x-1; 
        UC.SetNghb('l', Crd(r,siteB));

        // add the top neighbour
        siteB = ((y+1)%height)*width + x; 
        UC.SetNghb('u', Crd(r,siteB));

        // add the bottom neighbour
        if (y==0) siteB = (height-1)*width + x;
        else      siteB = (y     -1)*width + x;
        UC.SetNghb('d', Crd(r,siteB));
        
        // store it in the vector
        lattice.push_back(UC);
    }

    // remember the boundary spins
    for (auto x=0; x!=width; x++)
        boundary.push_back(make_pair(x, width*(height-1)+x));
}

/****************************************************************************************************
 * Ugly but useful debugging routines
 ***************************************************************************************************/
void Rectangle::print(){
    cout << "--- Rectangle " << r << " state ---" << endl << "    ";
    cout << endl;
    cout << "    replica: " << r << " width: " << width << " height: " << height << endl << "   ";
    for (int i=0; i!=N; i++){
        cout << setfill(' ') << setw(2) << i << ": "; 
        GetUnitCell(i).print();
        cout << endl << "   ";
    }
    cout << endl;
    cout << "    Boundary: " << endl << "    ";
   
    for (auto b=boundary.begin(); b!=boundary.end(); b++)
        cout << "("<< b->first << ", " << b->second << ") ";
    cout << endl;

};

 
#endif
