#ifndef SIMULATIONCELL_H
#define SIMULATIONCELL_H

//#include <vector>
#include <algorithm>  // find routine
#include <iostream>

#include "spins.h"
#include "hypercube.h"

using namespace std;


/****************************************************************************************************
 * The main data structure containing all the information about a simulation. 
 ***************************************************************************************************/

class SimulationCell{
    protected:
        int rN;                          // number of spins in a replica
        char u; char d;                  // orientation in the +1 dimension
        int dim;

        Spins*  spins;                   // the spins state
        vector<pair<int, int>> boundary; // spins connectivity at the boundary between replicas
                                         // as dictated by the region A
        vector<vector<int>> nghbs;       // for each index, the vector keeps track of the neighbours

    public:
        SimulationCell(int _X, int _Y, int _Z, Spins* _spins, vector<int>& _A);
       
        Spins& GetSpins(){ return *spins;};    // get access to the spins  
        int    GetSize(){  return 2*rN;};         // get the full size of simulation cell
        int    GetDim(){   return dim;};
        char   GetU(){     return u;};
        char   GetD(){     return d;};
        
        vector<vector<int>>&   GetLattice(){  return nghbs;};    // get access to the neighbours info 
        vector<pair<int,int>>& GetBoundary(){ return boundary;}; // get the boundary connectivity
        void print();    // ugly but useful debugging routines
};

/****************************************************************************************************
 * Initialize the spins state object and build a rectangular lattice connectivity data structures 
 * required to run cluster and normal spin updates on the simulation cell. The region A determines
 * the spins connectivity at the boundary between two replicas constituing the lattice.
 ***************************************************************************************************/
SimulationCell::SimulationCell(int _X, int _Y, int _Z, Spins* _spins, vector<int>& _A){
    spins = _spins;

    // Initialize two replicas that are not connected to each other at first. 
    // Those are temporary data structures helping to build the full lattice.
    HyperCube HC(_X, _Y, _Z);
    
    rN = HC.GetSize();   
    d  = HC.GetDown();
    u  = HC.GetUp();
    dim = HC.GetDim();    

    int bspin1; int tspin1;
    int bspin2; int tspin2;

    for (auto sindex=0; sindex!=rN; sindex++){
        nghbs.push_back(HC.GetSpinNghbs(sindex));
    }
    for (auto sindex=0; sindex!=rN; sindex++){
        nghbs.push_back(vector<int>());
        for (auto nghb=nghbs[sindex].begin(); nghb!=nghbs[sindex].end(); nghb++){
            nghbs[rN+sindex].push_back(*nghb + rN);
        }
    }
     
    for (auto i=0; i!=HC.GetBoundarySize(); i++){
        
        // determine the boundary spins that need to get rewired 
        HC.GetBoundaryPair(i, bspin1, tspin1); 
        bspin2 = bspin1 + rN;
        tspin2 = tspin1 + rN; 

        // if the spins aren't part of region A, link them intra-replicas
        if (find(_A.begin(), _A.end(), bspin1)==_A.end()) {
            nghbs[bspin1][DMap.at(d)] = tspin1;
            nghbs[tspin1][DMap.at(u)] = bspin1;
            nghbs[bspin2][DMap.at(d)] = tspin2;
            nghbs[tspin2][DMap.at(u)] = bspin2;
            
            boundary.push_back(make_pair(bspin1, tspin1)); // the order is important!
            boundary.push_back(make_pair(bspin2, tspin2)); // the order is important!
        }
        // otherwise, an inter-replicas links are created
        else{
            nghbs[bspin1][DMap.at(d)] = tspin2;
            nghbs[tspin1][DMap.at(u)] = bspin2;
            nghbs[bspin2][DMap.at(d)] = tspin1;
            nghbs[tspin2][DMap.at(u)] = bspin1;
            
            boundary.push_back(make_pair(bspin1, tspin2)); // the order is important!
            boundary.push_back(make_pair(bspin2, tspin1)); // the order is important!
        }
    }

}

/****************************************************************************************************
 * Print the main data structures of the class. Useful for debugging purposes. 
 ***************************************************************************************************/
void SimulationCell::print(){
    cout << endl << "--- Simulation cell state ---" << endl;
    cout << endl << "   Spins state:" << endl << "   ";
    spins->print();
    
    cout << endl << "   Neighbours state: " << endl;
    for (int i=0; i!=(2*rN); i++){
        cout << "   " << setfill(' ') << setw(2) << i << ": ";
        for (auto nghb=nghbs[i].begin(); nghb!=nghbs[i].end(); nghb++)
            cout << setfill(' ') << setw(2) << *nghb << " ";
        cout << endl;
    }
    cout << endl;
    
    cout << "   Boundary state: " << endl << "   ";
    for (auto b=boundary.begin(); b!=boundary.end(); b++)
        cout << "("<< b->first << ", " << b->second << ") ";
    cout << endl;
};


#endif
