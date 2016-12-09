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
        int r;                           // number of replicas
        int rN;                          // number of spins in a replica
        char u; char d;                  // orientation in the +1 dimension
        int dim;

        Spins*  spins;                   // the spins state
        vector<pair<int, int>> boundary; // spins connectivity at the boundary between replicas
                                         // as dictated by the region A
        vector<vector<int>> nghbs;       // for each index, the vector keeps track of the neighbours

        int  GetFlatCrd(int _r, int _sindex){ // get the index of a spin located at _sindex of replica r
            return _r*rN + _sindex;         
        }

    public:
        SimulationCell(int _X, int _Y, int _Z, Spins* _spins, vector<int>& _A);
       
        Spins& GetSpins(){ return *spins;};    // get access to the spins  
        int    GetSize(){  return r*rN;};         // get the full size of simulation cell
        int    GetDim(){   return dim;};
        char   GetU(){     return u;};
        char   GetD(){     return d;};
        
        vector<vector<int>>& GetLattice(){    return nghbs;};    // get access to the neighbours info 
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
    r=2;
    vector<HyperCube> replicas;
    replicas.push_back(HyperCube(0, _X, _Y, _Z)); 
    replicas.push_back(HyperCube(1, _X, _Y, _Z)); 
    
    rN  = replicas[0].GetSize();   
    dim = replicas[0].GetDim(); 
    
    u  = replicas[0].GetUp();
    d  = replicas[1].GetDown();
    

    int bspin1; int tspin1;
    int bspin2; int tspin2;


    for (auto i=0; i!=replicas[0].GetBoundary().size(); i++){
        
        // determine the boundary spins that need to get rewired 
        bspin1 = replicas[0].GetBoundary().at(i).first; 
        tspin1 = replicas[0].GetBoundary().at(i).second;
        bspin2 = replicas[1].GetBoundary().at(i).first; 
        tspin2 = replicas[1].GetBoundary().at(i).second;

        // if the spins aren't part of region A, link them intra-replicas
        if (find(_A.begin(), _A.end(), bspin1)==_A.end()) {
            replicas[0].GetUnitCell(bspin1).SetNghb(d, Crd(0, tspin1));
            replicas[0].GetUnitCell(tspin1).SetNghb(u, Crd(0, bspin1));
            replicas[1].GetUnitCell(bspin2).SetNghb(d, Crd(1, tspin2));
            replicas[1].GetUnitCell(tspin2).SetNghb(u, Crd(1, bspin2));

            boundary.push_back(make_pair(GetFlatCrd(0, bspin1), GetFlatCrd(0, tspin1)));
            boundary.push_back(make_pair(GetFlatCrd(1, bspin2), GetFlatCrd(1, tspin2)));
            
        }
        // otherwise, an inter-replicas links are created
        else{
            replicas[0].GetUnitCell(bspin1).SetNghb(d, Crd(1, tspin2));
            replicas[0].GetUnitCell(tspin1).SetNghb(u, Crd(1, bspin2));
            replicas[1].GetUnitCell(bspin2).SetNghb(d, Crd(0, tspin1));
            replicas[1].GetUnitCell(tspin2).SetNghb(u, Crd(0, bspin1));
         
            boundary.push_back(make_pair(GetFlatCrd(0, bspin1), GetFlatCrd(1, tspin2)));
            boundary.push_back(make_pair(GetFlatCrd(1, bspin2), GetFlatCrd(0, tspin1)));
        }
    }

    
    // store all the connectivity information into a flat vector
    for (auto _r=0; _r!=r; _r++){
        for (auto _s=0; _s!=replicas[_r].GetSize(); _s++){
            nghbs.push_back(vector<int>());
            for (auto nghb  = replicas[_r].GetUnitCell(_s).GetNghbs().begin(); 
                      nghb != replicas[_r].GetUnitCell(_s).GetNghbs().end(); 
                      nghb++)
                nghbs[GetFlatCrd(_r, _s)].push_back(GetFlatCrd(nghb->GetRep(), nghb->GetInd()));

        }
    }

    // At this point, the informaion contained in the replicas is redundant.  Destroy them.
    replicas.clear();

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
