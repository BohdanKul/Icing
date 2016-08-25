#ifndef SIMULATIONCELL_H
#define SIMULATIONCELL_H

//#include <vector>
#include <algorithm>  // find routine
#include <iostream>

#include "spins.h"
#include "lattice.h"

using namespace std;


/****************************************************************************************************
 * The main data structure containing all the information about a simulation. 
 ***************************************************************************************************/

class SimulationCell{
    protected:
        int width; int height; int N;    // geometrical quantities of the full system
        
        int r;                           // number of replicas
         
        int rep_width; int rep_height; int rep_N; // geometrical quantities of each replica

        Spins*  spins;                   // the spins state
        vector<pair<int, int>> boundary; // spins connectivity at the boundary between replicas
                                         // as dictated by the region A
        
        vector<vector<int>> nghbs;       // for each index, the vector keeps track of the neighbours

        int  GetFlatCrd(int _r, int _sindex){ // get the index of a spin located at _sindex of replica r
            return _r*rep_N + _sindex;         
        }

    public:
        SimulationCell(int _width, int _height, Spins* _spins, vector<int>& _A);
       
        Spins& GetSpins(){                     // get access to the spins  
            return *spins;
        };

        int GetSize(){                         // get the full size of simulation cell
            return N;
        };

        vector<vector<int>>& GetLattice(){     // get access to the neighbours info 
            return nghbs;
        };

        vector<pair<int,int>>& GetBoundary(){  // get the boundary connectivity
            return boundary;
        };

        // ugly but useful debugging routines
        void print(){
            cout << endl << "--- Simulation cell state ---" << endl;
            cout << endl << "    r: " << r << " width: " << width << " height: " << height << endl; 
            cout << endl << "   Spins state:" << endl << "   ";
            spins->print();
            
            cout << endl << "   Neighbours state: " << endl;
            for (int i=0; i!=N; i++){
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

};

/****************************************************************************************************
 * Initialize the spins state object and build a rectangular lattice connectivity data structures 
 * required to run cluster and normal spin updates on the simulation cell. The region A determines
 * the spins connectivity at the boundary between two replicas constituing the lattice.
 ***************************************************************************************************/
SimulationCell::SimulationCell(int _width, int _height, Spins* _spins, vector<int>& _A){
    
    spins = _spins;
    
    width  = _width;
    height = _height;
    N      = width*height;

    rep_width = width;
    rep_height = (int) height/2;
    rep_N = rep_width*rep_height;

    
    // initialize two replicas that are not connected to 
    // each other at first. Those are temporary data structures
    // helping to build the full lattice.
    vector<Rectangle> replicas;
    replicas.push_back(Rectangle(0, rep_width, rep_height, false)); 
    replicas.push_back(Rectangle(1, rep_width, rep_height, false)); 
    r=2;
   
    replicas[0].print();
    replicas[1].print();
   
    int bspin1; int tspin1;
    int bspin2; int tspin2;


    for (auto i=0; i!=width; i++){
        
        // determine the boundary spins that need to get rewired 
        bspin1 = replicas[0].GetBoundary().at(i).first;
        tspin1 = replicas[0].GetBoundary().at(i).second;
        cout << bspin1 << "," << tspin1 << endl; 

        bspin2 = replicas[1].GetBoundary().at(i).first;
        tspin2 = replicas[1].GetBoundary().at(i).second;
        //cout << bspin2 << "," << tspin2 << endl; 

        // if the spins aren't part of region A, link them intra-replicas
        if (find(_A.begin(), _A.end(), bspin1)==_A.end()) {
            replicas[0].GetUnitCell(bspin1).ResetNghb('b', Crd(0, tspin1));
            replicas[0].GetUnitCell(tspin1).ResetNghb('t', Crd(0, bspin1));
            replicas[1].GetUnitCell(bspin2).ResetNghb('b', Crd(1, tspin2));
            replicas[1].GetUnitCell(tspin2).ResetNghb('t', Crd(1, bspin2));

            boundary.push_back(make_pair(GetFlatCrd(0, bspin1), GetFlatCrd(0, tspin1)));
            boundary.push_back(make_pair(GetFlatCrd(1, bspin2), GetFlatCrd(1, tspin2)));
        }
        // otherwise, an inter-replicas links are created
        else{
            replicas[0].GetUnitCell(bspin1).ResetNghb('b', Crd(1, tspin2));
            replicas[0].GetUnitCell(tspin1).ResetNghb('t', Crd(1, bspin2));
            replicas[1].GetUnitCell(bspin2).ResetNghb('b', Crd(0, tspin1));
            replicas[1].GetUnitCell(tspin2).ResetNghb('t', Crd(0, bspin1));
         
            boundary.push_back(make_pair(GetFlatCrd(0, bspin1), GetFlatCrd(1, tspin2)));
            boundary.push_back(make_pair(GetFlatCrd(0, tspin1), GetFlatCrd(1, bspin2)));
        }
    }

    replicas[0].print();
    replicas[1].print();
    
    cout << endl << "--- Milestone 2: " << " about to build the boundary" << endl; 
    // store all the connectivity information into a flat vector
    for (auto _r=0; _r!=r; _r++){
        for (auto _s=0; _s!=rep_N; _s++){
            nghbs.push_back(vector<int>());
            for (auto nghb  = replicas[_r].GetUnitCell(_s).GetNghbs().begin(); 
                      nghb != replicas[_r].GetUnitCell(_s).GetNghbs().end(); 
                      nghb++)
                nghbs[GetFlatCrd(_r, _s)].push_back(GetFlatCrd(nghb->GetRep(), nghb->GetInd()));

        }
    }

    // At this point, the informaion contained in the replicas is redundant.
    // Destroy them.
    cout << endl << "--- Milestone 3: " << " done building the boundary" << endl; 
    replicas.clear();
}

/****************************************************************************************************
 *
 ***************************************************************************************************/


#endif
