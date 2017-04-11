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
        
        vector<vector<int>> nghbs;       // for each index, the vector keeps track of the neighbours

        void Init(int _width, int _height, Spins* _spins);

    public:
        SimulationCell(int _width, int _height, Spins* _spins);
       
        Spins& GetSpins(){                     // get access to the spins  
            return *spins;
        };

        int GetSize(){                         // get the full size of simulation cell
            return N;
        };

        vector<vector<int>>& GetLattice(){     // get access to the neighbours info 
            return nghbs;
        };


        void print();    // ugly but useful debugging routines
};



/****************************************************************************************************
 *
 ***************************************************************************************************/
SimulationCell::SimulationCell(int _width, int _height, Spins* _spins){
    vector<int> A;
    A.clear();
    Init(_width, _height, _spins);
}

/****************************************************************************************************
 * Initialize the spins state object and build a rectangular lattice connectivity data structures 
 * required to run cluster and normal spin updates on the simulation cell. The region A determines
 * the spins connectivity at the boundary between two replicas constituing the lattice.
 ***************************************************************************************************/
void SimulationCell::Init(int _width, int _height, Spins* _spins){

    spins = _spins;
    
    width  = _width;
    height = _height;
    N      = width*height;

    // initialize two replicas that are not connected to 
    // each other at first. Those are temporary data structures
    // helping to build the full lattice.
    vector<Rectangle> replicas;
    replicas.push_back(Rectangle(0, width, height, false)); 
   
   
    // store all the connectivity information into a flat vector
    for (auto _s=0; _s!=N; _s++){
        nghbs.push_back(vector<int>());
        for (auto nghb  = replicas[0].GetUnitCell(_s).GetNghbs().begin(); 
                  nghb != replicas[0].GetUnitCell(_s).GetNghbs().end(); 
                  nghb++)
            nghbs[_s].push_back(nghb->GetInd());

    }

    // At this point, the informaion contained in the replicas is redundant.
    // Destroy them.
    //cout << endl << "--- Milestone 3: " << " done building the boundary" << endl; 
    replicas.clear();

}

/****************************************************************************************************
 * Print the main data structures of the class. Useful for debugging purposes. 
 ***************************************************************************************************/
void SimulationCell::print(){
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
    
};


#endif
