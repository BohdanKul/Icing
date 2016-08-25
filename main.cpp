#include <vector>
#include <math.h>
#include <iostream>

#include "spins.h"
#include "simulationcell.h"
#include "lattice.h"
#include "subroutines.h"

using namespace std;

/****************************************************************************************************
*
 ****************************************************************************************************/
int main(int argc, char *argv[]){
    cout << "   Milestone 0: " <<  endl;
    
    float RN   = 0.2;
    float beta = 1.0;
    int  signJ = 1;

    int Nmeas = 100;
    int Nbin  = 100;

    /****************************************************************************************************
    *
     ****************************************************************************************************/
    vector<double> probTable = {0.5, exp(-4.0*beta), exp(-8.0*beta)};



    int width  = 3;
    int height = 6;
    //int N      = width*height;

    //Spins * spins;
    //spins = new Spins();
    Spins spins;
    for (int i=0; i!=width*height; i++){
        spins.Add(1);
    }
    
    spins.print();

    cout << endl << "   Milestone 1: " << " --spins initialized -- " << endl;
    vector<int> A;
    A.clear();
    SimulationCell SC(width, height, &spins, A);
    SC.print();

    vector<int> Cluster_Partition;
    Cluster_Partition.resize(SC.GetSize());


    int initSpin;
    int spinA; 
    int EF;
    cout << endl << "--- Milestone 4: " << " start MC " << endl;
    for (auto i=0; i!=Nmeas; i++){
        for (auto j=0; j!=Nbin; j++){
            cout << "--- " << i << " " << j << endl;
            initSpin = 0;
            
            // reset to default values 
            fill(Cluster_Partition.begin(), Cluster_Partition.end(), -1);
            
            cout << "    cluster partition is reset" << endl;

            // perform a Wolff cluster update
            TraceCluster(SC, Cluster_Partition, 0, initSpin, true);

            cout << "    cluster is flipped" << endl << "    ";
            for (auto c=Cluster_Partition.begin(); c!=Cluster_Partition.end(); c++)
                cout << *c << " ";
            cout << endl;
            
            SC.GetSpins().print();

            // perform single spin updates
            for (auto ispin=0; ispin!=SC.GetSize(); ispin++){
                EF    = GetEffectiveField(SC, ispin);
                spinA = SC.GetSpins().Get(ispin);

                if  (spinA*signJ*EF)
                    SC.GetSpins().Flip(ispin);   
                else{
                    RN = 0;
                    if (RN < probTable[(int) abs(EF)/2]) 
                        SC.GetSpins().Flip(ispin);
                    }
            }
        }
    }
    return 0;
}
