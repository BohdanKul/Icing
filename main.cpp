#include <vector>
#include <math.h>
#include <iostream>

#include <boost/random/uniform_real.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/variate_generator.hpp>

#include <boost/format.hpp>
#include <boost/filesystem.hpp>

#include "spins.h"
#include "lattice.h"
#include "simulationcell.h"
#include "subroutines.h"
#include "communicator.h"
#include "clusterbuilder.h"

using namespace std;
using namespace boost;

// --------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------
int main(int argc, char *argv[]){
    
    // Define physical constants --------------------------------------------------------
    float beta   = 1.0;
    int   signJ  = 1;
    int   width  = 3;
    int   height = 6;
    int   N      = width * height;
    int   reps   = 2;

    // Initialize random objects --------------------------------------------------------
    vector<double> probTable = {1, exp(-4.0*beta), exp(-8.0*beta)};
    
    //double Sum = 0;
    //for (auto pe=probTable.begin(); pe!=probTable.end(); pe++)
    //    Sum += *pe;
    //
    //for (auto pe=probTable.begin(); pe!=probTable.end(); pe++){
    //    *pe = *pe/Sum;
    //    cout << *pe << endl;
    //}

    int seed = 1;
    boost::mt19937 engine(seed);
    boost::uniform_int<uint64_t> UDInt(0, N-1); 
    boost::uniform_real<float>   UDReal(0,1);
    boost::variate_generator<boost::mt19937&, boost::uniform_int<uint64_t> > RandInt( engine, UDInt );
    boost::variate_generator<boost::mt19937&, boost::uniform_real<float> >   RandReal(engine, UDReal);
    
    
    // Initialize a spins state ---------------------------------------------------------
    Spins spins;
    for (int i=0; i!=width*height; i++){
        spins.Add((RandInt()%2)*2-1);
    }
    spins.print();

    
    // Initialize region A --------------------------------------------------------------
    vector<int> A;
    A.clear();
    A.push_back(0);
    //A.push_back();
    //A.push_back();
    
    // Initialize the simulation cell ---------------------------------------------------
    SimulationCell SC(width, height, &spins, A);
    //SC.print();
    

    // Initialize helper objects required to calculate the partition function. Those are simulation 
    // cells corresponding to the region A' and emty and region A. 
    vector<int> Ap;
    Ap.clear();
    Ap.push_back(0);
    Ap.push_back(1);
    
    
    SimulationCell SCP(width, height, &spins, Ap);
    //SCP.print();
    
    SimulationCell SCF(width, height, &spins);
    //SCF.print();

    // Initialize the cluster builders -------------------
    ClusterBuilder CB( SC,  beta, signJ, RandReal); // one for a MC cluster update
    ClusterBuilder CBF(SCF, beta, signJ, RandReal); // and one for the EE estimator
    

    
    // Initialize file objects ----------------------------------------------------------
    Communicator communicator(reps, width, height, beta, seed);
    string eHeader = "";
    eHeader += boost::str(boost::format("#%15s%16s")%"ET/N"%"Z_Ap/Z_A"); 
    *communicator.stream("estimator")<< eHeader << endl; 
    long ID = communicator.getId();


    // Start the main MC loop -----------------------------------------------------------
    int initSpin;
    int spinA; 
    int EF;
    float debug = false;
    int Nmeas = 2;
    int binSize  = 1;
    long  ET;
    float ZR ;
    float RN=0;
    int cn=0; int CNA; int CNP;
    for (auto i=0; i!=Nmeas; i++){
        ET = 0; ZR = 0;
        for (auto j=0; j!=binSize; j++){
            if (debug) cout << "--- " << i << " " << j << endl;
            
            if (debug) cout << "    cluster partition is reset" << endl;

            // perform a Wolff cluster update
            
            initSpin = RandInt() ; // pick randomly the initial spin
            CB.ResetPartition();
            CB.TraceCluster(0, initSpin, true);

            if (debug){
                cout << "    cluster is flipped " << endl; //<< Cluster << endl << "    ";
                SC.GetSpins().print();
            }
            // perform single spin updates
            int ispin = 0;
            for (auto k=0; k!=SC.GetSize(); k++){
                ispin = RandInt();
                EF    = GetEffectiveField(SC, ispin);
                spinA = SC.GetSpins().Get(ispin);
                
                RN = RandReal();
                //cout << "Spins #" << ispin << " = " << spinA << " field: " << EF << " sign: " << signJ << " RN = " << RN << " index: " << (int) abs(EF)/2 << endl; 
                if  (spinA*signJ*EF>0)
                    SC.GetSpins().Flip(ispin);   
                else{
                    if (RandReal() < probTable[(int) abs(EF)/2]) 
                        SC.GetSpins().Flip(ispin);
                    }
            }
            //SC.GetSpins().print();

            // Perform measurements --------------------------------------------------------
            ET += GetEnergy(SC);
            
            cn = 0;
            CBF.ResetPartition();
            for (auto bpair = SCF.GetBoundary().begin();
                      bpair!= SCF.GetBoundary().end();
                      bpair++){
                if (CBF.TraceCluster(cn, bpair->first,  false)) cn += 1; 
                if (CBF.TraceCluster(cn, bpair->second, false)) cn += 1;
            }
            CNA = CBF.MergeClusters( SC.GetBoundary());
            CNP = CBF.MergeClusters(SCP.GetBoundary());
            ZR += pow(2, CNP-CNA);
        }
        
        *communicator.stream("estimator") << str(format("%16.8E") %(signJ*ET/(1.0*binSize*N)));
        *communicator.stream("estimator") << str(format("%16.8E") %(ZR/(1.0*binSize)));
        *communicator.stream("estimator") << endl;    
        
        //cout << ID << ": Measurement taken" << endl;
    }
    return 0;
}
