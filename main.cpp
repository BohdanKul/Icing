#include <vector>
#include <math.h>
#include <iostream>

#include <boost/random/uniform_real.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/variate_generator.hpp>

#include <boost/program_options.hpp>

#include <boost/format.hpp>
#include <boost/filesystem.hpp>

#include "spins.h"
#include "simulationcell.h"
#include "subroutines.h"
#include "communicator.h"
#include "clusterbuilder.h"

using namespace std;
namespace boo = boost;
namespace po  = boost::program_options; 

// --------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------
int main(int argc, char *argv[]){

    double beta_c = log(1.0+sqrt(2.0))/2.0;
    po::options_description SimulationOptions("Simulation options");
    SimulationOptions.add_options()
        ("help,h",  "produce help message")
        ("seed, s", po::value<long>()->default_value(0),    "random generator seed")
        ("meas, m", po::value<int>(),                       "number of measurements (bins) to take")
        ("binsize", po::value<int>()->default_value(100),   "number of MC sweeps per bin")
        ("state",   po::value<string>()->default_value(""), "path to the state file")
        ;
   
    bool snake_tiling; 
    po::options_description PhysicalOptions("Physical options");
    PhysicalOptions.add_options()
        ("beta,   b", po::value<double>(),  "inverse temperature")
        ("T,      T", po::value<double>(),  "temperature")
        ("signJ,  J", po::value<int>()->default_value(-1),         "ferromagnetic (-1) or antiferromagnetic (+1) coupling")
        ("X, X",      po::value<int>(),                            "lattice X")
        ("Y, Y",      po::value<int>(),                            "lattice Y")
        ("Z, Z",      po::value<int>(),                            "lattice Z") 
        ("snake",     po::bool_switch(&snake_tiling)->default_value(false),"snake-like build up of region A") 
        ("Ax",        po::value<int>(),                            "region A width")
        ("Ay",        po::value<int>()->default_value(1),          "region A height")
        ("APx",       po::value<int>(),                            "region A prime width")
        ("APy",       po::value<int>()->default_value(1),          "region A prime height")
        ;

    po::options_description cmdLineOptions("Command line options");
    cmdLineOptions.add(SimulationOptions).add(PhysicalOptions);

    po::variables_map params;
    po::store(po::parse_command_line(argc, argv, cmdLineOptions), params);
    po::notify(params);
    
    // help
    if (params.count("help")){
        cout << cmdLineOptions << "\n";
        return 1;
    }

    // seed
    int Nmeas;
    if (not(params.count("meas"))){
        cerr << "Error: define the number of bins to take (meas)" << endl;
        return 1;
    }
    else
        Nmeas = params["meas"].as<int>();

    // seed
    int binSize = params["binsize"].as<int>();

    // seed
    long seed = 1;
    if (params.count("seed"))
        seed = params["seed"].as<long>();

   
    // beta
    double beta   = 0.0;
    if (not(params.count("beta")) and not(params.count("T"))){
        cerr << "Error: define the temperature (beta or T)" << endl;
        return 1;
    }
    else{
        if (params.count("beta")) beta = params["beta"].as<double>();
        else                      beta = 1.0/params["T"].as<double>();
    }

    // sign J
    int signJ;
    if (not(params.count("signJ"))){
        cerr << "Error: define the nature of the coupling (signJ)" << endl;
        return 1;
    }
    else{
        if ((params["signJ"].as<int>()!=+1) and (params["signJ"].as<int>()!=-1)){
            cerr << "Error: signJ can't take other values but +/-1" << endl;
            return 1;
        }
        else{
            signJ = params["signJ"].as<int>();
        }
    }

    // lattice
    int X;
    int Y;
    int Z;

    if (not(params.count("X")) or not(params.count("Y")) or not(params.count("Z"))){
        cerr << "Error: define lattice dimensions (X, Y, Z)" << endl;
        return 1;
    }
    else{
        X = params["X"].as<int>();
        Y = params["Y"].as<int>();
        Z = params["Z"].as<int>();
    }

    // region A
    int Ax;
    int Ay;
    int APx;
    int APy;
    if (not(params.count("Ax")) or not(params.count("APx"))){
        cerr << "Error: define the entangled regions (Ax, Ay, APx, APy)" << endl;
        return 1;
    }
    else{ 
        Ax  = params["Ax"].as<int>();
        Ay  = params["Ay"].as<int>();
        APx = params["APx"].as<int>();
        APy = params["APy"].as<int>();
    }



    // Define physical constants --------------------------------------------------------
    int   N      = X * Y * Z;
    int   reps   = 2;

    // Initialize random objects --------------------------------------------------------
    vector<double> probTable = {1, exp(-4.0*beta), exp(-8.0*beta), exp(-12.0*beta)};

    boo::mt19937 engine(seed);
    boo::uniform_int<uint64_t> UDInt(0, reps*N-1); 
    boo::uniform_real<double>  UDReal(0,1);
    boo::variate_generator<boo::mt19937&, boo::uniform_int<uint64_t>> RandInt( engine, UDInt );
    boo::variate_generator<boo::mt19937&, boo::uniform_real<double> > RandReal(engine, UDReal);
   
    // Initialize file objects ----------------------------------------------------------
    int NA; int NAP;
    if (snake_tiling){ NA  = Ax;    NAP = APx; }
    else{              NA  = Ax*Ay; NAP = APx*APy; }

    Communicator communicator(params["state"].as<string>(), reps, X, Y, Z, NA, NAP, beta, seed);
    string eHeader = "";
    eHeader += boo::str(boo::format("#%15s%16s%16s")%"ET/N"%"Z_Ap/Z_A"%"Z_Ap/Z_A2"); 
    *communicator.stream("estimator")<< eHeader << endl; 
    long ID = communicator.getId();

 
    // The spins object
    Spins spins;
    
    // If restarting from a previously saved state, load the RN generator engine and spins states.
    if (params["state"].as<string>() != ""){
        // Buffer variables 
        string        sBuf;
        stringstream ssBuf;

        // Reset the distributions' state
         RandInt.distribution().reset();
        RandReal.distribution().reset();
    
        // Load the engine    
        getline(*communicator.stream("state"), sBuf);
        ssBuf << sBuf;
        ssBuf >> engine; 
        
        // Load the spins state
        istringstream issBuf;
        getline(*communicator.stream("state"), sBuf);
        issBuf.str(sBuf);
        int ispin;
        while (issBuf >> ispin){
            spins.Add(ispin);
        }
    }
    // Otherwise, initialize the spins state to random values
    else{
        for (int i=0; i!=N*reps; i++) 
            spins.Add((RandInt()%2)*2-1);
    }
    
    
    // Initialize region A and A', dA = A - A' (A' is always assummed to be larger)
    vector<int> A;  A.clear(); 
    vector<int> Ap; Ap.clear();
    vector<int> dA; dA.clear();
    
    // For the snake-like build up of region A, spins are added strip by strip until the number
    // specified by Ax or Apx is reached
    if (snake_tiling){
        for (auto x=0; x!=Ax;  x++) A.push_back(x);
        for (auto x=0; x!=APx; x++) Ap.push_back(x);
    }
    // Otherwise, the spins which lie within the square boundaries defined by (Ax, Ay) pair are added.
    else{
        for (auto y=0; y!=Ay; y++){
            for (auto x=0; x!=Ax; x++){
                A.push_back(y*X + x);
            }
        }
        for (auto y=0; y!=APy; y++){
            for (auto x=0; x!=APx; x++){
                Ap.push_back(y*X + x);
            }
        }
    }
    
    // Remember the difference between the two regions
    for (auto a=Ap.begin(); a!=Ap.end(); a++){
        if (find(A.begin(), A.end(), *a)==A.end()){
            dA.push_back(*a);
        }
    }

    // Initialize the main simulation cell corresponding to the region A --------------------
    SimulationCell SC(X, Y, Z, &spins, A);
    //SC.print();
   
    // Also initilize the simulation cell based on the region A'. It is used as a A'-BCs builder
    SimulationCell SCP(X, Y, Z, &spins, Ap);
    //SCP.print();
    
    // Initialize the cluster builders -------------------
    ClusterBuilder CB(SC, SCP,  beta, signJ, RandReal); 
    

    // Start the main MC loop -----------------------------------------------------------
    int initSpin;
    int spinA; 
    int EF;
    double debug = false;
    long  ET;
    double ZR;
    double ZR2;
    double RN=0;
    
    // Set up MC update parameters
    const float fw = 1.0/16.0;  // fraction controlling the scaling of Wolff clusters per a MC sweep
    int Nclusters = int(fw*Z);  // number of Wolff clusters per a sweep
    if (Nclusters < 1){ Nclusters = 1;};
    
    const float fs = 1.0/8.00;  // fraction controlling the number of single spin updates per a MCsweep
    int Nss  = int(fs*SC.GetSize()); // number of single spin updates 

    for (auto i=0; i!=Nmeas; i++){
        ET = 0; ZR = 0; ZR2 = 0;
        for (auto j=0; j!=binSize; j++){
            if (debug) cout << "--- " << i << " " << j << endl;
            
            if (debug) cout << "    cluster partition is reset" << endl;

            // perform Wolff cluster updates
            CB.ResetPartition();
            for (auto k=0; k!= Nclusters; k++){
                initSpin = RandInt() ; // pick randomly the initial spin
                CB.FlipTraceCluster(0, initSpin);
            }
            //CB.CrumbTraceCluster(0, initSpin);
            //int t=0;
            //for (auto s=CB.GetPartition().begin(); s!=CB.GetPartition().end(); s++){
            //    if (*s != -1)
            //        SC.GetSpins().Flip(t);
            //    t += 1;
            //}

            if (debug){
                cout << "    cluster is flipped " << endl; //<< Cluster << endl << "    ";
                SC.GetSpins().print();
            }
            // perform single spin updates
            int ispin = 0;
            for (auto k=0; k!=Nss; k++){
                ispin = RandInt();
                EF    = GetEffectiveField(SC, ispin);
                spinA = SC.GetSpins().Get(ispin);
                
                RN = RandReal();
                if  (spinA*signJ*EF>0)
                    SC.GetSpins().Flip(ispin);   
                else{
                    if (RandReal() < probTable[(int) abs(EF)/2]) 
                        SC.GetSpins().Flip(ispin);
                    }
            }

            // Perform measurements --------------------------------------------------------
            //ET += GetEnergy(SC);
            //ZR += exp(-beta*signJ*( GetBoundaryEnergy(SCP) - GetBoundaryEnergy(SC)));
            
            ZR += exp(-beta*signJ*( GetBoundaryEnergy(SCP, dA) - GetBoundaryEnergy(SC, dA)));
            
            int bspin1; int bspin2;
            int tspin1; int tspin2;
            int NCA = 0;
            CB.ResetPartition();
            CB.ResetLinks();
            for (auto k=dA.begin(); k!=dA.end(); k++){
                bspin1 = SC.GetBoundary()[2*(*k)  ].first; 
                tspin1 = SC.GetBoundary()[2*(*k)  ].second;
                bspin2 = SC.GetBoundary()[2*(*k)+1].first; 
                tspin2 = SC.GetBoundary()[2*(*k)+1].second;
                
                if (CB.CrumbTraceCluster(NCA, bspin1)==true) NCA +=1; 
                if (CB.CrumbTraceCluster(NCA, tspin1)==true) NCA +=1;
                if (CB.CrumbTraceCluster(NCA, bspin2)==true) NCA +=1; 
                if (CB.CrumbTraceCluster(NCA, tspin2)==true) NCA +=1;
            }
            CB.ReconnectLinks(dA);
            CB.ResetPartition();

            int NCP = 0;
            for (auto k=dA.begin(); k!=dA.end(); k++){
                bspin1 = SCP.GetBoundary()[2*(*k)  ].first; 
                tspin1 = SCP.GetBoundary()[2*(*k)  ].second;
                bspin2 = SCP.GetBoundary()[2*(*k)+1].first; 
                tspin2 = SCP.GetBoundary()[2*(*k)+1].second;
                
                if (CB.EatCrumbs(NCP, bspin1)==true) NCP +=1; 
                if (CB.EatCrumbs(NCP, tspin1)==true) NCP +=1;
                if (CB.EatCrumbs(NCP, bspin2)==true) NCP +=1; 
                if (CB.EatCrumbs(NCP, tspin2)==true) NCP +=1;
            }

            ZR2 += pow(2, NCP-NCA);
        }
        
        *communicator.stream("estimator") << boo::str(boo::format("%16.8E") %(signJ*ET/(1.0*binSize*reps*N)));
        *communicator.stream("estimator") << boo::str(boo::format("%16.8E") %(ZR /(1.0*binSize)));
        *communicator.stream("estimator") << boo::str(boo::format("%16.8E") %(ZR2/(1.0*binSize)));
        *communicator.stream("estimator") << endl;    
        
    
        // Save the state. 
        // Reseting the distributions state doesn't interfer with the ongoing generation of RNs.
        // However, it is required for a proper restart. 
         RandInt.distribution().reset();
        RandReal.distribution().reset();
        communicator.reset("state");
        *communicator.stream("state") << engine << endl;
        
        for (auto k=0; k!=SC.GetSize(); k++)
            *communicator.stream("state") << SC.GetSpins().Get(k) << " ";
        *communicator.stream("state") << endl;    
    
        for (auto k=0; k!=SC.GetSize(); k++)
            *communicator.stream("state") << SC.GetSpins().Get(k) << " ";
        *communicator.stream("state") << endl;    
    
        //for (auto l=0; l!=2; l++){
        //    for (auto k=0; k!=N; k++)
        //        *communicator.stream("spins") << (SC.GetSpins().Get(l*N + k)+1)/2 << " ";
        //    *communicator.stream("spins") << endl;
        //}

        cout << ID << ": Measurement taken" << endl;
    }
    return 0;
}
