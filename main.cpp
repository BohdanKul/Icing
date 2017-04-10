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
#include "lattice.h"
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
        ("seed, s", po::value<long>()->default_value(0),   "random generator seed")
        ("meas, m", po::value<int>(),                      "number of measurements (bins) to take")
        ("binsize", po::value<int>()->default_value(100),  "number of MC sweeps per bin")
        ;
    
    po::options_description PhysicalOptions("Physical options");
    PhysicalOptions.add_options()
        //("beta,   b",         po::value<double>()->default_value(beta_c),  
        //                                          "inverse temperature")
        ("T",                 po::value<double>(), "temperature")
        ("signJ,  J",         po::value<int>()->default_value(1),   
                                                  "ferromagnetic (+1) or antiferromagnetic (-1) coupling")
        ("width,  w",         po::value<int>(),   "lattice width")
        ("height, h",         po::value<int>(),   "lattice height")
        ("A, A",              po::value<int>()->default_value(0),
                                                   "region A widht")
        ("Ap, P",             po::value<int>()->default_value(0),
                                                   "region A prime width")
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

   
    //// beta
    double beta   = 0.0;
    //if (not(params.count("beta"))){
    //    cerr << "Error: define the inverse temperature (beta)" << endl;
    //    return 1;
    //}
    //else 
    //    beta = params["beta"].as<double>();

    // temperature 
    double T   = 0.0;
    if (not(params.count("T"))){
        cerr << "Error: define a temperature (T)" << endl;
        return 1;
    }
    else 
        T = params["T"].as<double>();
    beta = 1.0/T;
    

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
    int width;
    int height;
    if (not(params.count("width")) or not(params.count("height"))){
        cerr << "Error: define lattice dimensions (width, height)" << endl;
        return 1;
    }
    else{
        width  = params["width"].as<int>();
        height = params["height"].as<int>();
    }

    // Define physical constants --------------------------------------------------------
    int   N      = width * height;
    int   reps   = 1;

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

    boo::mt19937 engine(seed);
    boo::uniform_int<uint64_t> UDInt(0, N-1); 
    boo::uniform_real<double>  UDReal(0,1);
    boo::variate_generator<boo::mt19937&, boo::uniform_int<uint64_t> > RandInt( engine, UDInt );
    boo::variate_generator<boo::mt19937&, boo::uniform_real<double> >   RandReal(engine, UDReal);
    
    
    // Initialize a spins state ---------------------------------------------------------
    Spins spins;
    for (int i=0; i!=width*height; i++){
        spins.Add((RandInt()%2)*2-1);
    }
    //spins.print();

    
    // Initialize region A and A', dA = A - A' (A' is always assummed to be larger)
    int A_size   = params["A"].as<int>();
    int Ap_size = params["Ap"].as<int>();
    vector<int> A;
    vector<int> Ap;
    vector<int> dA;
    A.clear(); Ap.clear(); dA.clear();
    for (auto i=0;      i!=A_size;  i++) A.push_back(i);
    for (auto i=0;      i!=Ap_size; i++) Ap.push_back(i);
    for (auto i=A_size; i!=Ap_size; i++) dA.push_back(i);
   
     
    // Initialize the main simulation cell corresponding to the region A --------------------
    SimulationCell SC(width, height, &spins, A);
    //SC.print();
   
   
    // Initialize the cluster builders -------------------
    ClusterBuilder CB(SC, beta, signJ, RandReal); 
    

    
    // Initialize file objects ----------------------------------------------------------
    Communicator communicator(reps, width, height, A_size, Ap_size, T, seed);
    string eHeader = "";
    eHeader += boo::str(boo::format("#%15s%16s%16s%16s%16s")%"e"%"e2"%"m"%"m2"%"|m|"); 
    *communicator.stream("estimator")<< eHeader << endl; 
    long ID = communicator.getId();


    // Start the main MC loop -----------------------------------------------------------
    int initSpin;
    int spinA; 
    int EF;
    double debug = false;
    long  ET; long E2T;           long _ET;
    long  MT; long M2T; long aMT; long _MT;
    double ZR;
    double ZR2;
    double RN=0;
    for (auto i=0; i!=Nmeas; i++){
        ET = 0; E2T = 0; 
        MT = 0; M2T = 0; 
       aMT = 0; 
        ZR = 0; ZR2 = 0;  
        for (auto j=0; j!=binSize; j++){
            if (debug) cout << "--- " << i << " " << j << endl;
            
            if (debug) cout << "    cluster partition is reset" << endl;

            // perform a Wolff cluster update
            
            initSpin = RandInt() ; // pick randomly the initial spin
           
            //CB.ResetPartition();
            //CB.FlipTraceCluster(0, initSpin);
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
            for (auto k=0; k!=SC.GetSize(); k++){
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
            GetEnergyMagnetization(SC, _ET, _MT);
            ET  += _ET;
            E2T += _ET*_ET;

            MT  += _MT;
            M2T += _MT*_MT;
            aMT += labs(_MT);

            //cout << " Magnetization: " << _MT << " cumulative: " << MT << endl;
            //for (auto k=0; k!=SC.GetSize(); k++){
            //     cout << SC.GetSpins().Get(k) << " ";
            //}
            //cout << endl;


        }
        
        *communicator.stream("estimator") << boo::str(boo::format("%16.8E") %(signJ*ET / ( 1.0*binSize*N   )));
        *communicator.stream("estimator") << boo::str(boo::format("%16.8E") %(     E2T / ( 1.0*binSize*N*N )));
        *communicator.stream("estimator") << boo::str(boo::format("%16.8E") %(      MT / ( 1.0*binSize*N   )));
        *communicator.stream("estimator") << boo::str(boo::format("%16.8E") %(     M2T / ( 1.0*binSize*N*N )));
        *communicator.stream("estimator") << boo::str(boo::format("%16.8E") %(     aMT / ( 1.0*binSize*N   )));
        *communicator.stream("estimator") << endl;    
        
        cout << ID << ": Measurement taken" << endl;
    }
    return 0;
}
