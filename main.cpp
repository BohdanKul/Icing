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
#include "unitcell.h"

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

    //double Sum = 0;
    //for (auto pe=probTable.begin(); pe!=probTable.end(); pe++)
    //    Sum += *pe;
    //
    //for (auto pe=probTable.begin(); pe!=probTable.end(); pe++){
    //    *pe = *pe/Sum;
    //    cout << *pe << endl;
    //}

    const double PI = 3.141592653589793;    
    boo::mt19937 engine(seed);
    boo::uniform_int<uint64_t> UDInt(0, N-1); 
    boo::uniform_real<double>  UDReal(0,1);
    boo::uniform_real<double>  angleUDReal(-PI,PI);
    boo::variate_generator<boo::mt19937&, boo::uniform_int<uint64_t> >  RandInt( engine,  UDInt );
    boo::variate_generator<boo::mt19937&, boo::uniform_real<double> >   RandReal(engine,  UDReal);
    boo::variate_generator<boo::mt19937&, boo::uniform_real<double> >   RandAngle(engine, angleUDReal);
    
    
    // Initialize a spins state ---------------------------------------------------------
    Spins spins;
    for (int i=0; i!=width*height; i++){
        //spins.Add(RandAngle());
        spins.Add(0);
    }
    spins.Set(2, -3*PI);   
    
    //spins.print();

    
  
     
    // Initialize the main simulation cell corresponding to the region A --------------------
    SimulationCell SC(width, height, &spins);
    //SC.print();
   
   
    
    // Initialize file objects ----------------------------------------------------------
    Communicator communicator(width, height, T, seed);
    string eHeader = "";
    eHeader += boo::str(boo::format("#%15s%16s%16s%16s")%"e"%"e2"%"rho_x"%"rho_y"); 
    *communicator.stream("estimator")<< eHeader << endl; 
    long ID = communicator.getId();


    // Start the main MC loop -----------------------------------------------------------
    int initSpin;
    double debug = false;
    double  ET; double E2T; double _ET;
    double RN=0;
    double prop_angle; double E1; double E2; double rho_x; double rho_y; double _rho_x; double _rho_y;
    int ispin = 0; int next; int current; double angle1; double angle2;
    double sum1; double sum2;
    for (auto i=0; i!=Nmeas; i++){
        ET = 0; E2T = 0; rho_y = 0; rho_x = 0; 
        for (auto j=0; j!=binSize; j++){
            if (debug) cout << "--- " << i << " " << j << endl;
            
            if (debug) cout << "    cluster partition is reset" << endl;

            // perform a Wolff cluster update
            
           
            if (debug){
                cout << "    cluster is flipped " << endl; //<< Cluster << endl << "    ";
                SC.GetSpins().print();
            }
            
            // perform single spin updates
            for (auto k=0; k!=SC.GetSize(); k++){
                ispin = RandInt();
                E1    = GetLocalEnergy(SC, ispin);
                
                prop_angle = SC.GetSpins().Get(ispin) + RandAngle()*(1.0-tanh(.2/T));
                E2    = GetLocalEnergy(SC, ispin, prop_angle);

                RN = RandReal();
                if  (E2<E1){
                    SC.GetSpins().Set(ispin, prop_angle);   
                }
                else{
                    if (RandReal() < exp(-(E2-E1)/T)){ 
                        SC.GetSpins().Set(ispin, prop_angle);   
                    }
                }
            }

            // Perform measurements --------------------------------------------------------
            _rho_y = 0; 
            for (auto x=0; x!=width; x++){
               sum1 = 0; sum2 = 0;
               current = x;
               next    = SC.GetLattice().at(x).at(DMap.at('u'));
               while (next!=x){
                     angle1 = SC.GetSpins().Get(current);
                     angle2 = SC.GetSpins().Get(next);
                     
                     sum1 += cos(angle1 - angle2);
                     sum2 += sin(angle1 - angle2);
                    
                     current = next;
                     next    = SC.GetLattice().at(next).at(DMap.at('u'));
               }
               angle1 = SC.GetSpins().Get(current);
               angle2 = SC.GetSpins().Get(next);
               
               sum1 += cos(angle1 - angle2);
               sum2 += sin(angle1 - angle2);
                
               _rho_y += sum1 - sum2*sum2/T; 
            }
            rho_y += _rho_y;

            _rho_x = 0; 
            for (auto y=0; y!=height; y++){
               sum1 = 0; sum2 = 0;
               current = y*width;
               next    = SC.GetLattice().at(y*width).at(DMap.at('r'));
               while (next!=y*width){
                     angle1 = SC.GetSpins().Get(current);
                     angle2 = SC.GetSpins().Get(next);
                     
                     sum1 += cos(angle1 - angle2);
                     sum2 += sin(angle1 - angle2);
                    
                     current = next;
                     next    = SC.GetLattice().at(next).at(DMap.at('r'));
               }
               angle1 = SC.GetSpins().Get(current);
               angle2 = SC.GetSpins().Get(next);
               
               sum1 += cos(angle1 - angle2);
               sum2 += sin(angle1 - angle2);
                
               _rho_x += sum1 - sum2*sum2/T; 
            }
            rho_x += _rho_x;

            _ET = GetEnergy(SC);
            ET  += _ET;
            E2T += _ET*_ET;

        }
        
        *communicator.stream("estimator") << boo::str(boo::format("%16.8E") %(ET / ( 1.0*binSize*N   )));
        *communicator.stream("estimator") << boo::str(boo::format("%16.8E") %(E2T/ ( 1.0*binSize*N*N )));
        *communicator.stream("estimator") << boo::str(boo::format("%16.8E") %(rho_x / ( 1.0*binSize*N   )));
        *communicator.stream("estimator") << boo::str(boo::format("%16.8E") %(rho_y/ ( 1.0*binSize*N )));
        *communicator.stream("estimator") << endl;    
        
        cout << ID << ": Measurement taken" << endl;
    }
    return 0;
}
