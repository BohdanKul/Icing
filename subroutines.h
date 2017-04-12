#include <math.h>

#ifndef SUBROUTINES_H
#define SUBROUTINES_H

#include "simulationcell.h"

using namespace std;



/****************************************************************************************************
* 
 ****************************************************************************************************/
double GetLocalEnergy2(SimulationCell& SC, int ispin, double iangle){
    
    double E=0;
    for (auto jspin =SC.GetLattice().at(ispin).begin(); jspin!=SC.GetLattice().at(ispin).end(); jspin++){
        E += -1.0*cos(iangle - SC.GetSpins().Get(*jspin));
        //cout << E << " " << cos(iangle - SC.GetSpins().Get(*jspin)) << " "; 
    }
    //cout << endl;    
    return E;  
}


/****************************************************************************************************
* 
 ****************************************************************************************************/
double GetLocalEnergy1(SimulationCell& SC, int ispin){
    
    double E = 0;
    double iangle = SC.GetSpins().Get(ispin);
    for (auto jspin =SC.GetLattice().at(ispin).begin(); jspin!=SC.GetLattice().at(ispin).end(); jspin++){
        E += -1.0*cos(iangle - SC.GetSpins().Get(*jspin));
        //cout << E << " " << cos(iangle - SC.GetSpins().Get(*jspin)) << " "; 
    }
    //cout << endl;    
     
    return E;  
}

/****************************************************************************************************
* Get the total energy 
 ****************************************************************************************************/
double GetEnergy(SimulationCell& SC){

    double energy = 0;

    int istate = 0;
    for (int ispin=0; ispin!=SC.GetSize(); ispin++){
        energy  += GetLocalEnergy1(SC, ispin); 
        //cout << GetLocalEnergy1(SC, ispin) << " " << endl;
    }
    //cout << endl;
    return energy = energy/2.0; // divide by 2 to compensate for the double counting  
}

#endif
