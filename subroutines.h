#include <math.h>

#ifndef SUBROUTINES_H
#define SUBROUTINES_H

#include "simulationcell.h"

using namespace std;


/****************************************************************************************************
* 
 ****************************************************************************************************/
int GetLocalEnergy(SimulationCell& SC, int ispin, double iangle){
    
    double E=0;
    for (auto jspin =SC.GetLattice().at(ispin).begin(); 
              jspin!=SC.GetLattice().at(ispin).end(); 
              jspin++)
        E += -1.0*cos(iangle - SC.GetSpins().Get(*jspin));
     
    return E;  
}


/****************************************************************************************************
* 
 ****************************************************************************************************/
int GetLocalEnergy(SimulationCell& SC, int ispin){
    
    double E = 0;
    double iangle = SC.GetSpins().Get(ispin);
    for (auto jspin =SC.GetLattice().at(ispin).begin(); 
              jspin!=SC.GetLattice().at(ispin).end(); 
              jspin++)
        E += -1.0*cos(iangle - SC.GetSpins().Get(*jspin));
     
    return E;  
}

/****************************************************************************************************
* Get the total energy 
 ****************************************************************************************************/
double GetEnergy(SimulationCell& SC){

    double energy = 0;

    int istate = 0;
    for (int ispin=0; ispin!=SC.GetSize(); ispin++){
        energy  += GetLocalEnergy(SC, ispin); 
    }

    return energy = energy/2.0; // divide by 2 to compensate for the double counting  
}

#endif
