#ifndef UNITCELL_H
#define UNITCELL_H

#include <map>
#include <iostream>
#include <iomanip>

using namespace std;
/****************************************************************************************************
 * Spin coordinate as define by its replica index r and its position within the replica ind
 ***************************************************************************************************/
class Crd{
    private:
        int r;   // replica index
        int ind; // index within the replica
    
    public:
        Crd(): r(-1), ind(-1){};

        Crd(int _r, int _ind){
            r   = _r;
            ind = _ind;
        };
   
        int GetRep(){  // get the replicas index
            return r;
        };
        
        int GetInd(){ // get the index with the replica
            return ind;
        };
};


// direction map - association between a relative direction and a number
const map<char, int> DMap = {
                                {'r', 0},
                                {'d', 1},
                                {'l', 2},
                                {'u', 3},
                            };

// the inverse map - from a number to a relative direction
const map<int, char> IDMap = {
                                {0, 'r'},
                                {1, 'd'},
                                {2, 'l'},
                                {3, 'u'}
                             };

// opposite direction map - associates with a each direction its opposite
const map<char, char> ODMap = {
                                {'l', 'r'},
                                {'r', 'l'},
                                {'u', 'd'},
                                {'d', 'u'}
                              };

// the same as the previous map but in terms of numbers
const map<int, int>  RDMap  = {
                                {0, DMap.at(ODMap.at(IDMap.at(0)))},
                                {1, DMap.at(ODMap.at(IDMap.at(1)))},
                                {2, DMap.at(ODMap.at(IDMap.at(2)))},
                                {3, DMap.at(ODMap.at(IDMap.at(3)))}
                              };

/****************************************************************************************************
 * A class containing information about the connectivity of a spin to its neighbours.
 ***************************************************************************************************/
class UnitCell{
    protected:
        vector<Crd>   nghbs; // all neighbours coordinates

    public:
        UnitCell(){
            nghbs.resize(4);
        }
        
        void SetNghb(char rel_pos, Crd _crd){  // add a new neigbour
            nghbs[DMap.at(rel_pos)] = _crd; 
        };
        
        vector<Crd>& GetNghbs(){              // return coordinates of all neighbours
            return nghbs;
        };

        // ugly but useful debugging routines
        void print(){
            int i=0;
            for (auto nghb=nghbs.begin(); nghb!=nghbs.end(); nghb++){
                cout  << IDMap.at(i) << " - ("<< setfill(' ') << setw(2) << nghb->GetRep() << ", "<< setfill(' ') << setw(2) << nghb->GetInd()  << ") ";
                i++;
            }
        };

};

#endif
