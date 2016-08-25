#ifndef SPINS_H
#define SPINS_H

#include <vector>
#include <iostream>

using namespace std;

/****************************************************************************************************
 * A container class holding information about the spins state in a flat vector
 ***************************************************************************************************/
class Spins{
    protected:
        vector<int> state;
        int N;

    public:
        Spins(){ N=0; };               // no hustle constructor
        
        Spins(int _N){                 // constructor with a pre-set vector size
            N=_N; 
            state.resize(N, 0);
        };

        int Get(int index){            // get the state of a spin at index
            return state[index];
        };

        int GetSize(){ return N;};     // get the size (not that it matters..)
        
        void Set(int index, int spin){ // set the spin value at index
            state[index] = spin;
        };

        void Add(int spin){            // add a new spin
            state.push_back(spin);
            N += 1;
        };

        void Flip(int index){          // flip the spin at index
            state[index] = -state[index];
        };

        // ugly but useful debugging routines
        void print(){
            cout << "--- Spins state ---" << endl << "   ";
            for (int i=0; i!=N; i++)
                cout << Get(i) << " ";
            cout << endl;
        };

};

#endif
