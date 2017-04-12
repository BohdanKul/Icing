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
        vector<double> state;
        int N;

    public:
        Spins(){ N=0; };               // no hustle constructor
        
        Spins(int _N){                 // constructor with a pre-set vector size
            N=_N; 
            state.resize(N, 0);
        };

        double Get(int index){            // get the state of a spin at index
            return state[index];
        };

        int GetSize(){ return N;};     // get the size (not that it matters..)
        
        void Set(int index, double spin){ // set the spin value at index
            state[index] = spin;
        };

        void Add(double spin){            // add a new spin
            state.push_back(spin);
            N += 1;
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
