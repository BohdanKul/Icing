#ifndef COMMUNICATOR_H
#define COMMUNICATOR_H


#include <vector>
#include <fstream>
#include <iostream>
#include <fstream>
#include <unordered_map>

using namespace std;

class Communicator
{
    public:
        Communicator(int _r, int _X, int _Y, int _Z, int _sizeA, int _sizeAp, float _beta, long _p);
        fstream* stream(string _fileName); 
        long     getId(){return id;};

        string dataName;
        string outDir;
        vector <string>  types;
        vector <fstream> files;   

    private:
        long id;
        long p;
        void GenerateId();
        unordered_map <string,fstream*> mFStreams;
}; 
#endif
