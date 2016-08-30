#include <boost/format.hpp>
#include <boost/filesystem.hpp>

#include "communicator.h"
#include <iostream>
#include <fstream>

using namespace std;
using namespace boost;

/****************************************************************************************************
 ***************************************************************************************************/
Communicator::Communicator(int _r, int _width, int _height, float _beta, long _p)
{
    p = _p;
     
    // Build the filename 
    dataName = str(format("%02d-%03d-%03d-b%06.3f") %_r %_width %_height %_beta);
    dataName +=str(format("-p%05d") %_p);

    types  = vector<string> {"estimator"};
    outDir = "OUTPUT"; 
    
    // Generate id
    GenerateId();
    string  fileName;
    for (auto type=types.begin(); type!=types.end(); type++){
        
        fileName = str(format("%s/%s-%s-%09d.dat") %outDir %*type %dataName %id);
        
        mFStreams[*type] = new fstream(fileName, ios_base::out | ios_base::app);
        
        if (!*mFStreams[*type]){
           cerr << "Unable to process file: " << fileName << endl;
           exit(EXIT_FAILURE); 
        }
    }
    
}

/****************************************************************************************************
 ***************************************************************************************************/
void Communicator::GenerateId()
{
    time_t seconds = long(time(NULL) - 39*365*24*60*60);
    id = long(seconds + p);
    string fName;
    fName = str(format("OUTPUT/estimator-%s-%09d.dat") % dataName %id);
    filesystem::path logPath(fName);

    while(filesystem::exists(logPath)) {
        id += 1;
        fName   = str(format("OUTPUT/estimator-%s-%09d.dat") % dataName %id);
        logPath = filesystem::path(fName);

    }

}

/****************************************************************************************************
 ***************************************************************************************************/
fstream* Communicator::stream(string _fileName)
{
    return mFStreams[_fileName];

}


