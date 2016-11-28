#include <boost/format.hpp>
#include <boost/filesystem.hpp>

#include "communicator.h"
#include <iostream>
#include <fstream>

using namespace std;
using namespace boost;

/****************************************************************************************************
 ***************************************************************************************************/
Communicator::Communicator(string rfName, int _r, int _X, int _Y, int _Z, int _sizeA, int _sizeAp, float _beta, long _p)
{
    p = _p;
     
    // Build the filename 
    dataName = str(format("%02d-%03d-%03d-%03d-%03d-%03d-b%06.3f") %_r %_X %_Y %_Z %_sizeA %_sizeAp %_beta);
    dataName +=str(format("-p%05d") %_p);

    types  = vector<string> {"estimator", "state"};
    outDir = "OUTPUT"; 
    
    // Generate id
    if  (rfName == ""){
        GenerateId();
    }
    else{
        string filename = string(find(rfName.rbegin(), rfName.rend(), '/').base(), rfName.end());
        id = atol(filename.substr(filename.length()-13, 9).c_str());
    }

    string  filename;
    for (auto type=types.begin(); type!=types.end(); type++){
        
        filename = str(format("%s/%s-%s-%09d.dat") %outDir %*type %dataName %id);
        
        // when loading from a state, open the state-file as a read-only
        if ((rfName!="") and (*type=="state")){
            mFStreams[*type] = new fstream(filename, ios_base::in);
        }
        // otherwise, write by appending
        else{
            mFStreams[*type] = new fstream(filename, ios_base::out | ios_base::app);
        }

        if (!*mFStreams[*type]){
           cerr << "Unable to process file: " << filename << endl;
           exit(EXIT_FAILURE); 
        }
    }
    sfile      = str(format("%s/%s-%s-%09d.dat") %outDir %"state" %dataName %id);
    sfile_copy = sfile + 'a';  
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
fstream* Communicator::stream(string _filetype)
{
    return mFStreams[_filetype];

}

/****************************************************************************************************
 ***************************************************************************************************/
void Communicator::reset(string _filetype)
{
    mFStreams[_filetype]->close();
    filesystem::copy_file(sfile, sfile_copy, filesystem::copy_option::overwrite_if_exists);
    mFStreams[_filetype]->open(sfile, fstream::out|fstream::trunc);
}


