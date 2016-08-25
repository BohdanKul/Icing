#ifndef LATTICE_H
#define LATTICE_H

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
        //Crd(int _r, int _ind): r(_r), ind(_ind){}; // init
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

/****************************************************************************************************
 * A class containing information about the connectivity of a spin to its neighbours.
 ***************************************************************************************************/
class UnitCell{
    protected:
        map<char, Crd> nghbs; // a map defining the relative neighbours position as up, down, etc
        vector<Crd>   vnghbs; // all neighbours coordinates

    public:
        void AddNghb(char rel_pos, Crd _crd){     // add a new neigbour
            //nghbs[rel_pos] = _crd;  
            nghbs.insert(make_pair(rel_pos, _crd));
            vnghbs.push_back(_crd); 
        };
        
        void DelNghb(char rel_pos){              // delete a neigbhour
            // check if the neigbour defined by 
            // its relative position already exists
            if (nghbs.find(rel_pos)!=nghbs.end()){

               // delete it in the map
               nghbs.erase(rel_pos);

               // since the position of this neigbour in the vnghbs 
               // is not know, the full vector must be reset
               vnghbs.clear();
               for (auto nghb=nghbs.begin(); nghb!=nghbs.end(); nghb++)
                    vnghbs.push_back(nghb->second);
            }
        };

        void ResetNghb(char rel_pos, Crd crd){ // reset the neigbour coordinates 
            DelNghb(rel_pos);
            AddNghb(rel_pos, crd);
        };

        vector<Crd>& GetNghbs(){              // return coordinates of all neighbours
            return vnghbs;
        };

        
        // ugly but useful debugging routines
        void print(){
            for (auto nghb=nghbs.begin(); nghb!=nghbs.end(); nghb++)
                cout  << nghb->first << " - ("<< setfill(' ') << setw(2) << nghb->second.GetRep() << ", "<< setfill(' ') << setw(2) << nghb->second.GetInd()  << ") ";
        };

};


/****************************************************************************************************
 * A rectangular lattice of bonds 
 ***************************************************************************************************/
class Rectangle{
    
    private:
        vector<UnitCell> lattice;         // a lattice of connections
        vector<pair<int, int>> boundary;  // a vector of spin pairs lying on the y-boundary
        int  width; int  height; int N;   // geometric parameters
        int  r;                           // replica index the lattice belongs to

        bool PBC; 

        void GetFullCrd(int index, int& x, int& y){ // convert a 1-d index to a 2-d coordinate
            y = (int) index / width;
            x =       index % width;
        };

     public:
        Rectangle(int _r, int _width, int _height, bool _PBC);  // constructs a rectangular lattice


        UnitCell& GetUnitCell(int sindex){                      // get a unit cell by its flat index
            return lattice[sindex];
        };
        
        vector<pair<int,int>>& GetBoundary(){                   // get spins pairs lying on the y-boundary
            return boundary;
        };

        // ugly but useful debugging routines
        void print(){
            cout << "--- Rectangle " << r << " state ---" << endl << "    ";
            cout << endl;
            cout << "    replica: " << r << " width: " << width << " height: " << height << endl << "   ";
            for (int i=0; i!=N; i++){
                cout << setfill(' ') << setw(2) << i << ": "; 
                GetUnitCell(i).print();
                cout << endl << "   ";
            }
            cout << endl;
            cout << "    Boundary: " << endl << "    ";
           
            for (auto b=boundary.begin(); b!=boundary.end(); b++)
                cout << "("<< b->first << ", " << b->second << ") ";
            cout << endl;
 
        };

 
};


/****************************************************************************************************
 * A 2-d lattice class
 ***************************************************************************************************/
Rectangle::Rectangle(int _r, int _width, int _height, bool _PBC){
    
    // store lattice dimensions
    width  = _width;
    height = _height;
    N      = width * height;
    r      = _r;
    PBC    = _PBC; 

    // build a nearest neighbours 2-d lattice 
    int x; int y; int siteB;
    for (auto sindex=0; sindex!=N; sindex++){
        GetFullCrd(sindex, x, y);
        
        // add a new spin unit cell 
        UnitCell UC;
        
        // add the right neighbour
        siteB = y*width+(x+1)%width; 
        UC.AddNghb('r', Crd(r,siteB));

        // add the left neighbour
        if (x==0) siteB = y*width + width-1;
        else      siteB = y*width +     x-1; 
        UC.AddNghb('l', Crd(r,siteB));

        // add the top neighbour
        siteB = ((y+1)%height)*width + x; 
        UC.AddNghb('t', Crd(r,siteB));

        // add the bottom neighbour
        if (y==0) siteB = (height-1)*width + x;
        else      siteB = (y     -1)*width + x;
        UC.AddNghb('b', Crd(r,siteB));
        
        // store it in the vector
        lattice.push_back(UC);
    }

    // remember the boundary spins
    for (auto x=0; x!=width; x++)
        boundary.push_back(make_pair(x, width*(height-1)+x));
}


#endif
