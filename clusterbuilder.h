#ifndef CLUSTERBUILDER_H
#define CLUSTERBUILDER_H

#include <set>  // used in Merge Clusters routine

#include <boost/random/uniform_real.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/variate_generator.hpp>

#include "simulationcell.h"
//#include "spins.h"
//#include "lattice.h"


using namespace std;

class ClusterBuilder{
    protected:
        double P;
        int signJ ;
        SimulationCell& SC;
        vector<int> Cluster_Partition;
        boost::variate_generator<boost::mt19937&, boost::uniform_real<double>>& RealRnd;
    public:
        ClusterBuilder(SimulationCell& _SC, double _beta, int _signJ, 
                       boost::variate_generator<boost::mt19937&, boost::uniform_real<double> >& _RealRnd);

        void ResetPartition();
            
        bool TraceCluster(int ncluster, int nspin, bool toFlip);
        int  MergeClusters(vector<pair<int, int>>& Boundary);
};


/****************************************************************************************************
 ***************************************************************************************************/
ClusterBuilder::ClusterBuilder(SimulationCell& _SC, double _beta, int _signJ, 
                      boost::variate_generator<boost::mt19937&, boost::uniform_real<double> >& _RealRnd):
    RealRnd(_RealRnd),
    SC(_SC)
    {
    P       = 1.0 - exp(-2.0*_beta);
    signJ   = _signJ;

    Cluster_Partition.resize(SC.GetSize());
    ResetPartition();
}

/****************************************************************************************************
 ***************************************************************************************************/
void ClusterBuilder::ResetPartition(){
    // reset to default values 
    fill(Cluster_Partition.begin(), Cluster_Partition.end(), -1);
}

/****************************************************************************************************
 * A recursive function that traces a cluster for a given simulation cell. The result is stored in 
 * the Cluster_Partition vector which should be initiated to -1. There is an option to flip the cluster
 * controlled with the bool variable toFlip. The returned value indicates whether a new cluster has
 * been traced. 
 ***************************************************************************************************/
bool ClusterBuilder::TraceCluster(int ncluster, int nspin, bool toFlip){
    
    if (Cluster_Partition[nspin] != -1) return false;
    
    Cluster_Partition[nspin] = ncluster;
    
    if (toFlip)
        SC.GetSpins().Flip(nspin);
    
    //cout << " Checking neighbours" << endl;
    for (auto nghb  = SC.GetLattice().at(nspin).begin(); 
              nghb != SC.GetLattice().at(nspin).end(); 
              nghb++){
              //cout << "   " << nspin << " -> " << *nghb << endl;
              //cout << "   (" << SC.GetSpins().Get(nspin) << ") -> (" << SC.GetSpins().Get(*nghb) << ") sign J = " << signJ << endl;
        if  ((Cluster_Partition[*nghb] == -1) and (SC.GetSpins().Get(*nghb) * SC.GetSpins().Get(nspin))==signJ){
            if (RealRnd() < P)
                TraceCluster(ncluster, *nghb, toFlip); 
        }
    }

    return true;
}

/****************************************************************************************************
 * Given a vector of cluster indices associated with spins as well as a vector of spin pairs attached 
 * by the boundary, the function  counts how many clusters exist on the boundary for such a connection.
 ***************************************************************************************************/
int ClusterBuilder::MergeClusters(vector<pair<int, int>>& Boundary){
     int bclus; int uclus;

     vector< set<int> > MergedClusters; 
     bool new_merge;
     
     cout << "   ----- Cluster Merge ----- " << endl << endl;
     cout << "   Boundary: " << endl;
     for (auto pspins=Boundary.begin(); pspins!=Boundary.end(); pspins++)
         cout << "   (" << pspins->first << "," << pspins->second << ")" << endl;
     
     cout << "   Cluster partition: " << endl;
     for (auto cp=Cluster_Partition.begin(); cp!=Cluster_Partition.end(); cp++)
         cout << *cp << " ";
     cout << endl;

     for (auto Spin_Pair=Boundary.begin(); Spin_Pair!=Boundary.end(); Spin_Pair++){
         bclus = Cluster_Partition[ Spin_Pair->first  ];
         uclus = Cluster_Partition[ Spin_Pair->second ];
         
         if (bclus != uclus){
             new_merge = true;
             
             for (auto mcluster = MergedClusters.begin(); mcluster != MergedClusters.end(); mcluster++){
                 if (mcluster->find(bclus) != mcluster->end()){
                     mcluster->insert(uclus); 
                     new_merge = false;
                     break;
                 }
                 if (mcluster->find(uclus) != mcluster->end()){
                     mcluster->insert(bclus); 
                     new_merge = false;
                     break;
                 }
             }
             
             if (new_merge){
                set<int> tset;
                tset.insert(bclus);
                tset.insert(uclus);
                MergedClusters.push_back(tset); 
             }
        }
    }
    
    int NC =  MergedClusters.size();
    
    cout << "   # clusters = " << NC << endl;

    return NC;
}


#endif
