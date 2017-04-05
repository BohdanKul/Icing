#ifndef CLUSTERBUILDER_H
#define CLUSTERBUILDER_H

#include <set>  // used in Merge Clusters routine
#include <list>

#include <boost/random/uniform_real.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/pending/disjoint_sets.hpp>
#include <boost/property_map/property_map.hpp>

#include "simulationcell.h"
#include "unitcell.h"
//#include "spins.h"


using namespace std;
namespace boo=boost;

class ClusterBuilder{
    protected:
        double P;                       // derived probability of a link activation
        int signJ ;                     // the sign of the Hamiltonian (1 for the ferromagnetic) 
        SimulationCell& SC;             // used for the spins state, lattice and BCs
        vector<int> Cluster_Partition;  // employed during cluster building to mark spins according to
                                        // to their cluster 
        vector<int> Links_State;        // employed in "Crumb" methods. It keeps the state of a link:
                                        // 0 - default 
                                        // 1 - unactivated
                                        // 2 - activated
        boost::variate_generator<boost::mt19937&, boost::uniform_real<double>>& RealRnd;
    public:
        ClusterBuilder(SimulationCell& _SC, double _beta, int _signJ, 
                       boost::variate_generator<boost::mt19937&, boost::uniform_real<double> >& _RealRnd);

        vector<int>& GetPartition(){return Cluster_Partition;};
        
        void ResetLinks();     // reset to default values 0
        void ResetPartition(); // reset to default values -1
        
        bool  FlipTraceCluster(int ncluster, int nspin);
        
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
    Links_State.resize(SC.GetSize()*4);
    ResetPartition();
    ResetLinks();
}

/****************************************************************************************************
 ***************************************************************************************************/
void ClusterBuilder::ResetPartition(){
    // reset to default values 
    fill(Cluster_Partition.begin(), Cluster_Partition.end(), -1);
}

/****************************************************************************************************
 ***************************************************************************************************/
void ClusterBuilder::ResetLinks(){
    // reset to default values 
    fill(Links_State.begin(), Links_State.end(), 0);
}

/****************************************************************************************************
 * A recursive function that traces a cluster for a given simulation cell. The result is stored in 
 * the Cluster_Partition vector which should be initiated to -1. There is an option to flip the cluster
 * controlled with the bool variable toFlip. The returned value indicates whether a new cluster has
 * been traced. 
 ***************************************************************************************************/
bool ClusterBuilder::FlipTraceCluster(int ncluster, int nspin){
    if (Cluster_Partition[nspin] != -1) return false;
    
    Cluster_Partition[nspin] = ncluster;
    
    for (auto nghb  = SC.GetLattice().at(nspin).begin(); nghb != SC.GetLattice().at(nspin).end(); nghb++){
        if  ((Cluster_Partition[*nghb] == -1) and (SC.GetSpins().Get(*nghb) * SC.GetSpins().Get(nspin))==-signJ){
            if (RealRnd() < P){
                FlipTraceCluster(ncluster, *nghb); 
            }
        }
    }

    SC.GetSpins().Flip(nspin);
    return true;
}

#endif
