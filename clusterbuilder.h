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
#include "hypercube.h"
//#include "spins.h"


using namespace std;
namespace boo=boost;

class ClusterBuilder{
    protected:
        double P;                       // derived probability of a link activation
        int signJ ;                     // the sign of the Hamiltonian (1 for the ferromagnetic) 
        int dim;
        SimulationCell& SC;             // used for the spins state, lattice and BCs
        SimulationCell& SCP;            // used for the BCs
        vector<int> Cluster_Partition;  // employed during cluster building to mark spins according to
                                        // to their cluster 
        vector<int> Links_State;        // employed in "Crumb" methods. It keeps the state of a link:
                                        // 0 - default 
                                        // 1 - unactivated
                                        // 2 - activated
        boost::variate_generator<boost::mt19937&, boost::uniform_real<double>>& RealRnd;
    public:
        ClusterBuilder(SimulationCell& _SC, SimulationCell& _SCP, double _beta, int _signJ, 
                       boost::variate_generator<boost::mt19937&, boost::uniform_real<double> >& _RealRnd);

        vector<int>& GetPartition(){return Cluster_Partition;};
        
        void ResetLinks();     // reset to default values 0
        void ResetPartition(); // reset to default values -1
        
        void ReconnectLinks(vector<int>& dA); // 

        bool  FlipTraceCluster(int ncluster, int nspin);
        bool CrumbTraceCluster(int ncluster, int nspin);
        bool EatCrumbs(int ncluster, int nspin);

        int  MergeClusters(vector<pair<int, int>>& Boundary); // !!! no longer used !!!
        
};


/****************************************************************************************************
 ***************************************************************************************************/
ClusterBuilder::ClusterBuilder(SimulationCell& _SC, SimulationCell& _SCP, double _beta, int _signJ, 
                      boost::variate_generator<boost::mt19937&, boost::uniform_real<double> >& _RealRnd):
    RealRnd(_RealRnd),
    SC(_SC),
    SCP(_SCP)
    {
    P       = 1.0 - exp(-2.0*_beta);
    signJ   = _signJ;

    Cluster_Partition.resize(SC.GetSize());
    dim = SC.GetDim();
    Links_State.resize(SC.GetSize()*2*dim);
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
        if  ((Cluster_Partition[*nghb] != ncluster) and (SC.GetSpins().Get(*nghb) * SC.GetSpins().Get(nspin))==-signJ){
        //if  ((Cluster_Partition[*nghb] == -1) and (SC.GetSpins().Get(*nghb) * SC.GetSpins().Get(nspin))==-signJ){
            if (RealRnd() < P){
                FlipTraceCluster(ncluster, *nghb); 
            }
        }
    }

    SC.GetSpins().Flip(nspin);
    return true;
}

/****************************************************************************************************
 * Trace a cluster leaving a crumb trail of activated links. (actually the state of all links touching
 * the cluster is kept). Based on those links, the EatCrumbs method can later retrace the cluster.
 ***************************************************************************************************/
bool ClusterBuilder::CrumbTraceCluster(int ncluster, int nspin){
    if (Cluster_Partition[nspin] != -1){
        return false;
    }
    Cluster_Partition[nspin] = ncluster;
    int i = 0;
    double RN;
    for (auto nghb  = SC.GetLattice().at(nspin).begin(); nghb != SC.GetLattice().at(nspin).end(); nghb++){
        // first see if a link can be created
        if (SC.GetSpins().Get(*nghb)*SC.GetSpins().Get(nspin) == -signJ){
           // then check whether the spin on the other link's end has been marked 
            RN = RealRnd();
            if (Cluster_Partition[*nghb] == -1){
                if (RN < P){
                    Links_State[2*dim*nspin + i]           = 2;
                    Links_State[2*dim*(*nghb)+RDMap.at(i)] = 2; 
                    CrumbTraceCluster(ncluster, *nghb); 
                }
                else{
                    Links_State[2*dim*nspin + i]           = 1;
                    Links_State[2*dim*(*nghb)+RDMap.at(i)] = 1; 
                }
            }
            // activate the link (it is possible that it wasn't activated before) 
            else if ((Cluster_Partition[*nghb] == ncluster) and (Links_State[2*dim*nspin+i] == 0)){ 
                if (RN < P){
                    Links_State[2*dim*nspin + i]           = 2;
                    Links_State[2*dim*(*nghb)+RDMap.at(i)] = 2; 
                }
                else{
                    Links_State[2*dim*nspin + i]           = 1;
                    Links_State[2*dim*(*nghb)+RDMap.at(i)] = 1; 
                }
            }
        }
        else{
             Links_State[2*dim*nspin + i]           = 1;
             Links_State[2*dim*(*nghb)+RDMap.at(i)] = 1; 
        }
        i += 1;
    }
    
    return true;
}




/****************************************************************************************************
 * Retrace clusters based on the links activated by CrumbTraceClusters method. It is designed to be
 * employed after applying the ReConnectLinks method. 
 ***************************************************************************************************/
bool ClusterBuilder::EatCrumbs(int ncluster, int nspin){
    if (Cluster_Partition[nspin] != -1){
        return false;
    }

    Cluster_Partition[nspin] = ncluster;
    int i = 0;
    for (auto nghb  = SCP.GetLattice().at(nspin).begin(); nghb != SCP.GetLattice().at(nspin).end(); nghb++){
        if  ((Cluster_Partition[*nghb] == -1) and (Links_State[2*dim*nspin+i]==2)){
                EatCrumbs(ncluster, *nghb); 
        }
        i += 1;
    }
    
    return true;
}

/****************************************************************************************************
 * Once a configuration state has been partitionned in the space of active and passive links with the 
 * help of CrumbTraceCluster method, producing Links_State datastructure, the following method is
 * employed to reconfigure the links in correspondance with the region A' (A'-A = dA).
 ***************************************************************************************************/
void ClusterBuilder::ReconnectLinks(vector<int>& dA){
    int bspin1; int bspin2;
    int tspin1; int tspin2;
    // spins in the region dA are assumed to be connected intra-replicas.
    for (auto i=dA.begin(); i!=dA.end(); i++){
        bspin1 = SC.GetBoundary()[2*(*i)  ].first; 
        tspin1 = SC.GetBoundary()[2*(*i)  ].second;
        bspin2 = SC.GetBoundary()[2*(*i)+1].first; 
        tspin2 = SC.GetBoundary()[2*(*i)+1].second;
        
        // Set the activation of bottom spins in the downward direction.
        // This corresponds to keeping the upward links on the top spins. 
        Links_State[2*dim*bspin1+DMap.at(SC.GetD())] = Links_State[2*dim*tspin2+DMap.at(SC.GetU())];
        Links_State[2*dim*bspin2+DMap.at(SC.GetD())] = Links_State[2*dim*tspin1+DMap.at(SC.GetU())];
        
        // Alternatively, set the activation in the reverse order. Either choice works.
        // Each one corresponds to setting the direction of the imaginary time axis.
        //Links_State[2*dim*tspin1+DMap.at('u')] = Links_State[2*dim*bspin2+DMap.at('d')];
        //Links_State[2*dim*tspin2+DMap.at('u')] = Links_State[2*dim*bspin1+DMap.at('d')];
    }

}

/****************************************************************************************************
 * Given a vector of cluster indices associated with spins as well as a vector of spin pairs attached 
 * by the boundary, the function  counts how many clusters exist on the boundary for such a connection.
 ***************************************************************************************************/
int ClusterBuilder::MergeClusters(vector<pair<int, int>>& Boundary){
    int bclus; int uclus;

    list< set<int> > MergedClusters; 
    bool new_cluster_b;
    bool new_cluster_u;
    list<set<int>>::iterator cluster_b;
    list<set<int>>::iterator cluster_u;

    for (auto Spin_Pair=Boundary.begin(); Spin_Pair!=Boundary.end(); Spin_Pair++){
        bclus = Cluster_Partition[ Spin_Pair->first  ];
        uclus = Cluster_Partition[ Spin_Pair->second ];
        
        if (bclus == uclus){
           new_cluster_b = true;

           for (auto mcluster = MergedClusters.begin(); mcluster != MergedClusters.end(); mcluster++){
               if (mcluster->find(bclus) != mcluster->end()){
                   new_cluster_b = false;
                   break;
               }
           } 
        
           if (new_cluster_b){
              set<int> tset;
              tset.insert(bclus);
              MergedClusters.push_back(tset); 
           }
        }
        else{
            new_cluster_b = true;
            new_cluster_u = true;

            for (auto mcluster = MergedClusters.begin(); mcluster != MergedClusters.end(); mcluster++){
                if (mcluster->find(bclus) != mcluster->end()){ new_cluster_b = false; cluster_b = mcluster;}
                if (mcluster->find(uclus) != mcluster->end()){ new_cluster_u = false; cluster_u = mcluster;}
               
                if ( (not(new_cluster_b)) and (not(new_cluster_u)) ) break;

             }
        
            if       ((new_cluster_b) and not(new_cluster_u)) cluster_u->insert(bclus);
            else if  ((new_cluster_u) and not(new_cluster_b)) cluster_b->insert(uclus);
            else if  (new_cluster_u and new_cluster_b){
                     set<int> tset;
                     tset.insert(bclus);
                     tset.insert(uclus);
                     MergedClusters.push_back(tset); 
            }
            else if  (cluster_b != cluster_u){
                     if (cluster_u->size() > cluster_b->size()){
                         cluster_u->insert(cluster_b->begin(), cluster_b->end());
                         //list<set<int>>::iterator it(cluster_b);
                         MergedClusters.erase(cluster_b);
                     }
                     else{
                         cluster_b->insert(cluster_u->begin(), cluster_u->end());
                         //list<set<int>>::iterator it(cluster_u);
                         MergedClusters.erase(cluster_u);
                     }

            }

        }
    }
    
    int NC =  MergedClusters.size();

    return NC;
}
    
/****************************************************************************************************
 * ----------------------------------------!!!BROKEN!!!----------------------------------------------
 *  ---------------------------------- use on your own risk -----------------------------------------
 * Given a vector of cluster indices associated with spins as well as a vector of spin pairs attached 
 * by the boundary, the function  counts how many clusters exist on the boundary for such a connection.
 ***************************************************************************************************/
//int ClusterBuilder::MergeClusters(vector<pair<int, int> > &Boundary){
//    typedef map<int, size_t> rank_t;
//    typedef map<int, int>  parent_t;
//    rank_t   rank_map;
//    parent_t parent_map;
//
//    boo::associative_property_map<   rank_t>   rank_pmap(  rank_map);
//    boo::associative_property_map< parent_t> parent_pmap(parent_map);
//    boo::disjoint_sets< boo::associative_property_map<rank_t>,
//                        boo::associative_property_map<parent_t> > ds(rank_pmap, parent_pmap);
//
//     cout << "   ----- Cluster Merge ----- " << endl << endl;
//     cout << "   Boundary: " << endl;
//     for (auto pspins=Boundary.begin(); pspins!=Boundary.end(); pspins++)
//         cout << "   (" << pspins->first << "," << pspins->second << ")" << endl;
//     
//     cout << "   Cluster partition: " << endl;
//     for (auto cp=Cluster_Partition.begin(); cp!=Cluster_Partition.end(); cp++)
//         cout << *cp << " ";
//     cout << endl;
//    
//     
//     vector<int> dv; 
//    for (auto i=0; i!=GetClustersN(); i++){
//        ds.make_set(i);
//        dv.push_back(i);
//    }
//    cout << "Before: " << ds.count_sets(dv.begin(), dv.end()) << endl;
//
//    for (auto sp=Boundary.begin(); sp!=Boundary.end(); sp++){
//        ds.union_set(sp->first, sp->second); 
//    }
//    cout << "After:  " << ds.count_sets(dv.begin(), dv.end()) << endl;
//
//    return ds.count_sets(dv.begin(), dv.end());
//}

#endif
