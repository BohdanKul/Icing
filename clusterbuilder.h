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
//#include "spins.h"
//#include "lattice.h"


using namespace std;
namespace boo=boost;

class ClusterBuilder{
    protected:
        double P;
        int signJ ;
        int nc;
        SimulationCell& SC;
        vector<int> Cluster_Partition;
        boost::variate_generator<boost::mt19937&, boost::uniform_real<double>>& RealRnd;
    public:
        ClusterBuilder(SimulationCell& _SC, double _beta, int _signJ, 
                       boost::variate_generator<boost::mt19937&, boost::uniform_real<double> >& _RealRnd);

        void ResetPartition();
        void TraceCluster(int nspin, bool toFlip);
        bool TraceCluster(int ncluster, int nspin, bool toFlip);
        int  MergeClusters(vector<pair<int, int>>& Boundary);
        int  MergeClusters2(vector<pair<int, int>>& Boundary);
        int  GetClustersN(){ return nc;};
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
    nc = 0;
}

/****************************************************************************************************
 
 
 ***************************************************************************************************/
void ClusterBuilder::TraceCluster(int nspin, bool toFlip){
     if (TraceCluster(nc, nspin, toFlip)) nc++;
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




/****************************************************************************************************
 * Given a vector of cluster indices associated with spins as well as a vector of spin pairs attached 
 * by the boundary, the function  counts how many clusters exist on the boundary for such a connection.
 ***************************************************************************************************/
int ClusterBuilder::MergeClusters(vector<pair<int, int>>& Boundary){
    int bclus; int uclus;

    list< set<int> > MergedClusters; 
    bool new_cluster_b;
    bool new_cluster_u;
    //set<int>* cluster_b;
    //set<int>* cluster_u;
    list<set<int>>::iterator cluster_b;
    list<set<int>>::iterator cluster_u;

    //cout << "   ----- Cluster Merge ----- " << endl << endl;
    //cout << "   Boundary: " << endl;
    //for (auto pspins=Boundary.begin(); pspins!=Boundary.end(); pspins++)
    //    cout << "   (" << pspins->first << "," << pspins->second << ")" << endl;
    //
    //cout << "   Cluster partition: " << endl;
    //for (auto cp=Cluster_Partition.begin(); cp!=Cluster_Partition.end(); cp++)
    //    cout << *cp << " ";
    //cout << endl;

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
    
    //int i = 0;
    //cout << "   Merged clusters: " << endl;
    //for (auto iset  = MergedClusters.begin();
    //          iset != MergedClusters.end();
    //          iset ++){
    //    cout << "   #" << i << ": ";
    //    for (auto ic  = iset->begin();
    //              ic != iset->end();
    //              ic ++){
    //        cout << *ic << " ";
    //    }
    //    cout << endl;
    //    i++;
    //}
    int NC =  MergedClusters.size();
    //cout << "   # clusters = " << NC << endl;

    return NC;
}


#endif
