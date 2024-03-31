/*
 * Filename:    mcspim.cpp
 * Description: execution entry for mcsp in-memory graph algorithms
 * Created:     04 March 2022
 * Authors:     Xinjie ZHOU
 */
#include "graph.hpp"

int main()
{
    /*unordered_map<int,unordered_multimap<NodeId,MCEdge>> HotPools;
    MCEdge e1(1, vector<EdgeWeight>(5,1));
    MCEdge e2(0, vector<EdgeWeight>(5,1));
    MCEdge e3(3, vector<EdgeWeight>(5,2));
    HotPools.insert({0,unordered_multimap<NodeId,MCEdge>()});
    HotPools[0].insert({0,e1});
    HotPools[0].insert({0,e3});
    HotPools[0].insert({1,e2});
    int count_i=-1;
    uint degree=0;

    for(auto it=HotPools[0].begin();it!=HotPools[0].end();it++){
        if(count_i <= 0) {
            if(count_i==0) cout<<endl;
            degree = HotPools[0].count(it->first);
            count_i = degree;
            cout<<it->first<<" "<< degree;
        }
        cout<<" "<<it->second.u;
        for(auto it2=it->second.w.begin();it2!=it->second.w.end();it2++){
            cout<<" "<<*it2;
        }
        count_i--;
    }
    cout<<endl;*/

    Timer tt;
    tt.start();
    gbstd::MCGraph mcg;
    string dt1 = "NY";//Beijing Louisiana NewYork California WUSA CTR USA Orkut UK2002
    mcg.dataset = dt1;
    cout << "Dataset: " << mcg.dataset << endl;
    cout << "Number of criteria: " << num_of_cri << endl;
    cout << "Query type: " << mcg.query_type << endl;
    cout << "Run times: " << run_times << endl;

//    mcg.MC_Evaluate(query_type);

    mcg.ReadGraph_W();
//    mcg.Dijkstra(705824, 263471);
    mcg.Dijkstra(17120, 223022);

    tt.stop();
    cout<<"\nTotal execution time: "<<tt.GetRuntime()<<" s."<<endl;
    return 0;
}
