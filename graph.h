/*
 * Filename:    graph.h
 * Created:     12 Feb 2022
 * Authors:     Xinjie ZHOU
 */

#ifndef GRAPH_H
#define GRAPH_H

//#include <boost/heap/binomial_heap.hpp>
//#include <boost/heap/priority_queue.hpp>
//#include <boost/heap/fibonacci_heap.hpp>
#include "head.h"

namespace gbstd{
    using namespace std;
    //Class of multi-criteria graph
    class MCGraph{
    public:
        /*** In-memory graph information ***/
        uint node_num = 0;      //the number of vertices
        unsigned long long edge_num = 0;      //the number of edges
        usint sc_i = 0;         //the index of criterion in processing
        string dataset;         //dataset name
        string criteria;            //the criteria(type) name of graph
        vector<string> mc_criteria; //the criteria(type) of graph
        vector<Node> Nodes;         //node vector for in-memory algorithms
        vector<vector<Edge>> Edges;     //adjacent edges list for in-memory algorithms
        vector<vector<MCEdge>> MCEdges; //adjacent mc-edges list for in-memory algorithms
        string QueryType;       //query type
        uint visited_number = 0;//the number of visited vertices
        unordered_set<NodeId> set_visited;  //the set of visited vertices
        double query_time = 0;  //query time, ms
        vector<vector<Distance>> mc_cost;   //vector of cost
        vector<vector<NodeId>> mc_pre;       //vector of predecessor id
        string query_type = "all";  // S M L all R

        /*** Basic Functions for Data Load ***/
        void ReadGraph_W();
        void MCReadGraph_Edges(const string& filename);  //function for reading edges
        void MCReadGraph_MCEdges(string filename);//read MCgraph edges
        void MCReadGraph_MCEdgesS(const string& filename);

        /*** In-memory SP algorithms ***/
        void Initiation();  //initiation function for graph setup
        Distance Dijkstra(NodeId node_start, NodeId node_end);
        Distance BiDijkstra(NodeId node_start, NodeId node_end);
//        Distance DijkstraCH(NodeId node_start, NodeId node_end);
        double MC_Multipass_one(const string& filename, bool ifIO, int algo_choice);//Function for MCSP by Dijkstra
        double MC_Multipass(const string& qtype, bool ifIO, int algo_choice);//Function for MCSP by Dijkstra, return average query time
        void MC_Evaluate(string qtype);//overall evaluation for all in-memory MCSP algorithms

        /*** I/O efficient SP algorithms ***/
//        void Initiation_IO();
//        Distance Dijkstra_IO(int node_start, int node_end);
//        Distance MC_OneStop_IO(int node_start, int node_end);
//        Distance MC_MultiStop_IO(int node_start, int node_end);
//        void slack_no_read(int cri_i, int item_id);//, vector<set<int>> & visited
//        void slack(int cri_i, int item_id);//function for slacking the node item_id of criteria cri_i
//        void one_stop_strategy(set<int> & set_remain, int & node_end);//unordered_
//        void multiple_stop_strategy(set<int> & set_remain, int & node_end);//, vector<pair<bool,int>> if_stop
//        double MC_OnePass_one(string filename, string qtype, int algo_choice);
//        double MC_OnePass(string qtype,int algo_choice);
//        void MC_Evaluate_IO(string qtype);//overall evaluation for all MCSP algorithms
//        bool JudgeEmpty();

        /*** Supportive functions for SP algorithms ***/
        void Pathshow(list<int>);//show the path sequence in console
        list<int> Dij_getPath(vector<NodeId> & pre, NodeId ID1,NodeId ID2);
        list<int> BiDij_getPath(vector<NodeId> & pre,vector<NodeId> & pre_b,NodeId ID1,NodeId ID2,NodeId terminate_id);

    };
}


#endif //GRAPH_H
