/*
 * Filename:    head.h
 * Description: basic head file for graph data processing(both in-memory and external-memory)
 * Created:     16 Sep 2021
 * Authors:     Xinjie ZHOU
 */

#ifndef HEAD_H_
#define HEAD_H_

#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
#include <vector>
#include <list>
#include <queue>
#include <set>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <tuple>
#include <algorithm>
#include <random>
#include <assert.h>
#include <bitset>
#include <sys/stat.h>
#include <sys/types.h>
#include "Timer.h"

using namespace std;

typedef int NodeId;
typedef int EdgeWeight;
typedef unsigned int uint;
typedef unsigned short int usint;
typedef unsigned short int RangeId;
//typedef unsigned int uint;
typedef unsigned int Distance;//long long
typedef tuple<int,int,int,int,int> MCWeight;

#define OneHop 1
#define MultiHops 2
#define BiMultiHops 3
#define OneHopNoIO 4
#define MultiHopsNoIO 5
#define BiMultiHopsNoIO 6
#define EMDijk 7
#define DijkstraIO 8
#define BiDijkstraIO 9
#define FORWARD 1
#define REVERSE 2
#define LONG_DIS 1
#define MEDIUM_DIS 2
#define SHORT_DIS 3
#define RANDOM_DIS 4
#define EXP1 1
#define EXP2 2
#define EXP3 3
#define EXP4 4
#define EXP5 5
#define EXP6 6
#define EXP7 7
#define INF INT32_MAX
#define Edge_SZ 16 //the memory size of EdgePairW, 3+1 int (16 bytes)
#define MCEdge_SZ 28 //the memory size of MCEdgeT, 7 int (28 bytes)
#define NUM_OF_CRITERIA 5 //the maximal number of criteria in the graph data


/*** Parameters for specifying ***/
const char* disk_conf = "disk=/Users/zhouxj/stxxl.tmp, 20 GiB, syscall unlink"; //disk configuration for stxxl
//For Windows: "disk=e:\stxxl.tmp, 20 GiB, wincall delete";
char* DataPath =  "/Users/zhouxj/Documents/1-Research/Datasets/"; //the data path of datasets
//For Windows: "E:/Datasets/";

///Basic information
int num_of_cri = 3;         //the number of criteria to be processed
int run_times = 10;        //50, 30, 20
string query_type = "R";  // S M L R
bool ifOptimal = false;     // flag of computing the optimal IO
bool ifMu = false;          //if test the parameter of Mu
uint PartitionNumber = 2048;    //default 2048, th number of partitions
uint PartitionSize = 256;//4096   //default 1024 KB, the memory size of one partition, unit(KB) 1024
const uint Block_SZ = 262144;   //default 262144 bytes(256KB), the memory size of io buffer block, unit(byte)
const uint WeightPowMax = 20;   //default 20, the power of maximal edge weight, 2^20=1,048,576(1.05 million)
const uint MinMemThreshold = 128; //default 256 MB, the minimal threshold of memory size for graph data
const uint MuForEM = 3;  //default 0.3, the memory portion for external vector
const uint MuForEM_Bi = 3;  //default 0.3, the memory portion for external vector
double alpha = 0.6;     // threshold of One-Hop algorithm
double alpha_multi = 0.6;  // threshold of Multi-Hop algorithm
double alpha_bi = 1;  // threshold of BiMulti-Hops algorithm
bool ifFile = true;


///Graph data memory
const uint MemGraph_EMDijk = 1326;//15013;//9837;//4893;//   2178;  //the memory size for graph data(i.e. adjacency lists), unit(MB)
const uint MemGraph_DijkIO = 1349;//15095;//9892;// 5036;//  2960;  // the memory size for graph data(i.e. adjacency lists), unit(MB)
const uint MemGraph_BiDijkIO = 1235; //15030;//9804;//5021; //     // the memory size for graph data(i.e. adjacency lists), unit(MB)
const uint MemGraph_NoIO = 1103;//14927;//9687;//  4948;// 2817;  // the memory size for graph data(i.e. adjacency lists), unit(MB)
const uint MemGraph_NoIO_Bi = 664;//14637;//9327;//4808;//  2594;//MinMemThreshold;  //the memory size for graph data(i.e. adjacency lists), unit(MB)
const uint MemGraph_IO = 1080; //14915;//9669;// 4945; //     // the memory size for graph data(i.e. adjacency lists), unit(MB)
const uint MemGraph_IO_Bi = 641; //14624;//9309;// 4805; //    //the memory size for graph data(i.e. adjacency lists), unit(MB)

///for EM_Dijk
unsigned long long HotPool_EdgeSZ = 0.9*MemGraph_EMDijk*1024*1024/(Edge_SZ*WeightPowMax);  //default 139810, the number of Edges in each in-memory hot pool, overall memory cost: BASIC_SZ
///for algorithms with io optimization
const uint BlockPerPage = 1; //default 4, the memory size of one page in external vector, 4*256KB=1MB
const uint PageNumberEMDijk = (MemGraph_EMDijk / (10 * WeightPowMax)) * 4;//the number of pages in memory
const uint PageNumberIO = (MemGraph_IO * MuForEM /10)*4;//32*4; //default 32, the number of pages in memory, overall memory cost for one vector is 32*1=32MB
const uint PageNumberIO_Bi = (MemGraph_IO * MuForEM_Bi / 10)*4;//32*4; //default 32, the number of pages in memory, overall memory cost for one vector is 32*1=32MB
uint Partition_N = 0.85*(MemGraph_IO-PageNumberIO/4)*1024/PartitionSize;// the number of partitions in memory/
uint Partition_N_Bi = 0.85*(MemGraph_IO_Bi-PageNumberIO_Bi/4)*1024/PartitionSize;// the number of partitions in memory, memory cost: BASIC_SZ

///for algorithms without io optimization
//const uint PageNumber_DijkIO = MemGraph_DijkIO*1024/PartitionSize;//the number of pages in stxxl vector for algorithms without io optimization, memory cost: BASIC_SZ
const uint PageNumber_DijkIO = (MemGraph_DijkIO*1024)/(BlockPerPage*Block_SZ/1024);//the number of pages in stxxl vector for algorithms without io optimization, memory cost: BASIC_SZ MemGraph_DijkIO
const uint PageNumber_BiDijkIO = (MemGraph_BiDijkIO*1024)/(BlockPerPage*Block_SZ/1024);//the number of pages in stxxl vector for algorithms without io optimization, memory cost: BASIC_SZ MemGraph_DijkIO
const uint PageNumber_NoIO = MemGraph_NoIO*1024/(BlockPerPage*Block_SZ/1024);//the number of pages in stxxl vector for algorithms without io optimization, memory cost: BASIC_SZ
const uint PageNumber_NoIO_Bi = MemGraph_NoIO_Bi*1024/(BlockPerPage*Block_SZ/1024);//the number of pages in stxxl vector for algorithms without io optimization, memory cost: BASIC_SZ

///External priority queue (overall memory cost for one priority queue is about 50MB in practice)
const uint MEMORY_FOR_PRIORITY_QUEUE = 32*1024*1024;//default 32 MB, internal memory for priority queue
const uint PRIORITY_QUEUE_MAX_SIZE = 32*1024*1024; //default 32*1024*1024(33,334,432), maximal number of items in memory, the actual number can be somewhat larger
const uint mem_for_pools = 16*1024*1024;    //default 16 MB, internal memory for prefetching
///External sorter
const uint MEMORY_FOR_SORTER = 256*1024*1024; //default 256 MB, internal memory for sorter

///Intermediate parameters
const uint Partition_MCEdgeSZ = PartitionSize*1024/MCEdge_SZ;//default 9362, the number of MCEdges that one partition can hold at most, 9362*28bytes=256KB
const uint Vector_EdgeSZ = Block_SZ/Edge_SZ;      //default 21845, the number of Edges of one page in external vector, one Edge is 12 bytes, 21845*12bytes=256KB
const uint Vector_MCEdgeSZ = Block_SZ/MCEdge_SZ;  //default 9362, the number of MCEdges of one page in external vector, one MCEdge is 28 bytes, 9362*28bytes=256KB


//Data structure for nodes in graph
struct Node {
    /*---Basic information for graph---*/
//    NodeId id=-1;                 //node ID
    uint degree;
    /*---Functions---*/
    Node():degree(0){}
//    Node(NodeId tid):id(tid){}
    ~Node(){}
};
//Data structure for edges in single criteria graph
struct Edge {
    NodeId u=0; //the node ID of end node
    EdgeWeight w=0; //the weight of edge
    Edge(){}
    Edge(NodeId u_, EdgeWeight w_) :u(u_), w(w_) {}//initiation function
    ~Edge() {}
};
//Data structure for edges in MC-graph
struct MCEdge {
    NodeId u=0;          //the node ID of end node
    std::vector<EdgeWeight> w=std::vector<EdgeWeight>();  //the weight of edge
    MCEdge(NodeId u_, std::vector<EdgeWeight> w_=std::vector<EdgeWeight>()) :u(u_), w(w_) {}//initiation function
    MCEdge() {}//u = -1; w.resize(num_of_cri);
    ~MCEdge() {w.clear();}
};
//Data structure for edges in MC-graph
struct MCEdgeT {
    NodeId ID1 = 0;     //the ID of source node
    NodeId ID2 = 0;     //the ID of target node
    MCWeight w;  //the weight of edge

    MCEdgeT() {}
    MCEdgeT(NodeId ID1_, NodeId ID2_, int w0,int w1,int w2,int w3, int w4):ID1(ID1_),ID2(ID2_) {
        get<0>(w) = w0; get<1>(w) = w1; get<2>(w) = w2;
        get<3>(w) = w3; get<4>(w) = w4;
    }
    MCEdgeT(NodeId ID1_, NodeId ID2_, MCWeight w_) :ID1(ID1_),ID2(ID2_), w(w_) {}//initiation function
    ~MCEdgeT() {}
    int getW(int index){
        assert(index<=4 && index >=0);
        switch (index) {
            case 0: return get<0>(w);
            case 1: return get<1>(w);
            case 2: return get<2>(w);
            case 3: return get<3>(w);
            case 4: return get<4>(w);
            default: return -1;
        }
    }
    void putW(int index, int newW){
        assert(index<=4 && index >=0);
        switch (index) {
            case 0: get<0>(w) = newW; break;
            case 1: get<1>(w) = newW; break;
            case 2: get<2>(w) = newW; break;
            case 3: get<3>(w) = newW; break;
            case 4: get<4>(w) = newW; break;
            default: break;
        }
    }
};
//Data structure for priority queue
struct VertexCost{
    NodeId id = 0;
    Distance cost = 0;
    VertexCost(){}
    VertexCost(NodeId id_,Distance cost_):id(id_),cost(cost_){}
    ~VertexCost(){}
};
//Data Structure for edge pairs
struct EdgePair{
    NodeId ID1=0;
    NodeId ID2=0;
    EdgePair(){}
    EdgePair(NodeId ID1_,NodeId ID2_):ID1(ID1_),ID2(ID2_){}
    ~EdgePair(){}
};
//Data structure for edge pairs with weight
struct EdgePairW{
    NodeId ID1=0;
    NodeId ID2=0;
    EdgeWeight w=0;
    EdgePairW(){}
    EdgePairW(NodeId ID1_,NodeId ID2_,EdgeWeight w_):ID1(ID1_),ID2(ID2_),w(w_){}
    ~EdgePairW(){}
    void update(NodeId ID1_,NodeId ID2_,EdgeWeight w_){
        ID1=ID1_;
        ID2=ID2_;
        w=w_;
    }
};
//Comparator for EdgePair
struct EdgePairComparator {
    bool operator () (const EdgePair& a, const EdgePair& b) const {
        if(a.ID1 < b.ID1){
            return true;
        }else{
            if(a.ID1==b.ID1){
                if(a.ID2 < b.ID2){
                    return true;
                }else{
                    return false;
                }
            }else{
                return false;
            }
        }
    }
    EdgePair min_value() const {
        EdgePair a;
        a.ID1=std::numeric_limits<NodeId>::min();
        a.ID2=std::numeric_limits<NodeId>::min();
        return a;
    }
    EdgePair max_value() const {
        EdgePair a;
        a.ID1=std::numeric_limits<NodeId>::max();
        a.ID2=std::numeric_limits<NodeId>::max();
        return a;
    }
};
//Comparator of EdgePairW, compare OD pair
struct EdgePairWComparator {//the less, the better
    bool operator () (const EdgePairW& a, const EdgePairW& b) const {
        if(a.ID1 < b.ID1){
            return true;
        }else{
            if(a.ID1==b.ID1){//if ID1 == ID1
                if(a.ID2 < b.ID2){
                    return true;
                }else {
                    if(a.ID2==b.ID2){//if ID2 == ID2
                        if(a.w < b.w){
                            return true;
                        }else{
                            return false;
                        }
                    }else {
                        return false;
                    }
                }
            }else{
                return false;
            }
        }
    }
    static EdgePairW min_value() {
        EdgePairW a;
        a.ID1=std::numeric_limits<NodeId>::min();
        a.ID2=std::numeric_limits<NodeId>::min();
        a.w=std::numeric_limits<EdgeWeight>::min();
        return a;
    }
    static EdgePairW max_value() {
        EdgePairW a;
        a.ID1=std::numeric_limits<NodeId>::max();
        a.ID2=std::numeric_limits<NodeId>::max();
        a.w=std::numeric_limits<EdgeWeight>::max();
        return a;
    }
};
struct PairComparator {
    bool operator () (const std::pair<NodeId,uint>& a, const std::pair<NodeId,uint>& b) const {//
        return a.second<b.second;
    }
    std::pair<NodeId,uint> min_value() const {
        std::pair<NodeId,uint> a;
        a.first=0;
        a.second=std::numeric_limits<uint>::min();
        return a;
    }
    std::pair<NodeId,uint> max_value() const {
        std::pair<NodeId,uint> a;
        a.first=0;
        a.second=std::numeric_limits<uint>::max();
        return a;
    }
};
struct PQCompareLess
{
    bool operator () (const VertexCost& a, const VertexCost& b) const
    { return (a.cost > b.cost); }

    VertexCost min_value() const
    { return VertexCost(0,std::numeric_limits<Distance>::max()); }
};
struct PQEdgePairCompareLess
{
    bool operator () (const pair<uint,EdgePair>& a, const pair<uint,EdgePair>& b) const
    { return (a.first > b.first); }

    pair<uint,EdgePair> min_value() const
    {   pair<uint,EdgePair> a;
        a.first=std::numeric_limits<uint>::max();
        a.second=EdgePair();
        return a; }
};

//created by Mengxuan, modified by Xinjie
namespace benchmark {

#define NULLINDEX 0xFFFFFFFF

    template<int log_k, typename id_t, typename k_t >//id,value
    class heap {
    public:

        // Expose types.
        typedef k_t key_t;
        typedef id_t node_t;

        // Some constants regarding the elements.
        //static const node_t NULLINDEX = 0xFFFFFFFF;
        static const node_t k = 1 << log_k;//equals k = 1*2^log_k, usually log_k = 2, so k = 4

        // A struct defining a heap element.
        struct element_t {
            key_t key;
            node_t element;

            element_t() : key(0), element(0) {}

            element_t(const key_t k, const node_t e) : key(k), element(e) {}
        };


        //public:

        // Constructor of the heap.
        heap(node_t n) : n(0), max_n(n), elements(n), position(n, NULLINDEX) {//n is the number of elements in current
            // state, max_n is the size of heap
        }

        heap(): n(0), max_n(0), elements(0), position(0, NULLINDEX) {}

        ~heap(){}

        // Risize the heap
        inline void resize(node_t a){
            n = 0; max_n = a;
            elements.resize(a);
            position.resize(a, NULLINDEX);
        }

        // Size of the heap.
        inline node_t size() const {
            return n;
        }

        // Heap empty?
        inline bool empty() const {
            return size() == 0;
        }

        // Extract min element.
        inline void extract_min(node_t &element, key_t &key) {
            assert(!empty());

            element_t &front = elements[0];

            // Assign element and key.
            element = front.element;
            key = front.key;

            // Replace elements[0] by last element.
            position[element] = NULLINDEX;
            --n;//n=n-1
            if (!empty()) {
                front = elements[n];//elements[n] is the top element
                position[front.element] = 0;//make its position valid, it is also the smallest one
                sift_down(0);
            }
        }

        inline key_t top_key() {//get the key, i.e. minimal cost
            assert(!empty());

            element_t &front = elements[0];

            return front.key;

        }

        inline node_t top_value() {//get the value, i.e. id number of minimal cost

            assert(!empty());

            element_t &front = elements[0];

            return front.element;
        }

        // Update an element of the heap.
        inline void update(const node_t element, const key_t key) {

            if (position[element] == NULLINDEX) {//if originally NULL
                element_t &back = elements[n];//add new element to position n
                back.key = key;
                back.element = element;
                position[element] = n;//set position id to n
                sift_up(n++);
            } else {//if already valid, update the value
                node_t el_pos = position[element];//position information
                element_t &el = elements[el_pos];//get the element
                if (key > el.key) {//update the element
                    el.key = key;
                    sift_down(el_pos);
                } else {
                    el.key = key;
                    sift_up(el_pos);
                }
            }
        }

        // Clear the heap.
        inline void clear() {
            for (node_t i = 0; i < n; ++i) {
                position[elements[i].element] = NULLINDEX;
            }
            n = 0;
        }

        // Cheaper erase.
        inline void erase(node_t v) {
            position[v] = NULLINDEX;
        }

        inline void clear_n() {
            n = 0;
        }

        // Test whether an element is contained in the heap.
        inline bool contains(const node_t element) const {
            return position[element] != NULLINDEX;
        }

        //return current elements information
        void get_elements(std::vector<std::pair<int,int>> &e_vector){
            std::pair<int,int> temp_pair;

            for(int i=0;i<n;i++){
                temp_pair.first = elements[i].key;
                temp_pair.second = elements[i].element;
                e_vector.push_back(temp_pair);
            }
        }

    protected:

        // Sift up an element.
        inline void sift_up(node_t i) {
            assert(i < n);
            node_t cur_i = i;
            while (cur_i > 0) {
                node_t parent_i = (cur_i - 1) >> log_k;//equals (cur_i - 1)/(2^log_k)
                if (elements[parent_i].key > elements[cur_i].key)//compare with parent node, if smaller, then swap
                    swap(cur_i, parent_i);
                else
                    break;
                cur_i = parent_i;
            }
        }

        // Sift down an element.
        inline void sift_down(node_t i) {
            assert(i < n);

            while (true) {
                node_t min_ind = i;
                key_t min_key = elements[i].key;

                node_t child_ind_l = (i << log_k) + 1;//equals i*2^log_k + 1
                node_t child_ind_u = std::min(child_ind_l + k, n);//equals min(child_ind_l+4,n)

                for (node_t j = child_ind_l; j < child_ind_u; ++j) {
                    if (elements[j].key < min_key) {
                        min_ind = j;
                        min_key = elements[j].key;
                    }
                }

                // Exchange?
                if (min_ind != i) {
                    swap(i, min_ind);
                    i = min_ind;
                } else {
                    break;
                }
            }
        }

        // Swap two elements in the heap.
        inline void swap(const node_t i, const node_t j) {
            element_t &el_i = elements[i];
            element_t &el_j = elements[j];

            // Exchange positions
            position[el_i.element] = j;
            position[el_j.element] = i;

            // Exchange elements
            element_t temp = el_i;
            el_i = el_j;
            el_j = temp;
        }

    private:

        // Number of elements in the heap.
        node_t n;

        // Number of maximal elements.
        node_t max_n;

        // Array of length heap_elements.
        std::vector<element_t> elements;

        // An array of positions for all elements.
        std::vector<node_t> position;
    };
}

#endif
