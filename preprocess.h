/*
 * Filename:    preprocess.h
 * Description: head file for graph data preprocessing
 * Created:     12 Jan 2022
 * Author:      Xinjie ZHOU
 */
#ifndef PREPROCESS_H
#define PREPROCESS_H

#include "head.h"
#include "graph.h"
#include "ioxxl.h"

namespace gbpre{
    using namespace std;
    using namespace ioxxl;

    RankSorter rankSorter(PairComparator(),MEMORY_FOR_SORTER);    //edge rank sorter

    template <size_t N, class T = unsigned long int>
    class bitset1;
    template <size_t N, class T = unsigned long int>
    class bitreference {
        typedef bitreference<N,T> Self;
        friend class bitset1<N,T>;
        bitreference();
        bitreference(bitset1<N,T>& ref, size_t pos):m_ref(ref),m_pos(pos){};
    public:
        Self& operator=(bool x) {
            m_ref.set(m_pos,x);
            return *this;
        };
        operator bool() const {
            return m_ref.test(m_pos);
        }
        Self& operator = (const Self& other) {
            this->operator =( bool(other) );
            return *this;
        }
    private:
        bitset1<N,T>& m_ref;
        const size_t m_pos;
    };
    template <size_t N, class T>
    class bitset1 {
        typedef bitset1<N,T> Self;
    public:
        typedef bitreference<N,T> reference;
        explicit bitset1(T val = 0):m_val(val){assert( N/8<=sizeof(T) );};
        Self& set() {//set the done bit to true, all the other to false
            m_val = 1;
            return *this;
        };
        Self& set(size_t pos, bool val) {//set certain pos to true or false
            if (val) m_val|=(1<<pos);
            else m_val &= ~(1<<pos);
            return *this;
        };
        bool test(size_t pos) const {//test certain pos
            return m_val&(1<<pos);
        }
        reference operator[](size_t pos) {
            assert(pos<N);
            return reference(*this,pos);
        };
        bool operator[](size_t pos) const {
            return test(pos);
        };
        Self& reset() {
            m_val = 0;
            return *this;
        };
        static size_t size() {
            return N;
        };
        bool operator==(const bitset1<N,T>& other) const {
            return m_val == other.m_val;
        };
        bool operator!=(const bitset1<N,T>& other) const {
            return !this->operator==(other);
        }
//    private:
        T m_val = 0;
    };

    class Graph_pre:public gbstd::MCGraph
    {
    public:
        /*---Graph information*/
        int num_criteria = 0;
//    vector<vector<MCEdge>> MCEdges2;//adjacent mc-edges list for in-memory algorithms
        vector<map<NodeId,Edge>> EdgesMap;    //adjacent edges list for in-memory algorithms, memory usage: 4*|E|
        vector<map<NodeId,MCEdge>> MCEdgesMap;//adjacent mc-edges list for in-memory algorithms, use map to make sure the consistency of undirected edges
//        vector<vector<MCEdge>> MCEdgesVec;//adjacent lists for partitioning
//        unordered_map< pair<int, int>, int > edge_map;//for removing the repeat edges
        pair<set<NodeId>,uint> max_component;//recording the vertices ids of the maximal connected component, <Nodes_of_component, Number_of_edges>, memory usage: |V|
//        unordered_map<NodeId,NodeId> nodes_map;  //the map from old vertex id to the new one, for component computation
//        map<NodeId,NodeId> nodes_map;  //the map from old vertex id to the new one, for component computation
        uint max_degree = 0;
        uint min_degree = INF;
        double ave_degree = 0;
//        map<int,int> old_to_new;//map the old id of vertices to the new one
//        map<pair<int,int>,vector<int>> edge_list;//new edge list: after partitioning
//        vector<int> partitions_info;//vector for storing the info of partitions, node number, start vertex id, end vertex id
        vector<vector<int>> cluster_to_node;//from cluster to vertices
//        vector<int> new_to_old;
        int base_max = 0;
        int base_min = INT32_MAX;
        vector<int> max_weights;    //the maximal edge weights of all criteria
        vector<int> min_weights;    //the minimal edge weights of all criteria
        uint PartitionMCEdge_NUM = PartitionSize*1024/MCEdge_SZ; //default 37449, the number of MCEdges in one partition
        string partMethod="DC";//DC: determine clustering; Metis

        string strategy = "wave";
//        double zeta[5] = {0.27273,0.27273,0.27273,0.09091,0.09091};
//        double zeta[5] = {0.5,0.3,0.2,0,0};
        double zeta[5] = {0.33,0.33,0.33,0,0};

        int threadnum=15;
        vector<Semaphore*> vSm;
        Semaphore* sm;

        /*** Basic Functions & Data Preprocessing ***/
        void Read(const string& file,int num_line);
        void ReadGraph_edges(const string& filename,bool ifBase);  //function for reading edges
        /*** Data Cleaning & Component Progressing ***/
//        void ReadCleaningW(const string& r_file);//deal with the invalid, repeat or missing edges
//        void Cleaning_Process();//deal with the repeat and missing edges
//        void Component_Process(Graph_pre & gr_t);//count the number of connected component and retain the largest one
//        template <class T>
//        void Component_DFS(T & EdgesMap2);   //function of finding the maximal connected component of the graph
////    void Component_ReadGraph(bool ifCoordinate,Graph_pre & gr_t);
//        void Component_ToMC(set<NodeId>& w_set,uint new_edge_num, Graph_pre & gr_t);
//        void Component_WriteEdge(const string& filename,set<NodeId>& w_set,uint new_edge_num,Graph_pre & gr_t);
        template <class T>
        pair<int, unsigned long long> DFS_CC(T & Edges, unordered_set<int> & set_A, set<int> & set_LCC, int nodenum);
        pair<int, unsigned long long> DFS_CC_MC(vector<vector<MCEdge>> & Edges, unordered_set<int> & set_A, set<int> & set_LCC, int nodenum);

        /// new data generation
        void BaseCriGenerate(uint n_num); //function for generating edge weights of basic criterion for graph
        void WriteGraph(string filename, set<int>& set_LCC, vector<unordered_map<int,vector<int>>>& NeighborMap, int cri_i);
        void MC_Generation(int num, set<int>& set_LCC, vector<unordered_map<int,vector<int>>>& NeighborMap);

        void GraphPartitioning();
        void ReadMCEdges(const string& filename);
        void IDMinusOne(string filename);

        /// For Other Partition Methods
        void PartitionPreprocess(string parMethod);
        void KaHyParConvertInput();
        void METISConvertInput();/// For Metis

        void PartitionPostprocess(string parMethod, int partitionNum);

        /*** Data Generation ***/
        void MC_DataGeneration(int num);//
        void MC_Generate(int num);
        void MC_WriteEdges(bool ifCluster, bool ifOriginal);
        void MC_WriteEdgesS(bool ifCluster, bool ifOriginal);
        void MC_WriteEdgesIO(bool ifCluster, bool ifOriginal);
        void MC_ToEdgesS();
//        void text_to_bin();
//        void text_to_bin_clusters();
        void ODpairGenerate(int times);
        void ODpairGenerate_Dis(int times, bool ifRead);
        void CoverageDijkstra(int node_start, vector<vector<pair<int,int>>>& queries, unsigned long long lower_bound, unsigned long long upper_bound);
        void CoverageDijkstraV(vector<int>& p, vector<vector<pair<int,int>>>& queries, unsigned long long lower_bound, unsigned long long upper_bound);
        void ODpairGenerate_Dis2(int times, bool ifRead);
        void ODpairGenerate_Remap();


        /*** Data Rearrangement ***/
        void Read_ToMCEdgesMap(const string& filename,bool ifBase);
        void Read_MCEdges(const string& filename);

        void Read_MCEdgesMap(const string& filename);
        void Read_MCEdgesMapS(const string& filename);
        bool DeterClustering(bool ifReuse);
        bool MinSpanTree(int node_start,EdgeSorter & mstEdges,vector<uint> & outDegree);
        void DFS_EulerTour(NodeId item_id,vector<EdgePair> & mstEdges,vector<uint>& indexID,vector<uint> & outDegree,VectorUint& EulerPath);//list<NodeId>
        void WriteClusters(RankSorter& rankSorter);//,map<NodeId,NodeId>& remap
        void WritePartitionInfo(vector<vector<NodeId>>& partitionNodes);
        void WritePartition(int cluster_i, vector<NodeId>& partitionNodes, bool ifbin);//map<NodeId,NodeId>& remap,
        void EMDIJk_preprocess(bool ifRead);
        int Weight_to_category(int w);


        /*** Data preprocessing for Multiplex networks ***/
        void Read_Multiplex(const string& filename);
        void Write_Multiplex(const string& w_file);
        void Preprocess_Multiplex();
        void Component_WriteMultiplex(set<NodeId>& w_set,uint new_edge_num);
    };
}


#endif //PREPROCESS_H
