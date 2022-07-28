/*
 * Filename:    emgraph.h
 * Description: head file for external-memory graph data processing
 * Created:     04 March 2021
 * Authors:     Xinjie ZHOU
 */

#ifndef MCSPS_EMGRAPH_H
#define MCSPS_EMGRAPH_H

#include <stxxl/vector>
#include <stxxl/priority_queue>
#include <stxxl/sorter>
#include "Timer.h"
#include "head.h"
#include "ioxxl.h"

namespace gbxxl{
    using namespace std;
    using namespace ioxxl;

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
        typedef gbxxl::bitreference<N,T> reference;
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
        bool operator==(const gbxxl::bitset1<N,T>& other) const {
            return m_val == other.m_val;
        };
        bool operator!=(const gbxxl::bitset1<N,T>& other) const {
            return !this->operator==(other);
        }
//    private:
        T m_val = 0;
    };

    /*** Collective Variables ***/
    uint bufferSize = Block_SZ;
    const uint sizeInt = sizeof(int);
    const uint intsPerBlock = bufferSize / sizeInt;//max number of int in buffer
    const uint MCEDGE_PER_BLK = 4096 / MCEdge_SZ;//maximal number of MCEdge per block

    ///struct
    //data structure for node index in disk, remain in memory
    struct nodeindex{
        int nodeid;//node id
        int offset;//used to locate the position of this node info in graph file
//        int degree;//node degree
    };
    //Class of LRU strategy
    class LRUCache {
    public:
        unsigned int capacity_ = 0;
        list<pair<int, int>> ls_;
        unordered_map<int, list<pair<int, int>>::iterator> hash_;

        LRUCache(){};
        LRUCache(int capacity) : capacity_(capacity) {}
        void clear(){
            ls_.clear();
            hash_.clear();
        }
        unsigned int size(){
            return ls_.size();
        }
        //function of erasing certain key
        bool erase(int key){
            if (hash_.find(key) == hash_.end())//if not found
                return false;
            else {//if found
                ls_.erase(hash_[key]);//erase the old element
                hash_.erase(key);
                return true;
            }
        }
        //function of popping the least important element
        int pop(){
            int pop_id = -1;
            if (!ls_.empty()) {//if not empty
                pop_id = ls_.back().first;
                hash_.erase(ls_.back().first);
                ls_.pop_back();
            }
            return pop_id;
        }
        void set_capacity(unsigned int capacity) {
            capacity_ = capacity;
        }
        //function of putting key to the back, i.e., the first to be evicted
        int push_back(int key){
            if(capacity_ == 0)
                return -1;
            if (hash_.find(key) == hash_.end())//if not found
                return -1;
            else {//if found
                int value = hash_[key]->second;
                ls_.erase(hash_[key]);//erase the old element
                ls_.emplace_back(key, value);//re-add the same element to the back
                hash_[key] = std::prev(ls_.end());  //the iterator of the last element of ls_
                return value;
            }
        }
        //function for finding the element
        int get(int key) {
            if (hash_.find(key) == hash_.end())//if not found
                return -1;
            else {//if found
                int value = hash_[key]->second;
                ls_.erase(hash_[key]);//erase the old element
                ls_.push_front(make_pair(key, value));//re-add the same element to the front
                hash_[key] = ls_.begin();
                return value;
            }
        }
        //function of updating info of elements
        int put(int key, int value) {
            if(capacity_ == 0)
                return -1;
            int temp_value=0,pop_id=-1;
            if (hash_.find(key) != hash_.end()){//if found, erase the element in list
                temp_value = hash_[key]->second;
                ls_.erase(hash_[key]);
            } else if (ls_.size() >= capacity_) {//if not found and exceed the capacity, delete the back element of list and map
                pop_id = ls_.back().first;
                hash_.erase(ls_.back().first);
                ls_.pop_back();
            }
            ls_.push_front(make_pair(key, value+temp_value));//add the element to the front of list
            hash_[key] = ls_.begin();//store info to map
            return pop_id;
        }
    };
    //Class of read buffer
    class ReadBuffer {
        int* buf;
        int ptr;
        FILE* infile;

    public:
        bool isend;
        int  nread;  // number of records read into memory
        int io_number = 0;
        double io_time = 0;//s

        ReadBuffer() {
            buf = (int*)malloc(bufferSize);
        };
        ReadBuffer(FILE* infile_t) {
            buf = (int*)malloc(bufferSize);
            open(infile_t);
        };
        ReadBuffer(FILE* infile_t, int offset) {
            buf = (int*)malloc(bufferSize);
            open(infile_t, offset);
        }
        ~ReadBuffer() {
            free(buf);
        };

        int read(int* p) {
            if (isend) return 0;
            *p = buf[ptr];

            if ((++ptr) == nread) {
                if (nread < intsPerBlock)//if cannot fulfill the buffer
                    isend = true;
                else {
                    nread = fread(buf, sizeInt, intsPerBlock, infile);
                    ++io_number;
                    if (nread == 0)//if none be read
                        isend = true;
                    ptr = 0;
                }
            }
            return ptr;
        };

        int read(int* p, int cnt) {
            int nread = 0;
            for (int i = 0; i < cnt; ++i) {
                if (this->read(p + i) == 0)
                    break;
                ++nread;
            }
            return nread;
        }

        void open(FILE* infile_t) {
            infile = infile_t;
            rewind(infile);
            nread = fread(buf, sizeInt, intsPerBlock, infile);
            ++io_number;
            ptr = 0;
            isend = false;
        };

        void open(FILE* infile_t, int offset) {//offset is the number of int
            infile = infile_t;
            rewind(infile);
            fseek(infile, offset * sizeInt, SEEK_SET);
            nread = fread(buf, sizeInt, intsPerBlock, infile);
            ++io_number;
            ptr = 0;
            isend = false;
        };

        void relocate(int offset){//relocate the index of infile
            rewind(infile);
            fseek(infile, offset * sizeInt, SEEK_SET);
            nread = fread(buf, sizeInt, intsPerBlock, infile);
            ++io_number;
            ptr = 0;
            isend = false;
        }
    };
    //Class of write buffer
    class WriteBuffer {
        int* buf;//the memory address of write buffer
    public:
        FILE* outfile;//the write file stream
        int ptr;//the pointer
        int io_number = 0;

    public:
        WriteBuffer() {
            buf = (int*)malloc(bufferSize); // change WriteBuffer
        };

        WriteBuffer(FILE* outfile_t) {
            buf = (int*)malloc(bufferSize);
            open(outfile_t);
        };

        WriteBuffer(FILE* outfile_t, int offset) {
            buf = (int*)malloc(bufferSize);
            open(outfile_t, offset);
        };

        ~WriteBuffer() {
            flush();
            free(buf);
        };

        void write(int* src, int cnt) {//src is source data address
            unsigned total = cnt;
            while (cnt > 0) {
                if (cnt < intsPerBlock - ptr) {//it indicates current block can load the data
                    memcpy(buf + ptr, src + (total - cnt), cnt * sizeInt);//copy cnt * sizeInt size memory data from src + (total - cnt) to buf + ptr
                    ptr += cnt;
                    cnt = 0;
                }
                else {//if current block can not hold the data
                    memcpy(buf + ptr, src + (total - cnt), (intsPerBlock - ptr) * sizeInt);
                    fwrite(buf, sizeInt, intsPerBlock, outfile);//write the content of buffer to outfile
                    ++io_number;
                    cnt -= intsPerBlock - ptr;
                    ptr = 0;
                }
            }
        };

        void open(FILE* outfile_t) {
            outfile = outfile_t;
            ptr = 0;
            rewind(outfile);
        }

        void open(FILE* outfile_t, int offset) {
            outfile = outfile_t;
            ptr = 0;
            rewind(outfile);
            fseek(outfile, offset * sizeInt, SEEK_SET);//get the location in io stream by offset
        }

        void flush() {
            if (ptr > 0) {//it indicates there is data(write buffer is not full) to be writen
                fwrite(buf, sizeInt, ptr, outfile);
                ++io_number;
                ptr = 0;
            }
        };
    };
    /*** HotPool ***/
    //hot pool struct for EM_Dijk
    struct HotPool {
        unsigned int capacityEdges = HotPool_EdgeSZ;//Capacity of how many edges can be stored in Edges_IM
        unordered_multimap<NodeId,EdgePairW> Edges_IM;  //in-memory multimap for EM_Dijk, memory cost: memSize/(2*WeightPowMax)
//        unordered_map<NodeId,vector<EdgePairW>> Edges_IM;  //in-memory multimap for EM_Dijk, memory cost: memSize/(2*WeightPowMax)
        vector<pair<NodeId,NodeId>> EM_Map;//the map from cluster id to index id of the first edge and the last edge of this cluster in external vector
        VectorEdges_EMDijk Edges_EM;//memory cost: 32 MB

        HotPool(){};
        ~HotPool(){clear();};
        void set_capacity(unsigned int a){
            capacityEdges = a;
        }
        void clear(){
            Edges_EM.clear();
            EM_Map.clear();
            Edges_IM.clear();
        }
    };
    /*struct HotPool {
        unsigned int capacityEdges = HotPool_EdgeSZ;//Capacity of how many edges can be stored in Edges_IM
        unordered_multimap<NodeId,EdgePairW> Edges_IM;  //in-memory multimap for EM_Dijk, memory cost: memSize/(2*WeightPowMax)
        vector<pair<NodeId,NodeId>> EM_Map;//the map from cluster id to index id of the first edge and the last edge of this cluster in external vector
        VectorEdges_EMDijk Edges_EM;//memory cost: 32 MB

        HotPool(){};
        ~HotPool(){clear();};
        void set_capacity(unsigned int a){
            capacityEdges = a;
        }
        void clear(){
            Edges_EM.clear();
            EM_Map.clear();
            Edges_IM.clear();
        }
    };*/
    //hot pool struct for MUL-Dijk: vector
    struct HotPool2 {
        uint capacityClusters = Partition_N;//Capacity of how many partitions can be store in memory, if one partition size is 1MBB, then 128 partitions equal 128MB
        LRUCache LRU_IMCluster;//LRU manager for clusters
        vector<vector<EdgePairW>> HotPools;  //storing the info of partitions in memory or temp file. <partition_id,<adjacency lists>>
        vector<bool> FlagIM;    //the flag of whether the partition is in memory

        HotPool2(){};
        HotPool2(uint pn,uint capacity){
            HotPools.assign(pn,vector<EdgePairW>());
            FlagIM.assign(pn, false);
            set_capacity(capacity);
        };
        ~HotPool2(){clear();};
        void set_capacity(uint cluster_n){//,uint edgesPerCluster_n
            capacityClusters = cluster_n;
            LRU_IMCluster.set_capacity(capacityClusters);
        }
        void clear(){
            LRU_IMCluster.clear();
            HotPools.clear();
            FlagIM.clear();
        }
    };
    //hot pool struct for MUL-Dijk: unordered_multimap
    struct HotPool3 {
        uint capacityClusters = Partition_N;//Capacity of how many partitions can be store in memory, if one partition size is 1MBB, then 128 partitions equal 128MB
        LRUCache LRU_IMCluster;//LRU manager for clusters
        vector<unordered_multimap<NodeId,Edge>> HotPools;  //storing the info of partitions in memory or temp file. <partition_id,<adjacency lists>>
        vector<bool> FlagIM;    //the flag of whether the partition is in memory

        HotPool3(){};
        HotPool3(uint pn,uint capacity){
            HotPools.assign(pn,unordered_multimap<NodeId,Edge>());
            FlagIM.assign(pn, false);
            set_capacity(capacity);
        };
        ~HotPool3(){clear();};
        void set_capacity(uint cluster_n){//,uint edgesPerCluster_n
            capacityClusters = cluster_n;
            LRU_IMCluster.set_capacity(capacityClusters);
        }
        void clear(){
            LRU_IMCluster.clear();
            HotPools.clear();
            FlagIM.clear();
        }
    };
    //hot pool struct for MUL-Dijk: unordered_multimap + stxxl vector
    template <class T>
    struct HotPool4{
        uint capacityClusters = 0;//Capacity of how many partitions can be store in memory, if one partition size is 1MBB, then 128 partitions equal 128MB
        LRUCache LRU_IMCluster;//LRU manager for clusters
//        vector<unordered_multimap<NodeId,MCEdgeT>> HotPools;  //storing the info of partitions in memory or temp file. <partition_id,<adjacency lists>> Edge
        vector<unordered_map<NodeId,vector<MCEdgeT>>> HotPools;  //storing the info of partitions in memory or temp file. <partition_id,<adjacency lists>> Edge
        vector<char> clusterStatus;    //the flag of where the partition is, D: in disk; I: in memory; S: in stxxl vector;
        T MCEdges_EM;//memory cost: 32 MB VectorMCEdgesEMTuple4

        HotPool4(){};
        HotPool4(uint pn,uint capacity){
//            HotPools.assign(pn,unordered_multimap<NodeId,MCEdgeT>());//Edge
            HotPools.assign(pn,unordered_map<NodeId,vector<MCEdgeT>>());//Edge
            clusterStatus.assign(pn, 'D');
            set_capacity(capacity);
            capacityClusters=capacity;
        };
        ~HotPool4(){clear();};
        void set_capacity(uint cluster_n){//,uint edgesPerCluster_n
            capacityClusters = cluster_n;
            LRU_IMCluster.set_capacity(capacityClusters);
        }
        void clear(){
            LRU_IMCluster.clear();
            HotPools.clear();
            clusterStatus.clear();
        }
    };
    //hot pool struct for proposed methods
    template <class T>
    struct MCHotPool{
        uint capacityClusters = 0;//Capacity of how many partitions can be store in memory, if one partition size is 1MBB, then 128 partitions equal 128MB Partition_N
        LRUCache LRU_IMCluster;//LRU manager for clusters
        vector<tuple<char,int,int>> clusterStatus;    //recording the status of clusters, <storage type, original size, current size>, for the first element: D: in disk; I: in memory; S: in stxxl vector;
        vector<vector<MCEdgeT>> HotPools;  //storing the info of partitions in memory or temp file. <partition_id,<adjacency lists>>
        T MCEdges_EM;    //the external vector, 32MB
//        VectorMCEdgesEMTuple_IO MCEdges_EM;    //the external vector, 32MB
//        set<int> im_partitions; //set for recording the partitions in memory, just for observation

        MCHotPool(uint capacity){
            capacityClusters = capacity;
            LRU_IMCluster.set_capacity(capacityClusters);
        }
        void init(uint p_number){
            HotPools.assign(p_number,vector<MCEdgeT>());
            clusterStatus.assign(p_number,make_tuple('D',0,0));
        }
        void set_capacity(uint cluster_n){//,uint edgesPerCluster_n
            capacityClusters = cluster_n;
            LRU_IMCluster.set_capacity(capacityClusters);
        }
        void clear(){
            LRU_IMCluster.clear();
            HotPools.clear();
            clusterStatus.clear();
            MCEdges_EM.clear();
//            im_partitions.clear();
        }
        ~MCHotPool(){}
    };
    //hot pool struct for OneHopNoIO
    struct MCHotPool2 {
        uint capacityClusters = Partition_N;//Capacity of how many partitions can be store in memory, if one partition size is 1MBB, then 128 partitions equal 128MB
        LRUCache LRU_IMCluster;//LRU manager for clusters
        vector<vector<MCEdgeCSR>> HotPools;  //storing the info of partitions in memory or temp file. <partition_id,<adjacency lists>>
        vector<bool> FlagIM;    //the flag of whether the partition is in memory

        MCHotPool2(){};
        MCHotPool2(uint pn,uint capacity){
            HotPools.assign(pn,vector<MCEdgeCSR>());
            FlagIM.assign(pn, false);
            set_capacity(capacity);
        };
        ~MCHotPool2(){clear();};
        void set_capacity(uint cluster_n){//,uint edgesPerCluster_n
            capacityClusters = cluster_n;
            LRU_IMCluster.set_capacity(capacityClusters);
        }
        void clear(){
            LRU_IMCluster.clear();
            HotPools.clear();
            FlagIM.clear();
        }
    };
    /// Class of external multi-criteria graph
    class EMMCGraph{
    public:
        /// basic variables
        uint memSize = 2048;    //Overall memory size, unit(MB)
        uint memGraph = 512;    //memory size for graph data, unit(MB)
        uint node_num = 0;      //the number of vertices
        uint edge_num = 0;      //the number of edges
        uint num_criteria = 0;   //the number of criteria
        usint sc_i = 0;         //the index of criterion in processing for multi-pass algorithms
        string dataset;         //dataset name
        vector<string> mc_criteria; //the criteria(type) of graph
        string QueryType;       //query type
        int max_degree;     //the maximal degree of the graph
        double ave_degree;  //the average degree of the graph
        vector<int> max_w;  //the maximal edge weight of different criteria
        list<int> list_sp;      //list to record nodes in the shortest path
//        vector<vector<MCEdge>> MCEdges; //adjacent mc-edges list for in-memory algorithms

        /// variables for performance
        uint visited_number = 0;//the number of visited vertices
        uint io_num = 0;    // the total number of I/O read
        uint p_num = 0;      //the number of partitions in memory
        int share_number = 0;//the number of shared visiting
        double query_time = 0;  //query time, ms
        double io_time = 0; //the total time of io, ms
//        unordered_set<int> set_visited;     //set of recording the visited vertices of all criteria
        unordered_set<int> set_readCluster; //the partition been read, used to compute the minimal required partition
        vector<Distance> mc_min_cost;//the min cost of each criterion
        vector<bool> mc_finished;   //whether this criterion finished

        /// variables for algorithms
        set<int> set_cri; //set to record whether the SP has been found in the criteria
        // Variables for EM_Dijk
        vector<usint> num_pool;
        vector<vector<bitset1<WeightPowMax+1,EdgeWeight>>> node_to_category;//vector for recording the edge category for vertices, the first bit is for done bit. In order to avoid repeat graph reading, we use vector to store node_to_category information for all criteria, which may increase the memory consumption by (num_of_cri-1)*|V|.
//        vector<vector<bitset<WeightPowMax+1>>> node_to_category;
        vector<vector<EdgeWeight>> pow_low;
        vector<tuple<bool,int,int>> EMEdgesIndex; //the state and index (im vector or external vector), <whether_to_evict,index_begin,index_end>
        vector<tuple<bool,bool,int,int>> EMEdgesIndex_Bi;   //the state and index for BiMultiHops algorithm
        // Variables for One-Pass algorithms
        vector<PriorityQueue*> EM_MC_PQueue;
        vector<PriorityQueue2*> EM_MC_PQueue_r;
        vector<int8_t> vertex_cri; //record the number that one criterion is processed
        vector<pair<int8_t,int8_t>> vertex_cri_Bi;

        /// variables for partition-based IO management
        uint partition_number = 0;       //the number of partitions

        vector<int> partitions_info;    //vector for storing the info of partitions, start vertex id
        vector<usint> node_to_cluster;//used to store the partition id of node
        vector<int> cluster_to_IO;  //map of recording visited partitions and corresponding #IO
        vector<bool> partitions_read;   //flag of whether a partition has been read
        int partition_left = 0;
        double alpha = 0.9;//1;//   // threshold of One-Hop algorithm
        double alpha_multi = 0.9;//1;// // threshold of Multi-Hop algorithm
        double alpha_bi = 0.6;//0.9;//  // threshold of BiMulti-Hops algorithm


        /// Common Functions
//        void MCReadGraph_MCEdges(string filename);  //one-off edges reading for MC-graph
        void ClusterInfoLoad(bool ifMap);     //function for reading partition information
        void CommonInitiation_IO();     //Initiation function for each round of IO-efficient algorithms

        void MC_Evaluate_IO(const string& qtype);  //Entry for io-efficient algorithm

        template <class T>
        void GraphReadCluster4(HotPool4<T> & mcHotPool,int target_p_id);//function of read partition for HotPool4

        /*** For One-Pass algorithms ***/
        double MC_OnePass(const string& qtype,int algo_choice);//Entry for one-pass algorithms
        double MC_OnePass_one(const string& filename, const string& qtype, int algo_choice);//Function for calling different types of one-pass algorithms
        void OnePassInitiation_IO(int algo_choice); //function of initiating One-Pass algorithms
        void OnePassClear_IO(int algo_choice);

        /// One-hop without io optimization
        void EM_MC_OneHop_NoIO_new(int node_start, int node_end);   //One-hop algorithm without io optimization

        /// Multi-hops without io optimization
        void EM_MC_MultiHop_NoIO_new(int node_start, int node_end);//Function for multi-hops algorithm with io optimization

        /// BiMulti-hops without io optimization
        void EM_MC_BiMultiHop_NoIO_new(int node_start, int node_end);

        /// One-hop with io optimization
        void EM_MC_OneHop(int node_start, int node_end);//One-hop algorithm with io optimization
        void GraphMCReadCluster_New(MCHotPool<VectorMCEdgesEMTuple_IO> & mcHotPool,NodeId ID1,int target_p_id,int algo_choice);//function of reading partition for one-pass algorithms with io optimization

        /// Multi-hops with io optimization
        void EM_MC_MultiHop(int node_start, int node_end);//Function of multi-stop algorithm with io optimization

        /// BiMulti-hops with io optimization
        void EM_MC_BiMultiHop(int node_start, int node_end);
        void GraphMCReadCluster_New_Bi(MCHotPool<VectorMCEdgesEMTuple_IO_Bi> & mcHotPool,int target_p_id,vector<int8_t>& cluster_to_bi,int8_t direction,vector<bool>& clusterRead);//function of reading partition for one-pass algorithms with io optimization NodeId& ID1,
//        void GraphMCReadCluster_New_Bi(MCHotPool & mcHotPool,NodeId& ID1,int target_p_id,vector<int8_t>& cluster_to_bi, bool flag_reverse);//function of reading partition for one-pass algorithms with io optimization

        static int Weight_to_category(int w);
        void EM_IORecord(Stats& a);
        int Minimal_IO();
        bool EM_JudgeEmpty(vector<PriorityQueue*> &em_mc_pqueue);
        bool EM_JudgeEmptyBi(vector<PriorityQueue*> &em_mc_pqueue, vector<PriorityQueue2*> &em_mc_pqueue_r);
        void InfoPrint() const;//function of printing basic information
        void MemoryCheck(int algo_choice, uint mem, bool &flag_exit);

        list<int> Dij_getPath(vector<NodeId> & pre, NodeId ID1,NodeId ID2);
        list<int> BiDij_getPath(vector<NodeId> & pre,vector<NodeId> & pre_b,NodeId ID1,NodeId ID2,NodeId terminate_id);
    };

}

#endif //MCSPS_EMGRAPH_H
