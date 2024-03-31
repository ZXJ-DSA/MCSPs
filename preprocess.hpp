/*
 * Filename:    preprocess.hpp
 * Description: functions of graph data preprocessing
 * Created:     05 March 2022
 * Authors:     Xinjie ZHOU
 */

#ifndef MCSPS_PREPROCESS_HPP
#define MCSPS_PREPROCESS_HPP
#include "preprocess.h"

namespace gbpre{
    //Function for initiating the class variables
    /*---Data load---*/
    //Function for general reading and print the first few lines
    void Graph_pre::Read(const string& file, int num_line)
    {
        int num_i=0;
        string lineStr;//data of line
        string temp;
        Timer tt;
        tt.start();
        //Open file
        ifstream inFile(file, ios::in);
        if (!inFile) {
            cout << "File opening failed." << endl;
            assert(inFile);
            exit(1);
        }
        cout << "Data loading..." << endl;
        while (getline(inFile, lineStr)) {//num_i < num_line
            cout << lineStr << endl;
            if(num_line > 0){
                ++num_i;
                if(num_i>num_line)
                    break;
            }
        }
        cout << "Data loaded." << endl;
        inFile.close();
        tt.stop();
        cout << "The time used for file loading:" << tt.GetRuntime() << " seconds" << endl;
    }
    //Function for reading graph edges data
    void Graph_pre::ReadGraph_edges(const string& filename,bool ifBase)
    {
        Timer tt;
        tt.start();
        string line_symbol;
        string temp_str;
        int num_line = 0;
        int num_show = 6;//number of lines to be printed
        int ID1, ID2, weight;

        //Open file
        ifstream inFile(filename, ios::in);
        if (!inFile) {
            cout << "File opening failed." << endl;
            assert(inFile);
            exit(1);
        }
        cout << dataset<<" "<< criteria<<" graph edges Data loading..." << endl;
        while (inFile) {//read each line to lineStr
            inFile >> line_symbol;
            if (line_symbol == "p") {//read graph basic information
                inFile >> temp_str >> this->node_num >> this->edge_num;
//                this->Nodes.assign(node_num,Node());//assign the space for graph
                EdgesMap.assign(node_num,map<NodeId,Edge>());
                cout << "Some edges:" << endl;
                //break;
            }
            else if (line_symbol == "a") {//read graph data
                inFile >> ID1 >> ID2 >> weight;
//                if(invalid_edges.find(make_pair(ID1,ID2))!=invalid_edges.end()){//if found in invalid set
//                    continue;
//                }
                if(ifBase){
                    base_max = max(base_max,weight);
                    base_min = min(base_min,weight);
                }
//                EdgesMap[ID1 - 1].insert({ID2-1,Edge(ID2 - 1, weight)});
                EdgesMap[ID1].insert({ID2,Edge(ID2, weight)});
                if (num_line < num_show) {
                    cout << ID1 << "\t" << ID2 << "\t" << weight << endl;
                    ++num_line;
                }
                else if (num_line == num_show) {
                    cout << "..." << endl;
                    ++num_line;
                }

            }else{
                getline(inFile, temp_str);
            }
        }
        cout << "Data loaded. ";
        inFile.close();
        tt.stop();
        cout << "The number used for data loading: " << tt.GetRuntime() << " s." << endl;
        cout << "Number of nodes:" << node_num << endl;
        cout << "Number of edges:" << edge_num << endl;
        if(ifBase){
            cout << "The maximal edge weight is " << base_max << " , while the minimal edge weight is "<< base_min << endl;
        }
        cout << "--------------------" << endl;
    }

    //Function for extra criteria data generation
    void Graph_pre::MC_DataGeneration(int num){//memory usage: 7*|E|
        /*** Extra criteria data generation ***/
        string r_Edge = string(DataPath)  + dataset + "/" + dataset + "_" + criteria + ".txt";
        string r_Edge_t = string(DataPath)  + dataset + "/" + dataset + "_Time.txt";
//        Read_MCEdgesMap(r_Edge);
        Read_ToMCEdgesMap(r_Edge,true);
        if(num_criteria == 2){
            Read_ToMCEdgesMap(r_Edge_t,false);
        }
        //Generate extra criteria data
        MC_Generate(num);
        //Write final MCEdges to disk, in edge-pair style
//        MC_WriteEdges(false,true);
        //Write final MCEdges to disk, in io style
//        MC_WriteEdgesIO(false);
    }
    //Function for data generation and writing, new version
    void Graph_pre::MC_Generate(int num){
        string w_file = DataPath + dataset + "/" + dataset + "_MCEdges.txt";
        int base;
        int temp_num, temp_int;
        int ID1,ID2;
        //normal distribution generator
        default_random_engine generator(rand()/1024);
        normal_distribution<double> norm_distribution(0.0, 0.2);
        cout <<"Extra criteria data generating..."<<endl;
        mc_criteria.emplace_back("Distance");
        if(num_criteria == 1){
            mc_criteria.emplace_back("Corr");
        }else if(num_criteria == 2){
            mc_criteria.emplace_back("Time");
        }
        mc_criteria.emplace_back("Corr");
//        mc_criteria.emplace_back("Corr");
        mc_criteria.emplace_back("Anti");
        mc_criteria.emplace_back("Indep");

        for(ID1=0; ID1<MCEdgesMap.size(); ID1++) {//for each node
            temp_num = 0;
            for(auto it=MCEdgesMap[ID1].begin();it!=MCEdgesMap[ID1].end();++it){//for each adjacent node
                ID2 = it->first;
                base = it->second.w[0];
                if(num_criteria == 1) {
                    //correlated data
                    temp_int = max(int((0.5 + norm_distribution(generator)) * base), base_min);
                    MCEdgesMap[ID1][ID2].w[1] = temp_int;
                    MCEdgesMap[ID2][ID1].w[1] = temp_int;
                    max_weights[1] = max(max_weights[1],temp_int);
                    min_weights[1] = min(min_weights[1],temp_int);
                }
                //correlated data
                temp_int = max(int((0.7 + norm_distribution(generator)) * base), base_min);
                MCEdgesMap[ID1][ID2].w[2] = temp_int;
                MCEdgesMap[ID2][ID1].w[2] = temp_int;
                max_weights[2] = max(max_weights[2],temp_int);
                min_weights[2] = min(min_weights[2],temp_int);
                //correlated data
//                temp_int = max(int((0.9 + norm_distribution(generator)) * base), base_min);
//                MCEdgesMap[ID1][ID2].w[3] = temp_int;
//                MCEdgesMap[ID2][ID1].w[3] = temp_int;
//                max_weights[3] = max(max_weights[3],temp_int);
//                min_weights[3] = min(min_weights[3],temp_int);
                //anti-correlated data
                temp_int = max(int(1000*(base_max - (1+norm_distribution(generator))*base) / base_max), base_min);
                MCEdgesMap[ID1][ID2].w[3] = temp_int;
                MCEdgesMap[ID2][ID1].w[3] = temp_int;
                max_weights[3] = max(max_weights[3],temp_int);
                min_weights[3] = min(min_weights[3],temp_int);
                //independent data, from base_min to 1000
                temp_int = base_min + int(1000*(rand()/double(RAND_MAX)) * (base_max - base_min) / base_max);
                MCEdgesMap[ID1][ID2].w[4] = temp_int;
                MCEdgesMap[ID2][ID1].w[4] = temp_int;
                max_weights[4] = max(max_weights[4],temp_int);
                min_weights[4] = min(min_weights[4],temp_int);
            }
        }
        num_criteria = 5;//update the number of criteria
        cout <<"Done."<<endl;
    }
    //Function for converting MCEdge data to MCEdgesS.txt
    void Graph_pre::MC_ToEdgesS() {
        string r_file = string(DataPath) + dataset + "/" + dataset + "_MCEdges_" + to_string(PartitionSize) + ".txt";
        string w_file = string(DataPath) + dataset + "/" + dataset + "_MCEdgesS_" + to_string(PartitionSize) + ".txt";

        Timer tt;
        tt.start();
        string lineStr;//data of line
        string line_symbol;
        string temp_str;
        int ID1, ID2, weight;
        int num_cri;

        //Open file
        ifstream inFile(r_file, ios::in);
        if (!inFile) {
            cout << "File opening failed." << endl;
            assert(inFile);
            exit(1);
        }
        ofstream outFile(w_file, ios::out);
        if (!outFile) {
            cout << "File opening failed." << endl;
            assert(outFile);
            exit(1);
        }
        cout << "Graph data converting..." <<endl;
        while (getline(inFile, lineStr)) {//read each line to lineStr
            if (lineStr != "") {
                istringstream sin(lineStr); //read each string to sin
                sin >> line_symbol;
                if (line_symbol == "p") {//read graph basic information
                    sin >> temp_str;
                    if(temp_str == "sp"){
                        sin >> node_num >> edge_num;
                        outFile << node_num <<" "<<edge_num<<endl;

                    }
                    else if(temp_str == "criteria"){
                        sin >> num_cri;
                        num_criteria = num_cri;
                        outFile<<num_cri;
                        for(int i=0;i<num_cri;++i){
                            sin >> temp_str;
                            outFile<<" "<<temp_str;
                        }
                        outFile<<endl;
                    }
                }
                else if (line_symbol == "a") {//read graph data
                    sin >> ID1 >> ID2;
                    outFile <<ID1<<" "<< ID2;

                    for(int i=0;i<num_cri;++i) {//get edge weights
                        sin >> weight;
                        outFile <<" "<< weight;
                    }
                    outFile << endl;
                }
            }
        }
        inFile.close();
        outFile.close();
        tt.stop();
        cout << "The time used for data converting:" << tt.GetRuntime() << " s." << endl;
        remove(r_file.c_str());
    }
    //Function of deterministic graph clustering, memory usage: about 4*|V|+5*|E|
    bool Graph_pre::DeterClustering(bool ifReuse){
        EdgeSorter edgeSorter(EdgePairComparator(),MEMORY_FOR_SORTER);//edge id sorter for MST
        VectorUint EulerPath;//external vector for recording the Euler tour path
//        list<NodeId> EulerPath;
        vector<uint> outDegree(node_num,0); //for undirected Euler tour, we only need out degree
        vector<uint> indexID(node_num,0);   //record the index ID of vertices in mst vector
        vector<int> vertexRank(node_num,-1);//vertex rank result
//        unordered_map<NodeId,NodeId> remap; //vertex id map after clustering
        uint mstEdgeNum=0;      //mst edge number
        bool ifSame = true;

        //statstics of IO
        Stats stats_b, stats_e;
        stats_b.start();
        int node_start = 0;//rand()%node_num;//source node for mst

        /// Computing MST
        std::cout<<"Starting MST computing.."<<std::endl;
        sc_i=0; //criterion for partitioning
        if(strategy == "pri"){
            cout<<"The primary criterion is: "<<sc_i<<endl;
        }
        bool mstExist = MinSpanTree(node_start,edgeSorter,outDegree);//complexity is O(|E|log(|E|))
        if(!mstExist){
            cout<<"Wrong for MST!"<<endl; exit(1);
        }
        edgeSorter.sort();//sort the edges in MST
        mstEdgeNum = edgeSorter.size();
//        VectorEdgePairs mstEdges2(mstEdgeNum);//vector of storing the edges of MST
//        vector<list<NodeId>> mstEdges(node_num);//vector of storing the edges of MST
        vector<EdgePair> mstEdges(mstEdgeNum);//MST edges
        assert(mstEdgeNum==2*node_num-2);
        stxxl::stream::materialize(edgeSorter,mstEdges.begin(),mstEdges.end());//store mst edges
        edgeSorter.clear();
        std::cout<<"MST computing is finished."<<std::endl;
        cout<<"MST edges: "<<endl;
        for(int j=0;j<5;++j){//50
            cout<<mstEdges[j].ID1<<" "<<mstEdges[j].ID2<<endl;
        }
        //compute the index id
        indexID[0] = 0;
        for(int i=1;i<node_num;++i){
            indexID[i]=indexID[i-1]+outDegree[i-1];
        }
        /// Euler tour computing
        //determine better node start
//        for(int i=0;i<node_num;++i){
//            if(outDegree[i]==1){
//                node_start = i;
//                break;
//            }
//        }
        node_start = 0;
        DFS_EulerTour(node_start,mstEdges,indexID,outDegree,EulerPath);//complexity is O(|E|)
        if(EulerPath.size() == mstEdgeNum+1){
            cout<<"Euler path found."<<endl;
        }else{
            cout<<"No Euler path found."<<endl;
        }
        /// List ranking
        uint i=0;
        cout<<"Euler path: "<<endl;
        for(auto it=EulerPath.begin();it!=EulerPath.end();++it){
            assert(*it < node_num);
            if(vertexRank[*it] == -1){
                vertexRank[*it] = i;
            }
            if(i<20){
                cout<<*it<<" ";
            }
            ++i;
        }
        cout<<endl;
        for(i=0;i<node_num;++i){
            if(vertexRank[i]==-1){
                cout<<"Wrong! vertexRank["<<i<<"] is "<<vertexRank[i]<<endl;
                exit(1);
            }
            rankSorter.push(make_pair(i,vertexRank[i]));//record the offset of each vertex in euler path
        }
        vertexRank.clear();
        rankSorter.sort();//sort
        // partition result for visualization
        /*uint temp_size=0;
        vector<vector<NodeId>> partitions;
        vector<NodeId> partitionNodes;
        while(!rankSorter.empty()){//if not full
            if(temp_size + MCEdgesMap[rankSorter->first].size() < PartitionMCEdge_NUM){
//                cout<<rankSorter.size()<<endl;
                partitionNodes.push_back(rankSorter->first);
                temp_size += MCEdgesMap[rankSorter->first].size();
                ++rankSorter;
            }else{
                partitions.push_back(partitionNodes);
                partitionNodes.clear();
                temp_size=0;
            }
        }
        partitions.push_back(partitionNodes);
        partitionNodes.clear();
        string file_w_p = "/Users/zhouxj/Python/MCSP/data/" + this->dataset + "_Partitions.txt";
        ofstream outFile(file_w_p, ios::out);
        if (!outFile) {
            cout << "Write File opening failed." << endl;
            assert(outFile);
        }
        outFile << "p all "<<partitions.size()<<endl;
        for(uint j=0;j<partitions.size();++j){
            outFile << "p partition "<<j<<endl;
            for(auto it=partitions[j].begin();it!=partitions[j].end();++it){
                outFile << "v "<<*it<<endl;
            }
        }
        outFile.close();*/
        /// Chopping into partitions
//        i=0;
//        while(!rankSorter.empty()){//map from old vertex id to the new one
//            nodes_map[rankSorter->first] = i;
//            if(ifSame){
//                if(rankSorter->first != i){//first is the vertex id; second is the rank of this vertex to the start of the Euler tour
//                    ifSame = false;
////                    cout << rankSorter->first << " "<<i<<endl;
//                }
//            }
//            ++rankSorter; ++i;
//        }
//        if(ifReuse){
//            string file_remap = string(DataPath) + this->dataset + "/" + this->dataset + "_Partition_remap_"+strategy+".txt";
//            ofstream outFile(file_remap, ios::out);
//            if (!outFile) {
//                cout << "Write File opening failed." << endl;
//                assert(outFile);
//                exit(1);
//            }
//            outFile<<"old_id new_id"<<endl;
//            for(int j=0;j<node_num;++j){
//                outFile<<j<<" "<<nodes_map[j]<<endl;//+1
//            }
//        }
        if(ifSame){
            cout << "There is no vertices remapping." << endl;
        }else{
            cout << "The id of vertices have been remapped!"<< endl;
        }
//        rankSorter.rewind();
        cout<<"Rank sorter size: "<<rankSorter.size()<<endl;
        cout<<"Sorted by vertex rank: "<<endl;
//        exit(0);
        WriteClusters(rankSorter);//generate clusters
        //io record
        stats_e = stats_b.get_stats();
        if (stats_e.m_noOfWrites > 0 || stats_e.m_noOfReads > 0) {
            cout << "There are " << stats_e.m_noOfWrites + stats_e.m_noOfReads << " io for stxxl." << endl;
            //cout<<stats_e<<endl;
        }
//        rankSorter.clear();
        return ifSame;
    }


    //Function of counting the minimal spanning tree, stl::priority queue
    bool Graph_pre::MinSpanTree(int node_start,EdgeSorter & mstEdges,vector<uint> & outDegree)//Lazy Prim's algorithm, complexity O(|E|log(|E|)). Eager Prim's algorithm use an indexed priority queue which can efficiently update and poll key-value pairs, further reduce the complexity to O(|E|log|V|).
    {
//        benchmark::heap<2, int, int> pqueue(node_num);
        //statstics of IO
        Stats stats_b, stats_e;
        stats_b.start();
//        PriorityQueueEdgePair pqueue(PQEdgePair_Pool);
        priority_queue<pair<uint,EdgePair>, vector<pair<uint,EdgePair>>, PQEdgePairCompareLess> pqueue;
        unsigned long long mstCost=0;
        pair<uint,EdgePair> item_;
        int item_id,item_dis;
        int temp_id,temp_dis;
        Distance temp_w=0;
        vector<bool> closed(node_num, false); //flag vector of whether closed

        uint mst_edgenum = node_num-1;//the number of undirected edges
        uint edgeCount = 0;
        double ave_w=0;
        unsigned long long searchFrontier=0;

        int aggregate=0;
        if(strategy == "ave"){
            aggregate=1;
            cout<<"Edge aggregation strategy: "<<strategy<<endl;
        }else if(strategy == "pri"){
            aggregate=2;
            cout<<"Edge aggregation strategy: "<<strategy<<endl;
        }else if(strategy == "wave"){
            aggregate=3;
            cout<<"Edge aggregation strategy: "<<strategy<<endl;
        }else{
            cout<<"Wrong edge aggregation strategy! "<<strategy<<endl;
            exit(1);
        }

        //Initiation of start node
        item_id = node_start;
        closed[item_id] = true;
        for (auto it = MCEdges[item_id].begin(); it != MCEdges[item_id].end(); ++it) {
            temp_id = it->u;
            if(closed[temp_id])
                continue;
            if(aggregate == 1){
                temp_w = ceil(accumulate(it->w.begin(),it->w.end(),0.0)/it->w.size());
            }else if(aggregate == 2){
                temp_w = it->w[sc_i];
            }else if(aggregate == 3){
                ave_w = 0;
                for(int i=0;i<num_criteria;++i){
                    ave_w += zeta[i] * it->w[i];
                }
                temp_w = ceil(ave_w);
            }else{
                cout<<"Wrong edge aggregation strategy! "<<strategy<<endl;
                exit(1);
            }
            pqueue.push(make_pair(temp_w,EdgePair(item_id,temp_id)));
        }

        //Iteration
        while (!pqueue.empty() && edgeCount!=mst_edgenum) {//for every node in pqueue
            if(searchFrontier<pqueue.size()){
                searchFrontier=pqueue.size();
            }
            item_ = pqueue.top();// top min item
            pqueue.pop();
            item_id = item_.second.ID2;
            if(closed[item_id])
                continue;
            //push the minimal edges into mst
            mstEdges.push(item_.second);
            mstEdges.push(EdgePair(item_.second.ID2,item_.second.ID1));
            //count the out degree too
            ++outDegree[item_.second.ID1];
            ++outDegree[item_.second.ID2];

            ++edgeCount;
            mstCost += temp_w;
            closed[item_id] = true;
            for (auto it = MCEdges[item_id].begin(); it != MCEdges[item_id].end(); ++it) {
                temp_id = it->u;
                if(closed[temp_id])
                    continue;
                if(aggregate == 1){
                    temp_w = ceil(accumulate(it->w.begin(),it->w.end(),0.0)/it->w.size());
                }else if(aggregate == 2){
                    temp_w = it->w[sc_i];
                }else if(aggregate == 3){
                    ave_w = 0;
                    for(int i=0;i<num_criteria;++i){
                        ave_w += zeta[i] * it->w[i];
                    }
                    temp_w = ceil(ave_w);
                }else{
                    cout<<"Wrong edge aggregation strategy! "<<strategy<<endl;
                    exit(1);
                }
                pqueue.push(make_pair(temp_w,EdgePair(item_id,temp_id)));
            }
        }
        cout<<"Size of spanning tree: "<<mstEdges.size()<<endl;
        cout<<"Largest search frontier size: "<<searchFrontier<<endl;
//        cout<<"MST cost is "<<mstCost<<endl;
        ///io record
        stats_e = stats_b.get_stats();
        if (stats_e.m_noOfWrites > 0 || stats_e.m_noOfReads > 0) {
            cout << "There are " << stats_e.m_noOfWrites + stats_e.m_noOfReads << " io for stxxl." << endl;
            //cout<<stats_e<<endl;
        }
        if(edgeCount!=mst_edgenum){
            cout<<"!!! MST does not exist! edgeCount: "<<edgeCount<<", mst_edgenum: "<<mst_edgenum<<endl;
//            exit(1);
            return false;
        } else{
            cout<<"Find MST. "<<edgeCount<<", mst_edgenum: "<<mst_edgenum<<endl;
            return true;
        }
    }

    //Function of DFS for Euler Tour
    void Graph_pre::DFS_EulerTour(NodeId node_start,vector<EdgePair> & mstEdges,vector<uint>& indexID,vector<uint> & outDegree,VectorUint& EulerPath)//list<NodeId>& EulerPath
    {
        stxxl::STACK_GENERATOR<NodeId>::result stack_A;//external stack
//        stack<uint> stack_A;
        NodeId item_id,temp_id;
        stack_A.push(node_start);
        int pre = -1;
        bool flag_add = false;
        while(!stack_A.empty()){
            item_id = stack_A.top();
            stack_A.pop();
            flag_add = false;
            while(outDegree[item_id]!=0) {
                //Select the next unvisited outgoing edge
                if(pre!=-1 && !flag_add){
                    stack_A.push(pre);
                    flag_add=true;
                }
                --outDegree[item_id];
                temp_id = mstEdges[indexID[item_id] + outDegree[item_id]].ID2;
                if(temp_id!=pre)
                    stack_A.push(temp_id);
            }
            pre = item_id;
            EulerPath.push_back(item_id);
//            cout<<item_id<<endl;
        }
    }
    void Graph_pre::WriteClusters(RankSorter& rankSorter)//unordered_map<NodeId,NodeId>& remap
    {
        vector<NodeId> partitionNodes;
        uint temp_size=0;
        uint index=0;
        int cluster_i=0;

//        partitions_info.clear();
//        partitions_info.emplace_back(0);
        cluster_to_node.clear();
        ///update edge number for each partition
//        PartitionMCEdge_NUM = edge_num / (PartitionNumber-1);//partition by overall partition number
        PartitionMCEdge_NUM = PartitionSize*1024 / MCEdge_SZ;//partition by average partition disk size
//        uint PartitionIntNum = PartitionSize*1024/sizeof(int);
        PartitionNumber = edge_num/PartitionMCEdge_NUM + 1;
        if(PartitionMCEdge_NUM < max_degree){
            cout << "!!!The partition size is so small that there will be partitions that exceed the average partition size." <<  endl;
        }
        cout << "!!The partition number is about " << PartitionNumber << " . The disk size of each partition is about " << (PartitionMCEdge_NUM*MCEdge_SZ/1024 + 1) << " KB." <<endl;
        cout << "Writing partitions into file..." <<endl;

        while ( !rankSorter.empty() ){
            bool ifExceed=false;
            while(!rankSorter.empty()){//if not full
//                if(cluster_i>=4734){
//                    cout<<"Flag 2 "<<cluster_i<<" "<<temp_size<<" "<<rankSorter->first<<endl;
//                }
                if(temp_size + MCEdges[rankSorter->first].size() <= PartitionMCEdge_NUM){//if the partition is not full after adding this vertex
//                    if(cluster_i>=4734){
//                        cout<<"Flag 3.1 "<<cluster_i<<" "<<MCEdges[rankSorter->first].size()<<endl;
//                    }
//                    cout<<rankSorter.size()<<endl;
                    partitionNodes.emplace_back(rankSorter->first);
                    temp_size += MCEdges[rankSorter->first].size();
                    if(index<5) {
                        cout<<rankSorter->first<<" "<<rankSorter->second<<endl;
                    }
                    ++rankSorter; ++index;
                }else{//if the partition is full after adding this vertex
//                    if(cluster_i>=4734){
//                        cout<<"Flag 3.2 "<<cluster_i<<" "<<MCEdges[rankSorter->first].size()<<" "<<partitionNodes.size()<<endl;
//                    }
                    if(partitionNodes.empty()){//if empty, which means the adjacency list of one single vertex is larger than file size limitation
                        cout << "Partition " << cluster_i << " exceed the average partition size! " << rankSorter.size()<< " "<<index<<endl;
                        partitionNodes.emplace_back(rankSorter->first);
                        ++rankSorter; ++index;
                        ifExceed = true;
                    }
                    break;
                }
            }
//            if(cluster_i>=4734){
//                cout<<"Flag 4 "<<cluster_i<<" "<<partitionNodes.size()<< endl;
//            }
//            partitions_info.emplace_back(index);
            cluster_to_node.push_back(partitionNodes);
            //if one partition is full write it to disk
            // bin file
            WritePartition(cluster_i, partitionNodes,  true);

            // txt file
            if(ifExceed){
                 WritePartition(cluster_i, partitionNodes, false);
            }

//            if(cluster_i>=4734){
//                cout<<"Flag 1 "<<cluster_i<<" "<<cluster_to_node.size()<<endl;
//            }

            partitionNodes.clear();
            ++cluster_i;
            temp_size = 0;
        }

        if(cluster_to_node.size() != PartitionNumber){
            cout << "!!!The actual number of partitions is " << cluster_to_node.size() << endl;
            PartitionNumber = cluster_to_node.size();
        }else{
            cout << "!!!The partition number is " << PartitionNumber  <<endl;
        }
        //write partition_info to disk
        WritePartitionInfo(cluster_to_node);
        cout << "Done." << endl;

        // write the remapped graph into disk
//        string w_file = DataPath + dataset + "/" + strategy +"/" +dataset + "_"+ to_string(PartitionSize) +".MCEdges";
//        int base;
//        int temp_num, temp_int;
//        int ID1,ID2;
//
//        ofstream outFile(w_file,ios::out);
//        if(!outFile.is_open()){//if exists
//            cout<<"Cannot open file "<<w_file<<endl;
//            exit(1);
//        }
//        outFile << node_num <<" "<<edge_num<< endl;
//        outFile<<NUM_OF_CRITERIA;
//
//        for (int i = 0; i < mc_criteria.size(); ++i) {
//            if(i == mc_criteria.size()-1)
//                outFile<<" "<<mc_criteria[i]<<endl;
//            else
//                outFile<<" "<<mc_criteria[i];
//        }
//
//        bool flag_double = true;
//        if(edge_num>1000000000){
//            flag_double=false;
//            cout<<"Write edge once!"<<endl;
//        }
//
//        for(int ID=0;ID<node_num;++ID){
//            for(auto it2=MCEdges[ID].begin();it2!=MCEdges[ID].end();++it2){
//                if(ID<it2->u){
//                    outFile<<ID<<" "<<it2->u;
//                    for(int i=0;i<it2->w.size();++i){
//                        outFile<<" "<<it2->w[i];
//                    }
//                    outFile<<endl;
//
//                    if(flag_double){
//                        outFile<<it2->u<<" "<<ID;
//                        for(int i=0;i<it2->w.size();++i){
//                            outFile<<" "<<it2->w[i];
//                        }
//                        outFile<<endl;
//                    }
//                }
//
//            }
//        }
//        outFile.close();
    }

    void Graph_pre::WritePartitionInfo(vector<vector<NodeId>>& partitionNodes){
        //write partition_info to disk
        string file_w_pf = string(DataPath) + this->dataset + "/" + strategy + "/" + this->dataset + "_Partitions_" +partMethod+"_"+ to_string(PartitionSize) +"_Info.txt";
        ofstream outFile2(file_w_pf, ios::out);
        if (!outFile2) {
            cout << "Write File opening failed." << endl;
            assert(outFile2);
            exit(1);
        }
        cout << "Writing partitions information into file..." << endl;
        outFile2 << node_num << " " <<edge_num <<endl;
        for (int i = 0; i < num_criteria; ++i) {
            if(i == num_criteria-1)
                outFile2<<mc_criteria[i]<<endl;
            else
                outFile2<<mc_criteria[i]<<" ";
        }
        outFile2 << partitionNodes.size() <<endl;
        for(int pid=0;pid<partitionNodes.size();++pid){
            outFile2 << partitionNodes[pid].size();
            for(auto it=partitionNodes[pid].begin();it!=partitionNodes[pid].end();++it){
                outFile2 << " " << *it;
            }
            outFile2 << endl;
        }
        outFile2.close();
    }

    void fwriteInt(int& x, int& read_number, int * write_buff, FILE * bfile){

        write_buff[read_number++]=x;
        if(read_number == VERTEX_PER_BLK){
            fwrite(write_buff,SZ_VERTEX,read_number,bfile);
            read_number = 0;
        }
    }

    void Graph_pre::WritePartition(int cluster_i, vector<NodeId>& partitionNodes, bool ifbin) {//map<NodeId,NodeId>& remap,
//        cout<<"Writing cluster "<<cluster_i<<endl;
//        if(cluster_i>=4734){
//            cout<<"Writing "<<cluster_i<<endl;
//        }
        if(ifbin){
            char file_bin_path[300];

            strcpy(file_bin_path, DataPath);
            strcat(file_bin_path, dataset.c_str());
            strcat(file_bin_path, "/");
            strcat(file_bin_path, strategy.c_str());
            strcat(file_bin_path, "/Partitions/");
            strcat(file_bin_path, dataset.c_str());
            strcat(file_bin_path, "_Partition_");
            strcat(file_bin_path, partMethod.c_str());
            strcat(file_bin_path, "_");
            strcat(file_bin_path, to_string(PartitionSize).c_str());
            strcat(file_bin_path, "_");
            strcat(file_bin_path, to_string(cluster_i).c_str());
            strcat(file_bin_path, ".bin");

            int * write_buff = (int *)malloc(sizeof(int)*VERTEX_PER_BLK);//4MB=4*1024*1024=4*1048576
            FILE * bfile = fopen(file_bin_path, "wb");
            if(!bfile){
                cout<<"Failed to open file "<<file_bin_path<<endl;
                exit(1);
            }
            int x;
            int read_number=0;
            int degree=0;

//            cout<<"Flag 2"<<endl;

            x=cluster_i;
            fwriteInt(x,read_number,write_buff,bfile);

//            cout<<"Flag 3"<<endl;

            x=partitionNodes.size();
            fwriteInt(x,read_number,write_buff,bfile);


            int item_id;
            for (uint j = 0; j < partitionNodes.size(); ++j) {
                item_id = partitionNodes[j];

                x=item_id;
                fwriteInt(x,read_number,write_buff,bfile);

//                cout<<"Flag 5 "<<item_id<<" "<<j<<" "<<partitionNodes.size()<<endl;
                x=MCEdges[item_id].size();
                fwriteInt(x,read_number,write_buff,bfile);

                for(auto it=MCEdges[item_id].begin();it!=MCEdges[item_id].end();++it){
                    x=it->u;
                    fwriteInt(x,read_number,write_buff,bfile);

                    for (auto it2 = it->w.begin(); it2 != it->w.end(); it2++) {
                        x=*it2;
                        fwriteInt(x,read_number,write_buff,bfile);
                    }
                }
            }

//            cout<<"Flag 6 "<<read_number<<endl;
            if(read_number>0){
                fwrite(write_buff,SZ_VERTEX,read_number,bfile);
                read_number = 0;
            }

            fclose(bfile);
            free(write_buff);
        }
        else{
            // txt file
            string file_w_p = string(DataPath) + this->dataset + "/"+strategy +"/Partitions/" + this->dataset + "_Partition_"+ partMethod+"_"+to_string(PartitionSize) + "_" +to_string(cluster_i);
            ofstream outFile(file_w_p, ios::out);
            if (!outFile) {
                cout << "Write File opening failed." << endl;
                assert(outFile);
                exit(1);
            }
            outFile << cluster_i<<" "<<partitionNodes.size()<<endl;//partition_id,partition_size,vertex_begin,vertex_end
            NodeId item_id;
            for (uint j = 0; j < partitionNodes.size(); ++j) {
                item_id = partitionNodes[j];
                outFile << item_id << " " << MCEdges[item_id].size();//partitions_info[cluster_i]+j+1
                for(auto it=MCEdges[item_id].begin();it!=MCEdges[item_id].end();++it){
                    outFile << " " << it->u;
                    for (auto it2 = it->w.begin(); it2 != it->w.end(); it2++) {
                        outFile << " " << *it2;
                    }
                }
                outFile << endl;
            }
            outFile.close();
        }
//        exit(0);
    }
    //function of converting edge weight to category
    int Graph_pre::Weight_to_category(int w){
        return floor(log2(w));
    }
    void Graph_pre::EMDIJk_preprocess(bool ifRead) {
        string file_r_graph = string(DataPath) + dataset+"/"+dataset+"_MCEdges_io_"+ to_string(PartitionSize) + ".txt";
        file_r_graph = string(DataPath) + this->dataset + "/" + this->dataset + "_MCEdgesS_"+ to_string(PartitionSize)+".txt";
        file_r_graph = string(DataPath) + this->dataset + "/" + this->dataset + "_MCEdgesS.txt";
        string file_r_graph1 = string(DataPath) + this->dataset + "/" + this->dataset + "_EM_Dijk_Index_"+ to_string(PartitionSize) + ".txt";

        vector<vector<bitset1<WeightPowMax+1,EdgeWeight>>> node_to_category;

        Timer tt;
        tt.start();

        //Construct the edge category index for all vertices
        int ID1,ID2;
        int degree,weight;
        int temp_int;
        int num_cri;
        string temp_cri;

        if(ifRead){//if we need to read graph data
            //Open file
            ifstream inFile(file_r_graph, ios::in);
            if (!inFile) {
                cout << "File opening failed." << endl;
                exit(1);
            }
            cout << "Graph Data loading for EM_Dijk preprocessing..." << endl;

            inFile >> node_num >> edge_num >> num_cri;
            cout<<"Node number: "<<node_num<<"\tEdge number: "<<edge_num<<endl;
            vector<int> temp_w;
            MCEdgesMap.assign(node_num,map<NodeId,MCEdge>());
            num_criteria = num_cri;
            mc_criteria.assign(num_cri,"");
            max_weights.assign(num_criteria,INT32_MIN);
            min_weights.assign(num_criteria,INT32_MAX);
            node_to_category.assign(num_criteria,vector<bitset1<WeightPowMax+1,EdgeWeight>>(node_num));//bitset<WeightPowMax+1>
            for(int i=0;i<num_cri;++i){
                inFile >> mc_criteria[i];
            }
            while(inFile) {
                inFile >> ID1 >> ID2;
                temp_w.clear();

                for(int i=0;i<num_cri;++i) {//get edge weights
                    inFile >> weight;
                    if(weight!=INF100){
                        max_weights[i] = max(max_weights[i],weight);
                        min_weights[i] = min(min_weights[i],weight);
                        node_to_category[i][ID1-1].set(Weight_to_category(weight)+1, true);//record edge category Weight_to_category(weight)+1
                    }
                    temp_w.push_back(weight);
                }
                MCEdgesMap[ID1-1][ID2-1]=MCEdge(ID2-1,temp_w);
//                MCEdgesMap[ID2-1][ID1-1]=MCEdge(ID1-1,temp_w);
            }
//
//            inFile >> node_num >> edge_num >> num_cri >> max_degree >> ave_degree;
//            for(int i=0;i<num_cri;++i){
//                inFile >> temp_cri;
//            }
//            max_weights.assign(num_cri,0);
//            assert(num_cri >= num_of_cri);
//            //get the maximal edge weights
//            for(int i=0;i<num_cri;++i){
//                inFile >> temp_int;
//                if(i<num_of_cri){
//                    max_weights[i]=temp_int;
//                }
//            }
//            // get node-to-category mapping
////            node_to_category.assign(NUM_OF_CRITERIA,vector<bitset1<WeightPowMax+1,EdgeWeight>>(node_num));//bitset<WeightPowMax+1>
//            for(int i=0;i<node_num;++i){
//                inFile>>ID1>>degree;
//                for(int j=0;j<degree;++j){
//                    inFile>>ID2;
//                    for(int k=0;k<num_cri;++k){
//                        inFile>>weight;
//                        if(k<num_of_cri){//only record the necessary information
//                            node_to_category[k][ID1-1].set(Weight_to_category(weight)+1, true);//record edge category Weight_to_category(weight)+1
//                        }
//                    }
//                }
//            }
            inFile.close();
            cout<<"Data load."<<endl;
        }else{//if no need to read
            node_to_category.assign(num_criteria,vector<bitset1<WeightPowMax+1,EdgeWeight>>(node_num));//bitset<WeightPowMax+1>
            rankSorter.rewind();
            ID1 = 0;
            while(!rankSorter.empty()){
                for(auto it1=MCEdgesMap[rankSorter->first].begin();it1!=MCEdgesMap[rankSorter->first].end();it1++){
                    for (int k=0; k< it1->second.w.size(); ++k) {
                        if(it1->second.w[k] != INF100){
                            node_to_category[k][ID1].set(Weight_to_category(it1->second.w[k])+1, true);//record edge category Weight_to_category(weight)+1
                        }
                    }
                }
                ++ID1; ++rankSorter;
            }
            rankSorter.clear();
        }
        /// generate category index file
        ofstream outFile(file_r_graph1, ios::out);
        if (!outFile) {
            cout << "Write File opening failed." << endl;
            assert(outFile);
            exit(1);
        }
        cout<<"Writing index of node_to_category to disk...\t";
        outFile << num_criteria << endl;
        for(int i=0;i<max_weights.size();++i){
            if(i<max_weights.size()-1){
                outFile << max_weights[i] << " ";
            }else{
                outFile << max_weights[i]<<endl;
            }
        }
        for(int i=0;i<node_num;++i){
            for(int j=0;j<node_to_category.size();++j){
                if(j<node_to_category.size()-1){
                    outFile<<node_to_category[j][i].m_val<<" ";
                }else{
                    outFile<<node_to_category[j][i].m_val<<endl;
                }
            }
        }
        outFile.close();
        cout<<"Done."<<endl;

        tt.stop();
        cout << "Done. The time used for data preprocessing is: " << tt.GetRuntime() <<" s."<<endl;
    }
    //Function for OD pairs generation
    void Graph_pre::ODpairGenerate(int times)
    {
        string RandomDis = string(DataPath) + dataset + "/" + strategy + "/" + dataset + "_OD_Random_"+ to_string(PartitionSize) +".txt";

        /*---OD pairs generation---*/
        int pairs = 0;
        int node_start, node_end;
        int temp = 0;
        ifstream inFile(RandomDis, ios::in);
        if(inFile.is_open()){
            inFile.close();
        }else{//if cannot open
            ofstream outFile(RandomDis, ios::out);
            if (!outFile) {
                cout << "File opening failed." << endl;
                assert(outFile);
            }
            cout << "OD pairs file generating..." << endl;
            outFile << times << endl;
            //generate random OD pairs
            pairs = 0;
            while (pairs < times) {
                node_start = rand() % node_num;
                node_end = rand() % node_num;
                while(node_end == node_start){
                    node_end = rand() % node_num;
                }
                outFile << node_start<< ' ' << node_end << endl;
                ++pairs;
            }

            outFile.close();
        }

        cout << "Finished." << endl;
    }
    //Function for OD pairs generation
    /*void Graph_pre::ODpairGenerate_Dis(int times, bool ifRead)//version 2: edge weight == 1
    {
        string file_r_mc = string(DataPath) + this->dataset + "/" + this->dataset + "_MCEdgesS_"+ to_string(PartitionSize)+".txt";

        string LongDis,MediumDis,ShortDis;
        LongDis = string(DataPath) + dataset + "/" + dataset + "_OD_LongDis_"+ to_string(PartitionSize) +".txt";
        MediumDis = string(DataPath) + dataset + "/" + dataset + "_OD_MediumDis_"+ to_string(PartitionSize) +".txt";
        ShortDis = string(DataPath) + dataset + "/" + dataset + "_OD_ShortDis_"+ to_string(PartitionSize) +".txt";
        benchmark::heap<2, NodeId, Distance> pqueue(node_num);
        NodeId item_id, temp_id;
        Distance item_dis, temp_dis;
        vector<bool> closed(node_num, false); //flag vector of whether closed
        vector<Distance> cost(node_num, INF);   //vector of cost
        int pairs = 0;
        int node_start, node_end;
        int temp = 0;
        vector<int> shortlist_S,shortlist_M,shortlist_L;
        uint low_bound = 2;//0.2;
        uint upper_bound = 8;//0.5;
        uint visited_edges = 0;

        ofstream outFile1(ShortDis, ios::out);
        if (!outFile1) {
            cout << "File opening failed." << endl;
            assert(outFile1);
        }
        ofstream outFile2(MediumDis, ios::out);
        if (!outFile2) {
            cout << "File opening failed." << endl;
            assert(outFile2);
        }
        ofstream outFile3(LongDis, ios::out);
        if (!outFile3) {
            cout << "File opening failed." << endl;
            assert(outFile3);
        }

        cout << "Coverage-based OD pairs files generating..." << endl;
        Timer tt;
        tt.start();

        outFile1 << times << endl;
        outFile2 << times << endl;
        outFile3 << times << endl;
        cout<<"Lower threshold: "<<low_bound<<"\t Upper threshold: "<<upper_bound<<endl;
        if(ifRead){
            Read_MCEdgesMapS(file_r_mc);

            cout<<"Lower bound: "<<low_bound<<"\tUpper bound: "<<upper_bound<<endl;
            pqueue.resize(node_num);
            //generate random OD pairs
            pairs = 0;
            while (pairs < times) {
                /// generate source node randomly
                node_start = rand() % node_num;
                shortlist_S.clear(); shortlist_M.clear(); shortlist_L.clear();
                closed.assign(node_num, false);
                cost.assign(node_num,INF);
                pqueue.clear();
                /// get the target node according to the ratio between visited number and graph node number
                cost[node_start] = 0;//cost of start node
                pqueue.update(node_start, 0);
                visited_number = 0;
                visited_edges = 0;
                //Iteration
                while (!pqueue.empty()) {//for every node in pqueue
                    pqueue.extract_min(item_id, item_dis);// top and delete min item
                    //relaxation
                    ++visited_number;//update #visited_number
                    for (auto it = MCEdgesMap[item_id].begin(); it != MCEdgesMap[item_id].end(); ++it) {//for MCEdgesMap
                        ++visited_edges;
                        temp_id = it->second.u;
                        if (closed[temp_id])//if closed
                            continue;
                        temp_dis = item_dis + 1;
                        if (cost[temp_id] > temp_dis) {//slack operation
                            cost[temp_id] = temp_dis;
                            pqueue.update(temp_id, temp_dis);
                        }
                    }

                    closed[item_id] = true;

                    if(item_dis > 0 && item_dis < low_bound){
                        shortlist_S.emplace_back(item_id);
                    }else if(item_dis >= low_bound && item_dis < upper_bound){
                        shortlist_M.emplace_back(item_id);
                    }else if(item_dis >= upper_bound){
                        shortlist_L.emplace_back(item_id);
                    }
                }
                if(pairs<10){
                    cout<<"Maximal hops: "<<item_dis<<endl;
                    cout<<"Size of lists: "<<shortlist_S.size()<<", "<<shortlist_M.size()<<", "<<shortlist_L.size()<<endl;
                    if(pairs == 0){
                        low_bound = item_dis/5;//node_num*0.1;
                        upper_bound = item_dis/2;//node_num*0.5;
                        cout<<"Change lower threshold to "<<low_bound<<" , upper threshold to "<<upper_bound<<endl;
                    }
                }

                //for short distance
                node_end = shortlist_S[rand() % shortlist_S.size()];//
                outFile1 << node_start+1 << ' ' << node_end+1 << endl;
                //for medium distance
                node_end = shortlist_M[rand() % shortlist_M.size()];
                outFile2 << node_start+1 << ' ' << node_end+1 << endl;
                //for long distance
                node_end = shortlist_L[rand() % shortlist_L.size()];
                outFile3 << node_start+1 << ' ' << node_end+1 << endl;
                ++pairs;
            }

        }else{
            //generate random OD pairs
            pairs = 0;
            cout<<"Lower bound: "<<low_bound<<"\tUpper bound: "<<upper_bound<<endl;
            while (pairs < times) {
                /// generate source node randomly
                node_start = rand() % node_num;
                shortlist_S.clear(); shortlist_M.clear(); shortlist_L.clear();
                closed.assign(node_num, false);
                cost.assign(node_num,INF);
                pqueue.clear();
                /// get the target node according to the ratio between visited number and graph node number
                cost[node_start] = 0;//cost of start node
                pqueue.update(node_start, 0);
                visited_number = 0;
                visited_edges = 0;
                //Iteration
                while (!pqueue.empty()) {//for every node in pqueue
                    pqueue.extract_min(item_id, item_dis);// top and delete min item
                    //relaxation
                    ++visited_number;//update #visited_number

                    for (auto it = MCEdgesMap[item_id].begin(); it != MCEdgesMap[item_id].end(); ++it) {//for MCEdgesMap
                        ++visited_edges;
                        temp_id = it->second.u;
                        if (closed[temp_id])//if closed
                            continue;
                        temp_dis = item_dis + 1;
                        if (cost[temp_id] > temp_dis) {//slack operation
                            cost[temp_id] = temp_dis;
                            pqueue.update(temp_id, temp_dis);
                        }
                    }
                    closed[item_id] = true;

                    if(item_dis > 0 && item_dis < low_bound){
                        shortlist_S.emplace_back(item_id);
                    }else if(item_dis >= low_bound && item_dis < upper_bound){
                        shortlist_M.emplace_back(item_id);
                    }else if(item_dis >= upper_bound){
                        shortlist_L.emplace_back(item_id);
                    }
                }
                if(pairs<10){
                    cout<<"Maximal hops: "<<item_dis<<endl;
                    cout<<"Size of lists: "<<shortlist_S.size()<<", "<<shortlist_M.size()<<", "<<shortlist_L.size()<<endl;
                    if(pairs == 0){
                        low_bound = item_dis/5;//node_num*0.1;
                        upper_bound = item_dis/2;//node_num*0.5;
                        cout<<"Change lower threshold to "<<low_bound<<" , upper threshold to "<<upper_bound<<endl;
                    }
                }
                //for short distance
                node_end = shortlist_S[rand() % shortlist_S.size()];
                outFile1 << node_start+1 << ' ' << node_end+1 << endl;
                //for medium distance
                node_end = shortlist_M[rand() % shortlist_M.size()];
                outFile2 << node_start+1 << ' ' << node_end+1 << endl;
                //for long distance
                node_end = shortlist_L[rand() % shortlist_L.size()];
                outFile3 << node_start+1 << ' ' << node_end+1 << endl;
                ++pairs;
            }
        }
        outFile1.close();
        outFile2.close();
        outFile3.close();
        tt.stop();
        cout << "Finished. The time used for queries generation is " << tt.GetRuntime() << " s." << endl;
    }*/
    void Graph_pre::ODpairGenerate_Dis(int times, bool ifRead)//version 1: based on Dijkstra's algorithm
    {
        string file_r_mc = string(DataPath) + this->dataset + "/" + this->dataset + "_MCEdgesS_"+ to_string(PartitionSize)+".txt";

        string LongDis,MediumDis,ShortDis;
        LongDis = string(DataPath) + dataset + "/" + strategy + "/" + dataset + "_OD_LongDis_"+ to_string(PartitionSize) +".txt";
        MediumDis = string(DataPath) + dataset + "/" + strategy + "/" +  dataset + "_OD_MediumDis_"+ to_string(PartitionSize) +".txt";
        ShortDis = string(DataPath) + dataset + "/" + strategy + "/" + dataset + "_OD_ShortDis_"+ to_string(PartitionSize) +".txt";
        benchmark::heap<2, NodeId, Distance> pqueue(node_num);
        NodeId item_id, temp_id;
        Distance item_dis, temp_dis;
        vector<bool> closed(node_num, false); //flag vector of whether closed
        vector<Distance> cost(node_num, INF);   //vector of cost
        int pairs = 0;
        int node_start, node_end;
        int temp = 0;
        vector<int> shortlist_S,shortlist_M,shortlist_L;
        uint low_bound = edge_num/10;//node_num*0.1;
        uint upper_bound = edge_num/2;//node_num*0.5;
        unsigned long long visited_edges = 0;

        ofstream outFile1(ShortDis, ios::out);
        if (!outFile1) {
            cout << "File opening failed." << endl;
            assert(outFile1);
        }
        ofstream outFile2(MediumDis, ios::out);
        if (!outFile2) {
            cout << "File opening failed." << endl;
            assert(outFile2);
        }
        ofstream outFile3(LongDis, ios::out);
        if (!outFile3) {
            cout << "File opening failed." << endl;
            assert(outFile3);
        }

        cout << "Coverage-based OD pairs files generating..." << endl;
        Timer tt;
        tt.start();

        outFile1 << times << endl;
        outFile2 << times << endl;
        outFile3 << times << endl;

        if(ifRead){
            Read_MCEdgesMapS(file_r_mc);
            low_bound = edge_num/10;//node_num*0.1;
            upper_bound = edge_num/2;//node_num*0.5;
            cout<<"Lower bound: "<<low_bound<<"\tUpper bound: "<<upper_bound<<endl;
            pqueue.resize(node_num);
            //generate random OD pairs
            pairs = 0;
            while (pairs < times) {
                /// generate source node randomly
                node_start = rand() % node_num;
                shortlist_S.clear(); shortlist_M.clear(); shortlist_L.clear();
                closed.assign(node_num, false);
                cost.assign(node_num,INF);
                pqueue.clear();
                /// get the target node according to the ratio between visited number and graph node number
                cost[node_start] = 0;//cost of start node
                pqueue.update(node_start, 0);
                visited_number = 0;
                visited_edges = 0;
                //Iteration
                while (!pqueue.empty()) {//for every node in pqueue
                    pqueue.extract_min(item_id, item_dis);// top and delete min item
                    //relaxation
                    ++visited_number;//update #visited_number
                    for (auto it = MCEdgesMap[item_id].begin(); it != MCEdgesMap[item_id].end(); ++it) {//for MCEdgesMap
                        ++visited_edges;
                        temp_id = it->second.u;
                        if (closed[temp_id])//if closed
                            continue;
                        temp_dis = item_dis + it->second.w[0]; // use first criterion
                        if (cost[temp_id] > temp_dis) {//slack operation
                            cost[temp_id] = temp_dis;
                            pqueue.update(temp_id, temp_dis);
                        }
                    }

                    closed[item_id] = true;

                    if(visited_edges >= 1000 && visited_edges < low_bound){
                        shortlist_S.emplace_back(item_id);
                    }else if(visited_edges >= low_bound && visited_edges < upper_bound){
                        shortlist_M.emplace_back(item_id);
                    }else if(visited_edges >= upper_bound){
                        shortlist_L.emplace_back(item_id);
                    }
                }

                //for short distance
                node_end = shortlist_S[rand() % shortlist_S.size()];
                outFile1 << node_start << ' ' << node_end << endl;
                //for medium distance
                node_end = shortlist_M[rand() % shortlist_M.size()];
                outFile2 << node_start << ' ' << node_end << endl;
                //for long distance
                node_end = shortlist_L[rand() % shortlist_L.size()];
                outFile3 << node_start << ' ' << node_end << endl;
                ++pairs;
            }

        }else{
            //generate random OD pairs
            pairs = 0;
            cout<<"Lower bound: "<<low_bound<<"\tUpper bound: "<<upper_bound<<endl;
            while (pairs < times) {
                /// generate source node randomly
                node_start = rand() % node_num;
                shortlist_S.clear(); shortlist_M.clear(); shortlist_L.clear();
                closed.assign(node_num, false);
                cost.assign(node_num,INF);
                pqueue.clear();
                /// get the target node according to the ratio between visited number and graph node number
                cost[node_start] = 0;//cost of start node
                pqueue.update(node_start, 0);
                visited_number = 0;
                visited_edges = 0;
                //Iteration
                while (!pqueue.empty()) {//for every node in pqueue
                    pqueue.extract_min(item_id, item_dis);// top and delete min item
                    //relaxation
                    ++visited_number;//update #visited_number

                    for (auto it = MCEdgesMap[item_id].begin(); it != MCEdgesMap[item_id].end(); ++it) {//for MCEdgesMap
                        ++visited_edges;
                        temp_id = it->second.u;
                        if (closed[temp_id])//if closed
                            continue;
                        temp_dis = item_dis + it->second.w[0];
                        if (cost[temp_id] > temp_dis) {//slack operation
                            cost[temp_id] = temp_dis;
                            pqueue.update(temp_id, temp_dis);
                        }
                    }
                    closed[item_id] = true;

                    if(visited_edges >= 1000 && visited_edges < low_bound){
                        shortlist_S.emplace_back(item_id);
                    }else if(visited_edges >= low_bound && visited_edges < upper_bound){
                        shortlist_M.emplace_back(item_id);
                    }else if(visited_edges >= upper_bound){
                        shortlist_L.emplace_back(item_id);
                    }
                }
                //for short distance
                node_end = shortlist_S[rand() % shortlist_S.size()];
                outFile1 << node_start << ' ' << node_end << endl;
                //for medium distance
                node_end = shortlist_M[rand() % shortlist_M.size()];
                outFile2 << node_start << ' ' << node_end << endl;
                //for long distance
                node_end = shortlist_L[rand() % shortlist_L.size()];
                outFile3 << node_start << ' ' << node_end << endl;
                ++pairs;
            }
        }
        outFile1.close();
        outFile2.close();
        outFile3.close();
        tt.stop();
        cout << "Finished. The time used for queries generation is " << tt.GetRuntime() << " s." << endl;
    }
    /*void Graph_pre::ODpairGenerate_Dis2(int times, bool ifRead){//version 1: based on Dijkstra's algorithm
        string file_r_mc = string(DataPath) + dataset + "/" + strategy + "/" +dataset + "_"+ to_string(PartitionSize)+".MCEdges";

        string LongDis,MediumDis,ShortDis;
        LongDis = string(DataPath) + dataset + "/" + strategy + "/" + dataset + "_OD_LongDis_"+ to_string(PartitionSize) +".txt";
        MediumDis = string(DataPath) + dataset + "/" + strategy + "/" +  dataset + "_OD_MediumDis_"+ to_string(PartitionSize) +".txt";
        ShortDis = string(DataPath) + dataset + "/" + strategy + "/" + dataset + "_OD_ShortDis_"+ to_string(PartitionSize) +".txt";
        benchmark::heap<2, NodeId, Distance> pqueue(node_num);
        NodeId item_id, temp_id;
        Distance item_dis, temp_dis;
        vector<bool> closed(node_num, false); //flag vector of whether closed
        vector<Distance> cost(node_num, INF);   //vector of cost
        int pairs = 0;
        int node_start, node_end;
        int temp = 0;
        vector<int> shortlist_S,shortlist_M,shortlist_L;
        uint low_bound = edge_num*0.1;//node_num*0.1;
        uint upper_bound = edge_num*0.4;//node_num*0.5;
        unsigned long long visited_edges = 0;

        ofstream outFile1(ShortDis, ios::out);
        if (!outFile1) {
            cout << "File opening failed." << endl;
            assert(outFile1);
        }
        ofstream outFile2(MediumDis, ios::out);
        if (!outFile2) {
            cout << "File opening failed." << endl;
            assert(outFile2);
        }
        ofstream outFile3(LongDis, ios::out);
        if (!outFile3) {
            cout << "File opening failed." << endl;
            assert(outFile3);
        }

        cout << "Coverage-based OD pairs files generating..." << endl;
        Timer tt;
        tt.start();

        outFile1 << times << endl;
        outFile2 << times << endl;
        outFile3 << times << endl;

        if(ifRead){
            ReadMCEdges(file_r_mc);
            cout<<"Lower bound: "<<low_bound<<"\tUpper bound: "<<upper_bound<<endl;
            //generate random OD pairs
            pairs = 0;
            while (pairs < times) {
                /// generate source node randomly
                node_start = rand() % node_num;
                shortlist_S.clear(); shortlist_M.clear(); shortlist_L.clear();
                closed.assign(node_num, false);
                cost.assign(node_num,INF);
                pqueue.clear();
                /// get the target node according to the ratio between visited number and graph node number
                cost[node_start] = 0;//cost of start node
                pqueue.update(node_start, 0);
                visited_number = 0;
                visited_edges = 0;
                //Iteration
                while (!pqueue.empty()) {//for every node in pqueue
                    pqueue.extract_min(item_id, item_dis);// top and delete min item
                    //relaxation
                    ++visited_number;//update #visited_number
                    for (auto it = MCEdges[item_id].begin(); it != MCEdges[item_id].end(); ++it) {//for MCEdgesMap
                        ++visited_edges;
                        temp_id = it->u;
                        if (closed[temp_id])//if closed
                            continue;
                        temp_dis = item_dis + it->w[0]; // use first criterion
                        if (cost[temp_id] > temp_dis) {//slack operation
                            cost[temp_id] = temp_dis;
                            pqueue.update(temp_id, temp_dis);
                        }
                    }

                    closed[item_id] = true;

                    if(visited_edges >= 10 && visited_edges < low_bound){
                        shortlist_S.emplace_back(item_id);
                    }else if(visited_edges >= low_bound && visited_edges < upper_bound){
                        shortlist_M.emplace_back(item_id);
                    }else if(visited_edges >= upper_bound){
                        shortlist_L.emplace_back(item_id);
                    }
                }

                //for short distance
                node_end = shortlist_S[rand() % shortlist_S.size()];
                outFile1 << node_start << ' ' << node_end << endl;
                //for medium distance
                node_end = shortlist_M[rand() % shortlist_M.size()];
                outFile2 << node_start << ' ' << node_end << endl;
                //for long distance
                node_end = shortlist_L[rand() % shortlist_L.size()];
                outFile3 << node_start << ' ' << node_end << endl;
                ++pairs;
            }

        }
        else{
            //generate random OD pairs
            pairs = 0;
            cout<<"Lower bound: "<<low_bound<<"\tUpper bound: "<<upper_bound<<endl;
            while (pairs < times) {
                /// generate source node randomly
                node_start = rand() % node_num;
                shortlist_S.clear(); shortlist_M.clear(); shortlist_L.clear();
                closed.assign(node_num, false);
                cost.assign(node_num,INF);
                /// get the target node according to the ratio between visited number and graph node number
                cost[node_start] = 0;//cost of start node
                pqueue.update(node_start, 0);
                visited_number = 0;
                visited_edges = 0;
                //Iteration
                while (!pqueue.empty()) {//for every node in pqueue
                    pqueue.extract_min(item_id, item_dis);// top and delete min item
                    //relaxation
                    ++visited_number;//update #visited_number

                    for (auto it = MCEdges[item_id].begin(); it != MCEdges[item_id].end(); ++it) {//for MCEdgesMap
                        ++visited_edges;
                        temp_id = it->u;
                        if (closed[temp_id])//if closed
                            continue;
                        temp_dis = item_dis + it->w[0];
                        if (cost[temp_id] > temp_dis) {//slack operation
                            cost[temp_id] = temp_dis;
                            pqueue.update(temp_id, temp_dis);
                        }
                    }
                    closed[item_id] = true;

                    if(visited_edges >= 10 && visited_edges < low_bound){
                        shortlist_S.emplace_back(item_id);
                    }else if(visited_edges >= low_bound && visited_edges < upper_bound){
                        shortlist_M.emplace_back(item_id);
                    }else if(visited_edges >= upper_bound){
                        shortlist_L.emplace_back(item_id);
                    }
                }
                //for short distance
                node_end = shortlist_S[rand() % shortlist_S.size()];
                outFile1 << nodes_map[node_start] << ' ' << nodes_map[node_end] << endl;
                //for medium distance
                node_end = shortlist_M[rand() % shortlist_M.size()];
                outFile2 << nodes_map[node_start] << ' ' << nodes_map[node_end] << endl;
                //for long distance
                node_end = shortlist_L[rand() % shortlist_L.size()];
                outFile3 << nodes_map[node_start] << ' ' << nodes_map[node_end] << endl;
                ++pairs;
            }
        }
        outFile1.close();
        outFile2.close();
        outFile3.close();
        tt.stop();
        cout << "Finished. The time used for queries generation is " << tt.GetRuntime() << " s." << endl;
    }*/
    //stxxl::priority queue
    /*void Graph_pre::CoverageDijkstra(int node_start, vector<vector<pair<int,int>>>& queries, unsigned long long lower_bound, unsigned long long upper_bound) {
        ioxxl::PriorityQueue pqueue(ioxxl::PQ_Pool);
//        priority_queue<VertexCost,vector<VertexCost>,PQCompareLess> pqueue;

        NodeId item_id, temp_id;
        Distance item_dis, temp_dis;
        vector<bool> closed(node_num, false); //flag vector of whether closed
        vector<Distance> cost(node_num, INF);   //vector of cost
        VertexCost item_;
        unsigned long long visited_edges = 0;
        vector<int> shortlist_S,shortlist_M,shortlist_L;
        pqueue.push(VertexCost(node_start, 0));
        cost[node_start]=0;


        //Iteration
        while (!pqueue.empty()) {//for every node in pqueue

            item_ = pqueue.top();// top min item
            item_id = item_.id;
            pqueue.pop();
            if (closed[item_id])
                continue;

            item_dis = cost[item_id];
            closed[item_id] = true;
            for (auto it = MCEdges[item_id].begin(); it != MCEdges[item_id].end(); ++it) {
                temp_id = it->u;
                ++visited_edges;
                if (closed[temp_id])
                    continue;
                temp_dis = item_dis + it->w[0];
                if (cost[temp_id] > temp_dis) {//slack operation
                    cost[temp_id] = temp_dis;
                    pqueue.push(VertexCost(temp_id, temp_dis));
                }
            }

            if(visited_edges >= 10 && visited_edges < lower_bound){
                if(item_id!=node_start){
                    shortlist_S.emplace_back(item_id);
                }
            }else if(visited_edges >= lower_bound && visited_edges < upper_bound){
                shortlist_M.emplace_back(item_id);
            }else if(visited_edges >= upper_bound){
                while(!pqueue.empty()){
                    pqueue.pop();
                }
                break;
            }
        }

        for(int id=0;id<node_num;++id){
            if(!closed[id]){
                shortlist_L.emplace_back(id);
            }
        }

//        cout<<"Bound: "<<lower_bound<<" "<<upper_bound<<" "<<visited_edges<<endl;
//        cout<<"Vector: "<<shortlist_S.size()<<" "<<shortlist_M.size()<<" "<<shortlist_L.size()<<endl;
        //for short distance
        vector<pair<int,int>> query;
        int node_end=0;
        node_end = shortlist_S[rand() % shortlist_S.size()];
        query.emplace_back(node_start,node_end);
        //for medium distance
        node_end = shortlist_M[rand() % shortlist_M.size()];
        query.emplace_back(node_start,node_end);
        //for long distance
        node_end = shortlist_L[rand() % shortlist_L.size()];
        query.emplace_back(node_start,node_end);
        sm->wait();
        queries.push_back(query);
        sm->notify();
    }*/

    //benchmakr::priority queue
    void Graph_pre::CoverageDijkstra(int node_start, vector<vector<pair<int,int>>>& queries, unsigned long long lower_bound, unsigned long long upper_bound) {
//        ioxxl::PriorityQueue pqueue(ioxxl::PQ_Pool);
//        priority_queue<VertexCost,vector<VertexCost>,PQCompareLess> pqueue;
        benchmark::heap<2, NodeId, Distance> pqueue(node_num);

        NodeId item_id, temp_id;
        Distance item_dis, temp_dis;
        vector<bool> closed(node_num, false); //flag vector of whether closed
        vector<Distance> cost(node_num, INF);   //vector of cost
        VertexCost item_;
        unsigned long long visited_edges = 0;
        vector<int> shortlist_S,shortlist_M,shortlist_L;
        pqueue.update(node_start, 0);
        cost[node_start]=0;


        //Iteration
        while (!pqueue.empty()) {//for every node in pqueue
            pqueue.extract_min(item_id, item_dis);// top and delete min item

            item_dis = cost[item_id];
            closed[item_id] = true;
            for (auto it = MCEdges[item_id].begin(); it != MCEdges[item_id].end(); ++it) {
                temp_id = it->u;
                ++visited_edges;
                if (closed[temp_id])
                    continue;
                temp_dis = item_dis + it->w[0];
                if (cost[temp_id] > temp_dis) {//slack operation
                    cost[temp_id] = temp_dis;
                    pqueue.update(temp_id, temp_dis);
                }
            }

            if(visited_edges >= 10 && visited_edges < lower_bound){
                if(item_id!=node_start){
                    shortlist_S.emplace_back(item_id);
                }
            }else if(visited_edges >= lower_bound && visited_edges < upper_bound){
                shortlist_M.emplace_back(item_id);
            }else if(visited_edges >= upper_bound){
                break;
            }
        }

        for(int id=0;id<node_num;++id){
            if(!closed[id]){
                shortlist_L.emplace_back(id);
            }
        }

//        cout<<"Bound: "<<lower_bound<<" "<<upper_bound<<" "<<visited_edges<<endl;
//        cout<<"Vector: "<<shortlist_S.size()<<" "<<shortlist_M.size()<<" "<<shortlist_L.size()<<endl;
        //for short distance
        vector<pair<int,int>> query;
        int node_end=0;
        node_end = shortlist_S[rand() % shortlist_S.size()];
        query.emplace_back(node_start,node_end);
        //for medium distance
        node_end = shortlist_M[rand() % shortlist_M.size()];
        query.emplace_back(node_start,node_end);
        //for long distance
        node_end = shortlist_L[rand() % shortlist_L.size()];
        query.emplace_back(node_start,node_end);
        sm->wait();
        queries.push_back(query);
        sm->notify();
    }

    void Graph_pre::CoverageDijkstraV(vector<int>& p, vector<vector<pair<int,int>>>& queries, unsigned long long lower_bound, unsigned long long upper_bound){
        int node_start;
        for(int i=0;i<p.size();++i){
            node_start = p[i];
//            cout<<"Vertex "<<i<<" "<< node_start <<endl;
            CoverageDijkstra(node_start, queries, lower_bound, upper_bound);
        }
    }

    //function of generating the coverage-based query sets
    void Graph_pre::ODpairGenerate_Dis2(int times, bool ifRead){//version 1: based on Dijkstra's algorithm
        string file_r_mc = string(DataPath) + dataset + "/" + strategy + "/" +dataset + "_"+ to_string(PartitionSize)+".MCEdges";
        file_r_mc = string(DataPath) + dataset + "/" + dataset + ".MCEdges";

        string LongDis,MediumDis,ShortDis;
        LongDis = string(DataPath) + dataset + "/" + strategy + "/" + dataset + "_OD_LongDis_"+ to_string(PartitionSize) +".txt";
        MediumDis = string(DataPath) + dataset + "/" + strategy + "/" +  dataset + "_OD_MediumDis_"+ to_string(PartitionSize) +".txt";
        ShortDis = string(DataPath) + dataset + "/" + strategy + "/" + dataset + "_OD_ShortDis_"+ to_string(PartitionSize) +".txt";
        int pairs = 0;
        int node_start, node_end;

        vector<vector<pair<int,int>>> queries;

        Timer tt;
        tt.start();

        ifstream inFile(ShortDis, ios::in);
        if(inFile.is_open() && !ifRead){
            inFile.close();
        }else{
            inFile.close();
            ofstream outFile1(ShortDis, ios::out);
            if (!outFile1) {
                cout << "File opening failed. " <<ShortDis<< endl;
                exit(1);
            }
            ofstream outFile2(MediumDis, ios::out);
            if (!outFile2) {
                cout << "File opening failed. " <<MediumDis<< endl;
                exit(1);
            }
            ofstream outFile3(LongDis, ios::out);
            if (!outFile3) {
                cout << "File opening failed. " <<LongDis<< endl;
                exit(1);
            }

            sm = new Semaphore(1);
            cout << "Coverage-based OD pairs files generating..." << endl;

            num_criteria = 5;

            if(ifRead){
                ReadMCEdges(file_r_mc);
            }
            unsigned long long lower_bound = edge_num*0.1;//node_num*0.1;
            unsigned long long upper_bound = edge_num*0.4;//node_num*0.5;
            cout<<"Lower bound: "<<lower_bound<<"\tUpper bound: "<<upper_bound<<endl;
            vector<int> sources;
            vector<vector<int>> ProcessID(threadnum, vector<int>());
            //generate random OD pairs
            for(int i=0;i<times;++i){
                node_start = rand() % node_num;
                sources.emplace_back(node_start);
                ProcessID[i%threadnum].emplace_back(node_start);
            }

            boost::thread_group thread;
            for(int i=0;i<ProcessID.size();i++){
//            cout<<"ProcessID "<<i<<endl;
                vector<int> p=ProcessID[i];
                thread.add_thread(new boost::thread(&Graph_pre::CoverageDijkstraV, this, boost::ref(ProcessID[i]), boost::ref(queries), lower_bound, upper_bound));
//            CoverageDijkstraV(ProcessID[i], queries, lower_bound, upper_bound);
            }
            thread.join_all();


            if(queries.size() != times){
                cout<<"Inconsistent! "<<queries.size()<<endl;
                exit(1);
            }

            outFile1 << times << endl;
            outFile2 << times << endl;
            outFile3 << times << endl;

            for(int i=0;i<times;++i){
                //for short distance
                outFile1 << queries[i][0].first << ' ' << queries[i][0].second << endl;
                //for medium distance
                outFile2 << queries[i][1].first << ' ' << queries[i][1].second << endl;
                //for long distance
                outFile3 << queries[i][2].first << ' ' << queries[i][2].second << endl;
            }

            outFile1.close();
            outFile2.close();
            outFile3.close();
        }


        tt.stop();
        cout << "Finished. The time used for queries generation is " << tt.GetRuntime() << " s." << endl;
    }
    //Function for OD pairs remapping
    //Function for OD pairs generation
    void Graph_pre::ODpairGenerate_Remap()
    {
        string qtype[3] = {"ShortDis","MediumDis","LongDis"};
        for(int i=0;i<3;++i){
            string file_r = string(DataPath) + dataset + "/ave/" + dataset + "_OD_"+ qtype[i]+"_"+to_string(PartitionSize) +".txt";
            string file_w = string(DataPath) + dataset + "/" + dataset + "_OD_"+ qtype[i]+"_"+to_string(PartitionSize) +".txt";
            ifstream inFile(file_r, ios::in);
            if (!inFile) {
                cout << "File opening failed." << endl;
                assert(inFile);
            }
            ofstream outFile(file_w, ios::out);
            if (!outFile) {
                cout << "File opening failed." << endl;
                assert(outFile);
            }
            int times;
            int ID1,ID2;
            inFile >> times;
            outFile << times<<endl;
            while(inFile){
                inFile >> ID1 >> ID2;
                outFile <<ID1<<" "<<ID2<<endl;
            }
            inFile.close();
            outFile.close();
            cout<<qtype[i]<<" remap writing done."<<endl;
        }
        cout << "Finished." << endl;
    }
    //function for generating edge weights of basic criterion for graph, memory usage is 2*|E|
    void Graph_pre::BaseCriGenerate(uint n_num){
        //normal distribution generator
        Timer tt;
        tt.start();
        string lineStr;
        string line_symbol;
        string temp_str;
        int ID1, ID2, weight=1;
        int temp_id = -1;
        uint temp_degree = 0;

        uint temp_num = 0;
        node_num = n_num;
        edge_num = 0;
        max_degree = 0;
        min_degree = INF;


        unsigned long long edgeNum=0;
        unordered_set<int> set_A; set<int> set_LCC;
        set_A.clear();

        vector<int> tempCri(num_criteria,0);
        vector<unordered_map<int,vector<int>>> pairGraph;
        vector<int> degrees;

        ifstream inFilen(DataPath + dataset + "/" + dataset + ".MCEdges",ios::in);
        if(inFilen.is_open()){//if exists
            cout<<"The MCEdges file exists!"<<endl;
            num_criteria = 5;//update the number of criteria
        }
        else{
            /// generate basic criteria
            if(num_criteria==1){//if there is only one existing criterion
                string r_file = string(DataPath) + dataset + "/" + dataset + ".original";
                string w_file = string(DataPath) + dataset + "/" + dataset + ".Distance";

                ifstream inFile2(w_file, ios::in);
                if (!inFile2) {
                    cout << "Generating the data of first criterion." << endl;
                    ifstream inFile(r_file, ios::in);
                    if (!inFile) {
                        cout << "File opening failed." << endl;
                        assert(inFile);
                        exit(1);
                    }

                    /// read graph and recording the degree of vertices
                    cout<<"Nominate node number: "<<n_num<<endl;

                    string line;
                    getline(inFile,line);
                    vector<string> re1;
                    boost::split(re1,line,boost::is_any_of(" \t"));
                    if(re1.size()==2){
                        node_num=stoi(re1[0]), edge_num=stoul(re1[1]); n_num=node_num;
                    }else if(re1.size()>=3){
                        if( re1[0]=="#"){
                            node_num=stoi(re1[1]), edge_num=stoul(re1[2]); n_num=node_num;
                        }else if(re1[0]=="%"){
                            cout<<"Node number is needed! " <<line<<endl;
                        }else{
                            cout<<"Wrong input! re.size: "<<re1.size()<<". " <<line<<endl;
                            exit(1);
                        }

                    }else{
                        cout<<"Wrong input! re.size: "<<re1.size()<<". " <<line<<endl;
                        exit(1);
                    }

                    cout<<"Original node number: "<<node_num<<" ; edge number: "<<edge_num<<endl;
                    degrees.assign(node_num,0);
                    pairGraph.assign(node_num,unordered_map<int,vector<int>>());

                    while(getline(inFile,line)){
                        if(line.empty())
                            continue;
                        re1.clear();
                        boost::split(re1,line,boost::is_any_of(" \t"));
                        if(re1.size()>=2){
                            ID1=stoi(re1[0]), ID2=stoi(re1[1]);
                        }
                        else{
                            cout<<"Wrong line. "<< line<<endl;
                            exit(1);
                        }

                        ID1-=1; ID2-=1;

                        if(ID1>=0 && ID1<node_num && ID2>=0 && ID2<node_num && ID1!=ID2){
                            if(pairGraph[ID1].find(ID2)==pairGraph[ID1].end()){//if not found
                                pairGraph[ID1].insert({ID2,tempCri});
                                pairGraph[ID2].insert({ID1,tempCri});
                                edgeNum+=2;
                                set_A.insert(ID1); set_A.insert(ID2);
                                degrees[ID1]++; degrees[ID2]++;
                            }

                            if(weight<=0){
                                cout<<"Negative edge weight "<<ID1<<" "<<ID2<<" "<<weight<<endl;
                            }
                        }else{
                            if(ID1!=ID2){
                                cout<<"Graph data is wrong! "<<ID1<<" "<<ID2<<" "<<weight<<endl;
                                exit(1);
                            }
                        }
                    }

                    inFile.close();
                    tt.stop();
                    cout << "The time for data loading: " << tt.GetRuntime() << " s." << endl;
                    temp_num = set_A.size();
                    if(temp_num != node_num){
                        cout << "\n!!! The valid node number is " << temp_num <<" , the valid edge number is "<<edgeNum<<endl;
                    }

                    max_degree = *max_element(degrees.begin(),degrees.end());
                    min_degree = *min_element(degrees.begin(),degrees.end());
                    cout << "Maximal degree: " << max_degree << " \tMinimal degree: " << min_degree <<endl;

                    /// generate the basic criterion data
                    int wei;
                    for(ID1=0;ID1<node_num;++ID1){
                        for(auto it=pairGraph[ID1].begin();it!=pairGraph[ID1].end();++it){
                            ID2 = it->first;
                            wei = degrees[ID1]+degrees[ID2];
                            weight = 2*max_degree/wei;
                            if(weight>0){
                                pairGraph[ID1][ID2][0] = weight;
                                base_max = max(base_max,weight);
                                base_min = min(base_min,weight);
                            }else{
                                cout<<"Wrong edge weight. "<<weight<<endl;
                            }
                        }
                    }

                    pair<int,unsigned long long> nums = DFS_CC(pairGraph, set_A, set_LCC, node_num);
                    node_num=nums.first, edge_num=nums.second;
                    cout<<"Node number: "<<node_num<<" ; edge number: "<<edge_num<<endl;

                    WriteGraph(w_file,set_LCC, pairGraph, 0);
                }
                else{//if the basic criterion exists
                    cout<<"The basic criterion data exists!"<<endl;
                    string line;
                    getline(inFile2,line);
                    vector<string> re1;
                    boost::split(re1,line,boost::is_any_of(" \t"));
                    if(re1.size()==2){
                        node_num=stoi(re1[0]), edge_num=stoul(re1[1]);
                    }else{
                        cout<<"Wrong input! re.size: "<<re1.size()<<". " <<line<<endl;
                        exit(1);
                    }
                    cout<<"Node number: "<<node_num<<" ; edge number: "<<edge_num<<endl;
                    degrees.assign(node_num,0);
                    pairGraph.assign(node_num,unordered_map<int,vector<int>>());

                    while(getline(inFile2,line)){
                        if(line.empty())
                            continue;
                        re1.clear();
                        boost::split(re1,line,boost::is_any_of(" \t"));
                        if(re1.size()==3){
                            ID1=stoi(re1[0]), ID2=stoi(re1[1]), weight=stoi(re1[2]);
                        }
                        else{
                            cout<<"Wrong line. "<< line<<endl;
                            exit(1);
                        }


                        if(ID1>=0 && ID1<node_num && ID2>=0 && ID2<node_num && ID1!=ID2 && weight>0){
                            if(pairGraph[ID1].find(ID2)==pairGraph[ID1].end()){//if not found
                                tempCri[0]=weight;
                                pairGraph[ID1].insert({ID2,tempCri});
                                pairGraph[ID2].insert({ID1,tempCri});
                                base_max = max(base_max,weight);
                                base_min = min(base_min,weight);
                                edgeNum+=2;
                                set_LCC.insert(ID1); set_LCC.insert(ID2);
                                degrees[ID1]++; degrees[ID2]++;
                            }

                        }else{
                            if(ID1!=ID2){
                                cout<<"Graph data is wrong! "<<ID1<<" "<<ID2<<" "<<weight<<endl;
                                exit(1);
                            }
                        }
                    }
                    inFile2.close();

                    max_degree = *max_element(degrees.begin(),degrees.end());
                    min_degree = *min_element(degrees.begin(),degrees.end());
                    cout << "Maximal degree: " << max_degree << " \tMinimal degree: " << min_degree <<endl;
                }

            }
            else if(num_criteria==2){//if there is two existing criteria, i.e., road network
                string r_file = string(DataPath) + dataset + "/" + dataset + "_Distance.gr";
                string w_file = string(DataPath) + dataset + "/" + dataset + ".Distance";
                string r_file2 = string(DataPath) + dataset + "/" + dataset + "_Time.gr";
                string w_file2 = string(DataPath) + dataset + "/" + dataset + ".Time";


                ifstream inFilew(w_file, ios::in);
                ifstream inFilew2(w_file2, ios::in);
                if (!inFilew.is_open() || !inFilew2.is_open()) {
                    cout<<"Generating data..."<<endl;
                    ifstream inFile(r_file, ios::in);
                    if (!inFile) {
                        cout << "File opening failed." << r_file<<endl;
                        exit(1);
                    }
                    ifstream inFile2(r_file2, ios::in);
                    if (!inFile2) {
                        cout << "File opening failed." << r_file2<<endl;
                        exit(1);
                    }

                    /// read graph and recording the degree of vertices
                    cout<<"Nominate node number: "<<n_num<<endl;
                    string line,line2;
                    vector<string> re1, re2;
                    int weight2=1;

                    degrees.assign(node_num,0);
                    pairGraph.assign(node_num,unordered_map<int,vector<int>>());

                    getline(inFile,line);
                    re1.clear();
                    boost::split(re1,line,boost::is_any_of(" \t"));
                    while(re1[0]=="c" || re1[0]=="p"){
                        getline(inFile,line);
                        re1.clear();
                        boost::split(re1,line,boost::is_any_of(" \t"));
                    }

                    getline(inFile2,line2);
                    re2.clear();
                    boost::split(re2,line2,boost::is_any_of(" \t"));
                    while(re2[0]=="c" || re2[0]=="p"){
                        getline(inFile2,line2);
                        re2.clear();
                        boost::split(re2,line2,boost::is_any_of(" \t"));
                    }


                    if(re1[0]=="a" && re1.size()>=4 && re2[0]=="a" && re2.size()>=4){
                        ID1=stoi(re1[1]), ID2=stoi(re1[2]), weight=stoi(re1[3]);
                        if(ID1!=stoi(re2[1]) || ID2!=stoi(re2[2])){
                            cout<<"Inconsistent! "<<ID1<<" "<<ID2<<" "<<stoi(re2[1])<<" "<<stoi(re2[2])<<endl;
                            exit(1);
                        }else{
                            weight2=stoi(re2[3]);
                            ID1-=1; ID2-=1;
                            if(ID1>=0 && ID1<node_num && ID2>=0 && ID2<node_num && ID1!=ID2){
                                if(pairGraph[ID1].find(ID2)==pairGraph[ID1].end()){//if not found
                                    tempCri[0]=weight, tempCri[1]=weight2;
                                    pairGraph[ID1].insert({ID2,tempCri});
                                    pairGraph[ID2].insert({ID1,tempCri});
                                    base_max = max(base_max,weight);
                                    base_min = min(base_min,weight);
                                    edgeNum+=2;
                                    set_A.insert(ID1); set_A.insert(ID2);
                                    degrees[ID1]++; degrees[ID2]++;
                                }

                                if(weight<=0){
                                    cout<<"Negative edge weight "<<ID1<<" "<<ID2<<" "<<weight<<endl;
                                }
                            }else{
                                if(ID1!=ID2){
                                    cout<<"Graph data is wrong! "<<ID1<<" "<<ID2<<" "<<weight<<endl;
                                    exit(1);
                                }
                            }
                        }
                    }

                    while(getline(inFile,line) && getline(inFile2,line2)){
                        if(line.empty() || line2.empty()){
                            if(!line.empty()){
                                cout<<"Line 1 is not empty! "<<line<<endl;
                            }
                            if(!line2.empty()){
                                cout<<"Line 2 is not empty! "<<line2<<endl;
                            }
                            continue;
                        }

                        re1.clear();
                        boost::split(re1,line,boost::is_any_of(" \t"));
                        re2.clear();
                        boost::split(re2,line2,boost::is_any_of(" \t"));
                        if(re1[0]=="a" && re1.size()>=4 && re2[0]=="a" && re2.size()>=4){
                            ID1=stoi(re1[1]), ID2=stoi(re1[2]), weight=stoi(re1[3]);
                            if(ID1!=stoi(re2[1]) || ID2!=stoi(re2[2])){
                                cout<<"Inconsistent! "<<ID1<<" "<<ID2<<" "<<stoi(re2[1])<<" "<<stoi(re2[2])<<endl;
                                exit(1);
                            }else{
                                weight2=stoi(re2[3]);
                                ID1-=1; ID2-=1;
                                if(ID1>=0 && ID1<node_num && ID2>=0 && ID2<node_num && ID1!=ID2){
                                    if(pairGraph[ID1].find(ID2)==pairGraph[ID1].end()){//if not found
                                        tempCri[0]=weight, tempCri[1]=weight2;
                                        pairGraph[ID1].insert({ID2,tempCri});
                                        pairGraph[ID2].insert({ID1,tempCri});
                                        base_max = max(base_max,weight);
                                        base_min = min(base_min,weight);
                                        edgeNum+=2;
                                        set_A.insert(ID1); set_A.insert(ID2);
                                        degrees[ID1]++; degrees[ID2]++;
                                    }

                                    if(weight<=0){
                                        cout<<"Negative edge weight "<<ID1<<" "<<ID2<<" "<<weight<<endl;
                                    }
                                }else{
                                    if(ID1!=ID2){
                                        cout<<"Graph data is wrong! "<<ID1<<" "<<ID2<<" "<<weight<<endl;
                                        exit(1);
                                    }
                                }
                            }

                        }
                        else{
                            cout<<"Wrong line. "<< line<<" "<<line2<<endl;
                            exit(1);
                        }

                    }

                    inFile.close();
                    tt.stop();
                    cout << "The time for data loading: " << tt.GetRuntime() << " s." << endl;
                    temp_num = set_A.size();
                    if(temp_num != node_num){
                        cout << "\n!!! The valid node number is " << temp_num <<" , the valid edge number is "<<edgeNum<<endl;
                    }

                    max_degree = *max_element(degrees.begin(),degrees.end());
                    min_degree = *min_element(degrees.begin(),degrees.end());
                    cout << "Maximal degree: " << max_degree << " \tMinimal degree: " << min_degree <<endl;


                    pair<int,unsigned long long> nums = DFS_CC(pairGraph, set_A, set_LCC, node_num);
                    node_num=nums.first, edge_num=nums.second;
                    cout<<"Node number: "<<node_num<<" ; edge number: "<<edge_num<<endl;

                    WriteGraph(w_file, set_LCC, pairGraph, 0);
                    WriteGraph(w_file2, set_LCC, pairGraph, 1);
//                exit(0);
                }
                else{//if the files exist
                    cout<<"The basic criterion data exists!"<<endl;
                    string line;
                    getline(inFilew,line);
                    vector<string> re1;
                    boost::split(re1,line,boost::is_any_of(" \t"));
                    if(re1.size()==2){
                        node_num=stoi(re1[0]), edge_num=stoul(re1[1]);
                    }else{
                        cout<<"Wrong input! re.size: "<<re1.size()<<". " <<line<<endl;
                        exit(1);
                    }
                    cout<<"Node number: "<<node_num<<" ; edge number: "<<edge_num<<endl;
                    degrees.assign(node_num,0);
                    pairGraph.assign(node_num,unordered_map<int,vector<int>>());

                    while(getline(inFilew,line)){
                        if(line.empty())
                            continue;
                        re1.clear();
                        boost::split(re1,line,boost::is_any_of(" \t"));
                        if(re1.size()==3){
                            ID1=stoi(re1[0]), ID2=stoi(re1[1]), weight=stoi(re1[2]);
                        }
                        else{
                            cout<<"Wrong line. "<< line<<endl;
                            exit(1);
                        }


                        if(ID1>=0 && ID1<node_num && ID2>=0 && ID2<node_num && ID1!=ID2 && weight>0){
                            if(pairGraph[ID1].find(ID2)==pairGraph[ID1].end()){//if not found
                                tempCri[0]=weight;
                                pairGraph[ID1].insert({ID2,tempCri});
                                pairGraph[ID2].insert({ID1,tempCri});
                                base_max = max(base_max,weight);
                                base_min = min(base_min,weight);
                                edgeNum+=2;
                                set_LCC.insert(ID1); set_LCC.insert(ID2);
                                degrees[ID1]++; degrees[ID2]++;
                            }

                        }else{
                            if(ID1!=ID2){
                                cout<<"Graph data is wrong! "<<ID1<<" "<<ID2<<" "<<weight<<endl;
                                exit(1);
                            }
                        }
                    }
                    inFilew.close();


                    /// For Time criterion
                    getline(inFilew2,line);
                    re1.clear();
                    boost::split(re1,line,boost::is_any_of(" \t"));
                    if(re1.size()==2){
//                    node_num=stoi(re1[0]), edge_num=stoul(re1[1]);
                    }else{
                        cout<<"Wrong input! re.size: "<<re1.size()<<". " <<line<<endl;
                        exit(1);
                    }

                    while(getline(inFilew2,line)){
                        if(line.empty())
                            continue;
                        re1.clear();
                        boost::split(re1,line,boost::is_any_of(" \t"));
                        if(re1.size()==3){
                            ID1=stoi(re1[0]), ID2=stoi(re1[1]), weight=stoi(re1[2]);
                        }
                        else{
                            cout<<"Wrong line. "<< line<<endl;
                            exit(1);
                        }

                        if(ID1>=0 && ID1<node_num && ID2>=0 && ID2<node_num && ID1!=ID2 && weight>0){
                            if(pairGraph[ID1].find(ID2)!=pairGraph[ID1].end()){//if found
                                tempCri[1]=weight;
                                pairGraph[ID1][ID2][1]=weight;
                                pairGraph[ID2][ID1][1]=weight;
                                edgeNum+=2;
                            }else{
                                cout<<"Not found! "<<ID1<<" "<<ID2<<endl;
                                exit(1);
                            }

                        }else{
                            if(ID1!=ID2){
                                cout<<"Graph data is wrong! "<<ID1<<" "<<ID2<<" "<<weight<<endl;
                                exit(1);
                            }
                        }
                    }
                    inFilew.close();

                    max_degree = *max_element(degrees.begin(),degrees.end());
                    min_degree = *min_element(degrees.begin(),degrees.end());
                    cout << "Maximal degree: " << max_degree << " \tMinimal degree: " << min_degree <<endl;

                }

            }

            /// generate extra criteria
            MC_Generation(NUM_OF_CRITERIA,set_LCC,pairGraph);
        }

    }
    //function of writing the graph with single criterion into disk
    void Graph_pre::WriteGraph(string filename, set<int>& set_LCC, vector<unordered_map<int,vector<int>>>& NeighborMap, int cri_i) {
        ofstream outFile(filename, ios::out);
        if (!outFile) {
            cout << "Write File opening failed. " << filename<<endl;
            assert(outFile);
            exit(1);
        }
        cout<<"Writing graph into disk..."<<endl;
        Timer tt;
        tt.start();
        outFile << node_num <<" "<<edge_num<< endl;

        bool flag_double = true;
        if(edge_num>1000000000){
            flag_double=false;
            cout<<"Write edge once!"<<endl;
        }

        int ID=0;
        map<int,int> IDMap;//map from old ID to new ID

        for(auto it=set_LCC.begin();it!=set_LCC.end();++it){
            IDMap.insert({*it,ID});
            ++ID;
        }
        assert(ID==node_num);

        for(auto it=set_LCC.begin();it!=set_LCC.end();++it){
            ID=*it;
            for(auto it2=NeighborMap[ID].begin();it2!=NeighborMap[ID].end();++it2){
                if(IDMap.find(it2->first) != IDMap.end()){
                    if(flag_double || IDMap[ID]<IDMap[it2->first]){
                        outFile<<IDMap[ID]<<" "<<IDMap[it2->first]<<" "<<it2->second[cri_i]<<endl;
                    }

                }else{
                    cout<<"Wrong! "<<ID<<" "<<it2->first<<" "<<it2->second[cri_i]<<endl;
                }

            }
        }

        outFile.close();
        tt.stop();
        cout<<"Write Done. Time: "<<tt.GetRuntime()<<" s."<< endl;

    }
    //Function for data generation and writing, new version
    void Graph_pre::MC_Generation(int num, set<int>& set_LCC, vector<unordered_map<int,vector<int>>>& NeighborMap){
        string w_file = DataPath + dataset + "/" + dataset + ".MCEdges";
        int base;
        int temp_num, temp_int;
        int ID1,ID2;

        ifstream inFile(w_file,ios::in);
        if(inFile.is_open()){//if exists
            cout<<"The MCEdges file exists!"<<endl;
            num_criteria = 5;//update the number of criteria
        }
        else{//if not exists
            cout<<"The MCEdges file does not exist!"<<endl;
            //normal distribution generator
            default_random_engine generator(0);
            normal_distribution<double> norm_distribution(0.0, 0.1);

            mc_criteria.emplace_back("Distance");
            if(num_criteria == 1){
                mc_criteria.emplace_back("Corr");
            }else if(num_criteria == 2){
                mc_criteria.emplace_back("Time");
            }
            mc_criteria.emplace_back("Corr");
//        mc_criteria.emplace_back("Corr");
            mc_criteria.emplace_back("Anti");
            mc_criteria.emplace_back("Indep");

            max_weights.assign(NUM_OF_CRITERIA,0);
            min_weights.assign(NUM_OF_CRITERIA,INF);

            ofstream outFile(w_file, ios::out);
            if (!outFile) {
                cout << "Write File opening failed. " << w_file<<endl;
                assert(outFile);
                exit(1);
            }
            cout<<"Writing muti-criteria graph into disk..."<<endl;
            Timer tt;
            tt.start();
            outFile << node_num <<" "<<edge_num<< endl;
            outFile<<NUM_OF_CRITERIA;

            for (int i = 0; i < mc_criteria.size(); ++i) {
                if(i == mc_criteria.size()-1)
                    outFile<<" "<<mc_criteria[i]<<endl;
                else
                    outFile<<" "<<mc_criteria[i];
            }


            bool flag_double = true;
            if(edge_num>1000000000){
                flag_double=false;
                cout<<"Write edge once!"<<endl;
            }

            int ID=0;
            map<int,int> IDMap;//map from old ID to new ID

            for(auto it=set_LCC.begin();it!=set_LCC.end();++it){
                IDMap.insert({*it,ID});
                ++ID;
            }
            assert(ID==node_num);

            vector<int> tempW(NUM_OF_CRITERIA,0);

            for(auto it=set_LCC.begin();it!=set_LCC.end();++it){
                ID=*it;
                for(auto it2=NeighborMap[ID].begin();it2!=NeighborMap[ID].end();++it2){
                    if(IDMap.find(it2->first) != IDMap.end()){
                        tempW[0] = it2->second[0];
                        base = tempW[0];
                        max_weights[0] = max(max_weights[0],temp_int);
                        min_weights[0] = min(min_weights[0],temp_int);
                        //criterion 2: correlated
                        if(num_criteria==1){
                            temp_int = max(int((2 + norm_distribution(generator)) * base), base_min);
                            tempW[1]=temp_int;
                            max_weights[1] = max(max_weights[1],temp_int);
                            min_weights[1] = min(min_weights[1],temp_int);
                        }else if(num_criteria==2){
                            tempW[1]=it2->second[1];
                            max_weights[1] = max(max_weights[1],temp_int);
                            min_weights[1] = min(min_weights[1],temp_int);
                        }
                        //criterion 3: correlated
                        temp_int = max(int((1 + norm_distribution(generator)) * base), base_min);
                        tempW[2] = temp_int;
                        max_weights[2] = max(max_weights[2],temp_int);
                        min_weights[2] = min(min_weights[2],temp_int);
                        //criterion 4: anti-correlated
                        temp_int = max(int(1000*(base_max - (1+norm_distribution(generator))*base) / base_max), base_min);
                        tempW[3] = temp_int;
                        max_weights[3] = max(max_weights[3],temp_int);
                        min_weights[3] = min(min_weights[3],temp_int);
                        //criterion 5: independent, from base_min to 1000
                        temp_int = base_min + int(1000*(rand()/double(RAND_MAX)) * (base_max - base_min) / base_max);
                        tempW[4] = temp_int;
                        max_weights[4] = max(max_weights[4],temp_int);
                        min_weights[4] = min(min_weights[4],temp_int);

                        if(IDMap[ID]<IDMap[it2->first]){
                            outFile<<IDMap[ID]<<" "<<IDMap[it2->first];
                            for(int i=0;i<tempW.size();++i){
                                outFile<<" "<<tempW[i];
                            }
                            outFile<<endl;

                            if(flag_double){
                                outFile<<IDMap[it2->first]<<" "<<IDMap[ID];
                                for(int i=0;i<tempW.size();++i){
                                    outFile<<" "<<tempW[i];
                                }
                                outFile<<endl;
                            }
                        }

                    }else{
                        cout<<"Wrong! "<<ID<<" "<<it2->first<<" "<<it2->second[0]<<endl;
                    }

                }
            }

            outFile.close();
            tt.stop();
            cout<<"Write Done. Time: "<<tt.GetRuntime()<<" s."<< endl;
            NeighborMap.clear();
            num_criteria = 5;//update the number of criteria
        }

    }

    void Graph_pre::PartitionPreprocess(string parMethod) {
        string r_file = DataPath + dataset + "/" + dataset + ".MCEdges";
        num_criteria = 5;
        ReadMCEdges(r_file);

        if(parMethod=="KaHyPar"){
            cout<<"Partition method: "<<parMethod<<endl;
            KaHyParConvertInput();
        }else if(parMethod == "Metis"){
            cout<<"Partition method: "<<parMethod<<endl;
            METISConvertInput();
        }else{
            cout<<"Wrong partition method!"<<endl;
            exit(1);
        }

    }

    void Graph_pre::KaHyParConvertInput(){
        string w_file = DataPath + dataset + "/" +strategy+"/"+ dataset + ".mcKaHyPar";
        ofstream outFile(w_file,ios::out);
        if(!outFile.is_open()){
            cout<<"Cannot open "<<w_file<<endl;
            exit(1);
        }

        bool flag_double = true;
        if(edge_num>1000000000){
            flag_double=false;
            cout<<"Read edge once!"<<endl;
        }

        int aggregate=0;
        if(strategy == "ave"){
            aggregate=1;
            cout<<"Edge aggregation strategy: "<<strategy<<endl;
        }else if(strategy == "pri"){
            aggregate=2;
            cout<<"Edge aggregation strategy: "<<strategy<<endl;
        }else if(strategy == "wave"){
            aggregate=3;
            cout<<"Edge aggregation strategy: "<<strategy<<endl;
        }else{
            cout<<"Wrong edge aggregation strategy! "<<strategy<<endl;
            exit(1);
        }

        outFile<<node_num<<" "<<edge_num<<" 1"<<endl;
        int temp_id,temp_w;
        for(int id=0;id<node_num;++id){
            for (auto it = MCEdges[id].begin(); it != MCEdges[id].end(); ++it) {
                temp_id = it->u;
                if(aggregate == 1){
                    temp_w = ceil(accumulate(it->w.begin(),it->w.end(),0.0)/NUM_OF_CRITERIA);
                }else if(aggregate == 2){
                    temp_w = it->w[sc_i];
                }else if(aggregate == 3){
                    int ave_w = 0;
                    for(int i=0;i<num_criteria;++i){
                        ave_w += zeta[i] * it->w[i];
                    }
                    temp_w = ceil(ave_w);
                }else{
                    cout<<"Wrong edge aggregation strategy! "<<strategy<<endl;
                    exit(1);
                }
                outFile<<temp_w<<" "<<id+1<<" "<<temp_id+1<<endl;
                if(!flag_double){
                    outFile<<temp_w<<" "<<temp_id+1<<" "<<id+1<<endl;
                }
            }
        }
        outFile.close();
    }


    void Graph_pre::METISConvertInput(){
        string w_file = DataPath + dataset + "/" +strategy+"/"+ dataset + ".mcMetis";
        ofstream outFile(w_file,ios::out);
        if(!outFile.is_open()){
            cout<<"Cannot open "<<w_file<<endl;
            exit(1);
        }

        bool flag_double = true;
        if(edge_num>1000000000){
            flag_double=false;
            cout<<"Read edge once!"<<endl;
        }

        int aggregate=0;
        if(strategy == "ave"){
            aggregate=1;
            cout<<"Edge aggregation strategy: "<<strategy<<endl;
        }else if(strategy == "pri"){
            aggregate=2;
            cout<<"Edge aggregation strategy: "<<strategy<<endl;
        }else if(strategy == "wave"){
            aggregate=3;
            cout<<"Edge aggregation strategy: "<<strategy<<endl;
        }else{
            cout<<"Wrong edge aggregation strategy! "<<strategy<<endl;
            exit(1);
        }

        outFile<<node_num<<" "<<edge_num/2<<" 001"<<endl;
        int temp_id,temp_w;
        for(int id=0;id<node_num;++id){
            for (auto it = MCEdges[id].begin(); it != MCEdges[id].end(); ++it) {
                temp_id = it->u;
                if(aggregate == 1){
                    temp_w = ceil(accumulate(it->w.begin(),it->w.end(),0.0)/NUM_OF_CRITERIA);
                }else if(aggregate == 2){
                    temp_w = it->w[sc_i];
                }else if(aggregate == 3){
                    int ave_w = 0;
                    for(int i=0;i<num_criteria;++i){
                        ave_w += zeta[i] * it->w[i];
                    }
                    temp_w = ceil(ave_w);
                }else{
                    cout<<"Wrong edge aggregation strategy! "<<strategy<<endl;
                    exit(1);
                }
                outFile<<temp_id+1<<" "<<temp_w<<" ";
            }
            outFile<<endl;
        }
        outFile.close();
    }

    void Graph_pre::PartitionPostprocess(string parMethod, int partitionNum){
        string r_file = DataPath + dataset + "/" + dataset + ".MCEdges";
        num_criteria = 5;
        ReadMCEdges(r_file);
        int partition_i=-1;

        cluster_to_node.assign(partitionNum,vector<NodeId>());
        string line;
        if(partMethod=="DC"){
            cout<<"Determine clustering."<<endl;
            DeterClustering(false);
        }
        else if(parMethod=="Metis"){
            cout<<"Postprocess of METIS partition result."<<endl;
            string filename=DataPath + dataset + "/" + strategy+"/" + dataset+".mcMetis.part."+to_string(partitionNum);
            ifstream inFile(filename, ios::in);
            if(!inFile.is_open()){
                cout<<"Cannot open "<<filename<<endl;
                exit(1);
            }

            NodeId id=0;
            while (getline(inFile, line)){
                if(line.empty())
                    continue;
                partition_i=stoi(line);
                if(partition_i>=0 && partition_i<partitionNum){
                    cluster_to_node[partition_i].emplace_back(id);
                }else{
                    cout<<"Wrong partition id!"<<partition_i<<endl;
                    exit(1);
                }
                ++id;
            }
            inFile.close();
        }
        else if(parMethod=="KaHyPar"){
            cout<<"Postprocess of METIS partition result."<<endl;
        }else{
            cout<<"Wrong partition method!"<<endl;
            exit(1);
        }
        // convert to desired output
        for(int pid=0;pid<partitionNum;++pid){
            WritePartition(pid, cluster_to_node[pid], true);
        }
        //write partition_info to disk
        WritePartitionInfo(cluster_to_node);
        cout << "Done." << endl;

    }

    void Graph_pre::GraphPartitioning(){
        string r_file = DataPath + dataset + "/" + dataset + ".MCEdges";
        ReadMCEdges(r_file);

        DeterClustering(false);

    }
    //function of reading graph data into MCEdges
    void Graph_pre::ReadMCEdges(const string& filename){
        Timer tt;
        tt.start();
        string temp_str;
        int ID1, ID2, weight;
        MCEdge temp_edge;
        vector<int> temp_w;
        int temp_id = -1;
        uint temp_degree = 1;
        unordered_set<int> set_A; set<int> set_LCC;
        set_A.clear();

        max_weights.assign(NUM_OF_CRITERIA,INT32_MIN);
        min_weights.assign(NUM_OF_CRITERIA,INT32_MAX);

        //Open file
        ifstream inFile(filename, ios::in);
        if (!inFile) {
            cout << "File opening failed. " << filename<< endl;
            exit(1);
        }

        cout << "MCEdges Data loading..." << endl;
        string line;
        getline(inFile,line);
        vector<string> re1;
        boost::split(re1,line,boost::is_any_of(" \t"));
        if(re1.size()==2){
            node_num=stoi(re1[0]), edge_num=stoul(re1[1]);
        }else{
            cout<<"Wrong input! re.size: "<<re1.size()<<". " <<line<<endl;
            exit(1);
        }
        cout<<"Node number: "<<node_num<<" ; edge number: "<<edge_num<<endl;

        MCEdges.assign(node_num,vector<MCEdge>());

        bool flag_double = true;
        if(edge_num>1000000000){
            flag_double=false;
            cout<<"Read edge once!"<<endl;
        }

        getline(inFile,line);
        re1.clear();
        boost::split(re1,line,boost::is_any_of(" \t"));
        mc_criteria.clear();
        int tempNum_criteria=stoi(re1[0]);
        if(re1.size()==tempNum_criteria+1){
            for(int i=0;i<tempNum_criteria;++i){
                temp_str = re1[i+1];
                cout << temp_str << " " ;
                mc_criteria.push_back(temp_str);
            }
        }else{
            cout<<"Wrong input! re.size: "<<re1.size()<<". " <<line<<endl;
            exit(1);
        }
        cout<<endl;

        num_criteria=tempNum_criteria;
        tempNum_criteria+=2;

        while(getline(inFile,line)){
            if(line.empty())
                continue;
            re1.clear();
            boost::split(re1,line,boost::is_any_of(" \t"));
            if(re1.size()==tempNum_criteria){
                ID1=stoi(re1[0]), ID2=stoi(re1[1]);
            }
            else{
                cout<<"Wrong line. "<< line<<endl;
                exit(1);
            }


            if(ID1>=0 && ID1<node_num && ID2>=0 && ID2<node_num && ID1!=ID2){
                set_A.insert(ID1), set_A.insert(ID2);
                temp_w.clear();
                for(int i=0;i<num_criteria;++i) {//get edge weights
                    weight=stoi(re1[i+2]);
                    max_weights[i] = max(max_weights[i],weight);
                    min_weights[i] = min(min_weights[i],weight);
                    temp_w.push_back(weight);
                }
                MCEdges[ID1].emplace_back(MCEdge(ID2,temp_w));
                if(!flag_double){
                    MCEdges[ID2].emplace_back(MCEdge(ID1,temp_w));
                }

            }else{
                if(ID1!=ID2){
                    cout<<"Graph data is wrong! "<<ID1<<" "<<ID2<<" "<<weight<<endl;
                    exit(1);
                }
            }
        }
        inFile.close();

//        pair<int,unsigned long long> nums = DFS_CC_MC(MCEdges, set_A, set_LCC, node_num);
//        if(node_num!=nums.first){
//            cout<<"Wrong for node_num! "<<node_num<<" "<<nums.first;
//            node_num=nums.first, edge_num=nums.second;
//            cout<<"Node number: "<<node_num<<" ; edge number: "<<edge_num<<endl;
//            exit(1);
//        }


        vector<int> degrees(node_num,0);
        ave_degree=0;
        for(int i=0;i<node_num;++i){
            degrees[i]=MCEdges[i].size();
            ave_degree+=degrees[i];
        }
        max_degree = *max_element(degrees.begin(),degrees.end());
        ave_degree /= node_num;
        cout<<"Maximal degree: "<<max_degree<<" ; average degree: "<<int(ave_degree)<<endl;

        tt.stop();
        cout << "The time used for data loading:" << tt.GetRuntime() << " s." << endl;
        inFile.close();
    }

    void Graph_pre::IDMinusOne(string filename) {
        Timer tt;
        tt.start();
        string temp_str;
        int ID1, ID2, weight;
        MCEdge temp_edge;
        vector<int> temp_w;
        int temp_id = -1;
        uint temp_degree = 1;


        //Open file
        ifstream inFile(filename, ios::in);
        if (!inFile) {
            cout << "File opening failed. " << filename<< endl;
            exit(1);
        }

        ofstream outFile(filename+".new", ios::out);
        if(!outFile.is_open()){
            cout << "File opening failed. " << filename+".new"<< endl;
            exit(1);
        }

        cout << "MCEdges Data loading..." << endl;
        string line;
        getline(inFile,line);
        vector<string> re1;
        boost::split(re1,line,boost::is_any_of(" \t"));
        if(re1.size()==2){
            node_num=stoi(re1[0]), edge_num=stoul(re1[1]);
        }else{
            cout<<"Wrong input! re.size: "<<re1.size()<<". " <<line<<endl;
            exit(1);
        }
        cout<<"Node number: "<<node_num<<" ; edge number: "<<edge_num<<endl;

        outFile<<node_num<<" "<<edge_num<<endl;


        getline(inFile,line);
        re1.clear();
        boost::split(re1,line,boost::is_any_of(" \t"));
        mc_criteria.clear();
        int tempNum_criteria=stoi(re1[0]);
        outFile<<tempNum_criteria;

        if(re1.size()==tempNum_criteria+1){
            for(int i=0;i<tempNum_criteria;++i){
                temp_str = re1[i+1];
                cout << temp_str << " " ;
                outFile<<" "<<temp_str;
                mc_criteria.push_back(temp_str);
            }
        }else{
            cout<<"Wrong input! re.size: "<<re1.size()<<". " <<line<<endl;
            exit(1);
        }
        cout<<endl;
        outFile<<endl;

        num_criteria=tempNum_criteria;
        tempNum_criteria+=2;

        while(getline(inFile,line)){
            if(line.empty())
                continue;
            re1.clear();
            boost::split(re1,line,boost::is_any_of(" \t"));
            if(re1.size()==tempNum_criteria){
                ID1=stoi(re1[0]), ID2=stoi(re1[1]);
            }
            else{
                cout<<"Wrong line. "<< line<<endl;
                exit(1);
            }

            outFile<<ID1-1<<" "<<ID2-1;
//            outFile<<ID1-1<<" "<<ID2-1<<" "<<stoi(re1[2]);
            for(int i=0;i<num_criteria;++i) {//get edge weights
                weight=stoi(re1[i+2]);
                if(i>0){
                    continue;
                }
                if(weight==100000000){
                    outFile<<" "<<INF100;
                }else{
                    outFile<<" "<<weight;
                }

            }
            outFile<<endl;
        }
        inFile.close();
        outFile.close();
    }

    //function of checking the connectivity
    template <class T>
    pair<int, unsigned long long> Graph_pre::DFS_CC(T & Edges, unordered_set<int> & set_A, set<int> & set_LCC, int nodenum) {
        /// DFS for connected component
        stack<int> stack_A;
//    set<int> set_A;//nodes waiting for visited
        unordered_set<int> set_B;//nodes visited for current component
        set_B.clear();
        int item_id,temp_id;
        vector<bool> flag_visited(nodenum,false);
        bool flag_finish = false;
        unsigned long long temp_num = 0;
        int component_i = 0;
        pair<unordered_set<int>,unsigned long long> LCC;
        vector<int> CCs;//the vertex size of each connected component

//    for(int i=0;i<nodenum;++i){
//        set_A.insert(i);
//    }
        int seed = *set_A.begin();
        stack_A.push(seed);
        set_A.erase(seed);
        set_B.insert(seed);
        flag_visited[seed] = true;
        //Get the connected components by DFS
        while(!set_A.empty()) {//if not finish
            temp_num = 0;
            while (!stack_A.empty()) {
                item_id = stack_A.top();
                stack_A.pop();
                for (auto it = Edges[item_id].begin(); it != Edges[item_id].end(); ++it) {
                    temp_id = it->first;
                    temp_num += 1;
                    if (!flag_visited[temp_id]) {//if not visited
                        stack_A.push(temp_id);
                        set_A.erase(temp_id);
                        set_B.insert(temp_id);
                        flag_visited[temp_id] = true;
                    }
                }
            }
            if (set_B.size() > LCC.first.size()) {
                LCC.first.clear();
                LCC.first = set_B;
                LCC.second = temp_num;// /2
            }
            assert(!set_B.empty());
            CCs.push_back(set_B.size());
//        if(!set_B.empty() && set_B.size() < mcc.first.size()){
//            cout<<"Smaller connected component with vertex size "<<set_B.size()<<": ";
//            for(auto it=set_B.begin();it!=set_B.end();++it){
//                cout<<*it<<" ";
//            }
//            cout<<"; degree: ";
//            for(auto it=set_B.begin();it!=set_B.end();++it){
//                cout<<Edges[*it].size()<<" ";
//            }
//            cout<<endl;
//        }
            ++component_i;
            set_B.clear();
            if (!set_A.empty()) {
                stack_A.push(*set_A.begin());
                set_B.insert(*set_A.begin());
                flag_visited[*set_A.begin()] = true;
                set_A.erase(*set_A.begin());
            } else {
                break;
            }
        }
        if(component_i==1){
            cout<<"This graph has only one connected component. ";
            cout<<"Nodes size of graph: "<< LCC.first.size() << " ; ";
            cout<<"Edges size of graph: "<< LCC.second << endl;
        }else{
            cout<<"!!! This graph has "<< component_i <<" connected component!"<<endl;
            cout<<"Nodes size of the largest connected component is: "<<LCC.first.size()<<endl;
            cout<<"Edges size of the largest connected component is: "<<LCC.second<<endl;
        }
        for(auto it=LCC.first.begin();it!=LCC.first.end();++it){
            set_LCC.insert(*it);
        }
        std::sort(CCs.begin(), CCs.end());

        return make_pair(LCC.first.size(),LCC.second);
    }

    pair<int, unsigned long long> Graph_pre::DFS_CC_MC(vector<vector<MCEdge>> & Edges, unordered_set<int> & set_A, set<int> & set_LCC, int nodenum) {
        /// DFS for connected component
        stack<int> stack_A;
//    set<int> set_A;//nodes waiting for visited
        unordered_set<int> set_B;//nodes visited for current component
        set_B.clear();
        int item_id,temp_id;
        vector<bool> flag_visited(nodenum,false);
        bool flag_finish = false;
        unsigned long long temp_num = 0;
        int component_i = 0;
        pair<unordered_set<int>,unsigned long long> LCC;
        vector<int> CCs;//the vertex size of each connected component

//    for(int i=0;i<nodenum;++i){
//        set_A.insert(i);
//    }
        int seed = *set_A.begin();
        stack_A.push(seed);
        set_A.erase(seed);
        set_B.insert(seed);
        flag_visited[seed] = true;
        //Get the connected components by DFS
        while(!set_A.empty()) {//if not finish
            temp_num = 0;
            while (!stack_A.empty()) {
                item_id = stack_A.top();
                stack_A.pop();
                for (auto it = Edges[item_id].begin(); it != Edges[item_id].end(); ++it) {
                    temp_id = it->u;
                    temp_num += 1;
                    if (!flag_visited[temp_id]) {//if not visited
                        stack_A.push(temp_id);
                        set_A.erase(temp_id);
                        set_B.insert(temp_id);
                        flag_visited[temp_id] = true;
                    }
                }
            }
            if (set_B.size() > LCC.first.size()) {
                LCC.first.clear();
                LCC.first = set_B;
                LCC.second = temp_num;// /2
            }
            assert(!set_B.empty());
            CCs.push_back(set_B.size());
//        if(!set_B.empty() && set_B.size() < mcc.first.size()){
//            cout<<"Smaller connected component with vertex size "<<set_B.size()<<": ";
//            for(auto it=set_B.begin();it!=set_B.end();++it){
//                cout<<*it<<" ";
//            }
//            cout<<"; degree: ";
//            for(auto it=set_B.begin();it!=set_B.end();++it){
//                cout<<Edges[*it].size()<<" ";
//            }
//            cout<<endl;
//        }
            ++component_i;
            set_B.clear();
            if (!set_A.empty()) {
                stack_A.push(*set_A.begin());
                set_B.insert(*set_A.begin());
                flag_visited[*set_A.begin()] = true;
                set_A.erase(*set_A.begin());
            } else {
                break;
            }
        }
        if(component_i==1){
            cout<<"This graph has only one connected component. ";
            cout<<"Nodes size of graph: "<< LCC.first.size() << " ; ";
            cout<<"Edges size of graph: "<< LCC.second << endl;
        }else{
            cout<<"!!! This graph has "<< component_i <<" connected component!"<<endl;
            cout<<"Nodes size of the largest connected component is: "<<LCC.first.size()<<endl;
            cout<<"Edges size of the largest connected component is: "<<LCC.second<<endl;
        }
        for(auto it=LCC.first.begin();it!=LCC.first.end();++it){
            set_LCC.insert(*it);
        }
        std::sort(CCs.begin(), CCs.end());

        return make_pair(LCC.first.size(),LCC.second);
    }

    //Affiliate variables and functions
    typedef int vertex_id;
    const long int BLK_SZ = 4194304*8;  //1024*1024*4 bytes(i.e. 4MB) sizeof blk for I/O
    const int SZ_VERTEX = sizeof(int);
    const int VERTEX_PER_BLK = BLK_SZ/SZ_VERTEX;//number of int readable, 1024*1024*8
    float ioTime; //time of io
    clock_t read_time = 0;
    clock_t write_time = 0;
    long no_fread = 0;
    long no_fwrite = 0;

    int fread_c(void * ptr, size_t size, size_t count,FILE * stream){

        Timer rt;
        rt.start();

        clock_t fr_bg = clock();
        int returnV = fread (ptr,size,count,stream);

        clock_t fr_end = clock();

        read_time += (fr_end-fr_bg);

        rt.stop();
        ioTime += rt.GetRuntime();
        ++no_fread;

        return returnV;
    }
    int fwrite_c ( const void * ptr, size_t size, size_t count, FILE * stream ){

        Timer rt;
        rt.start();

        clock_t fw_bg = clock();

        int returnV = fwrite(ptr,size,count,stream);

        clock_t fw_end = clock();

        rt.stop();
        ioTime += rt.GetRuntime();

        write_time += (fw_end - fw_bg);

        ++no_fwrite;

        return returnV;
    }
    //function for reading single criterion data from disk to MCEdgesMap
    void Graph_pre::Read_ToMCEdgesMap(const string& filename,bool ifBase){
        Timer tt;
        tt.start();
        string line_symbol;
        string temp_str;
        int num_line = 0;
        int num_show = 6;//number of lines to be printed
        int ID1, ID2, weight;
        vector<EdgeWeight> temp_w;
        int temp_max = INT32_MIN;
        int temp_min = INT32_MAX;

        //Open file
        ifstream inFile(filename, ios::in);
        if (!inFile) {
            cout << "File opening failed." << endl;
            assert(inFile);
            exit(1);
        }
        if(ifBase){
            cout << dataset<<" "<< criteria<<" graph edges Data loading..." << endl;
            while (inFile) {//read each line to lineStr
                inFile >> line_symbol;
                if (line_symbol == "p") {//read graph basic information
                    inFile >> temp_str;
                    if(temp_str=="sp"){
                        inFile >> this->node_num >> this->edge_num;
                        cout<<"Node number: "<<node_num<<"\tEdge number: "<<edge_num<<endl;
                        if(ifBase){
                            MCEdgesMap.assign(node_num,map<NodeId,MCEdge>());
                        }
                        cout << "Some edges:" << endl;
                    }else{
                        getline(inFile, temp_str);
                    }

                    //break;
                }
                else if (line_symbol == "a") {//read graph data
                    inFile >> ID1 >> ID2 >> weight;
//                if(invalid_edges.find(make_pair(ID1,ID2))!=invalid_edges.end()){//if found in invalid set
//                    continue;
//                }
                    temp_max = max(temp_max,weight);
                    temp_min = min(temp_min,weight);

                    temp_w.clear();
                    temp_w.emplace_back(weight);
                    temp_w.emplace_back(0);temp_w.emplace_back(0);temp_w.emplace_back(0);temp_w.emplace_back(0);
//                    MCEdgesMap[ID1 - 1].insert({ID2-1,MCEdge(ID2 - 1, temp_w)});
                    MCEdgesMap[ID1].insert({ID2,MCEdge(ID2, temp_w)});
                    if (num_line < num_show) {
                        cout << ID1 << "\t" << ID2 << "\t" << weight << endl;
                        ++num_line;
                    }
                    else if (num_line == num_show) {
                        cout << "..." << endl;
                        ++num_line;
                    }

                }else{
                    getline(inFile, temp_str);
                }
            }

        }else{
            cout << dataset<<" Time graph edges Data loading..." << endl;
            while (inFile) {//read each line to lineStr
                inFile >> line_symbol;
                if (line_symbol == "p") {//read graph basic information
                    inFile >> temp_str;
                    if(temp_str=="sp"){
                        inFile >> this->node_num >> this->edge_num;
                        cout<<"Node number: "<<node_num<<"\tEdge number: "<<edge_num<<endl;
                        cout << "Some edges:" << endl;
                    }else{
                        getline(inFile, temp_str);
                    }

                    //break;
                }
                else if (line_symbol == "a") {//read graph data
                    inFile >> ID1 >> ID2 >> weight;
//                if(invalid_edges.find(make_pair(ID1,ID2))!=invalid_edges.end()){//if found in invalid set
//                    continue;
//                }
                    temp_max = max(temp_max,weight);
                    temp_min = min(temp_min,weight);

//                    MCEdgesMap[ID1 - 1][ID2-1].w[1] = weight;
                    MCEdgesMap[ID1][ID2].w[1] = weight;
                    if (num_line < num_show) {
                        cout << ID1 << "\t" << ID2 << "\t" << weight << endl;
                        ++num_line;
                    }
                    else if (num_line == num_show) {
                        cout << "..." << endl;
                        ++num_line;
                    }

                }else{
                    getline(inFile, temp_str);
                }
            }
        }



        cout << "Data loaded. ";
        inFile.close();
        tt.stop();
        cout << "The number used for data loading: " << tt.GetRuntime() << " s." << endl;
        cout << "Number of nodes:" << node_num << endl;
        cout << "Number of edges:" << edge_num << endl;
        cout << "The maximal edge weight is " << temp_max << " , while the minimal edge weight is "<< temp_min << endl;

        cout << "--------------------" << endl;
//        num_criteria = 1;

        if(ifBase){
            base_max = temp_max;
            base_min = temp_min;
            max_weights.assign(5,INT32_MIN);
            min_weights.assign(5,INT32_MAX);
            max_weights[0] = base_max;
            min_weights[0] = base_min;
        }else{
            max_weights[1] = temp_max;
            min_weights[1] = temp_min;
        }

    }
    //function for reading MCEdges from disk
    void Graph_pre::Read_MCEdges(const string& filename){
        Timer tt;
        tt.start();
        string lineStr;//data of line
        string line_symbol;
        int num_line = 0;
        int num_show = 6;//number of lines to be printed
        string temp_str;
        int ID1, ID2, weight;
        int num_cri;
        MCEdge temp_edge;
        vector<int> temp_w;
        int temp_id = -1;
        uint temp_degree = 1;

        max_weights.assign(NUM_OF_CRITERIA,INT32_MIN);
        min_weights.assign(NUM_OF_CRITERIA,INT32_MAX);

        //Open file
        ifstream inFile(filename, ios::in);
        if (!inFile) {
            cout << "File opening failed." << endl;
            assert(inFile);
            exit(1);
        }
        cout << "MCEdges Data loading..." << endl;
        while (getline(inFile, lineStr)) {//read each line to lineStr
            if (lineStr != "") {
                istringstream sin(lineStr); //read each string to sin
                sin >> line_symbol;
                if (line_symbol == "p") {//read graph basic information
                    sin >> temp_str;
                    if(temp_str == "sp"){
                        sin >> node_num >> edge_num;
                        MCEdges.assign(node_num,vector<MCEdge>());//allocate space for MCEdges
                    }
                    else if(temp_str == "criteria"){
                        sin >> num_cri;
                        num_criteria = num_cri;
                        cout << "There are " << num_cri <<" criteria in this graph." << endl;
                        cout << "The criteria includes: ";
                        for(int i=0;i<num_cri;++i){
                            sin >> temp_str;
                            cout << temp_str << " " ;
                            if(i<num_cri)
                                mc_criteria.push_back(temp_str);
                        }
                        cout << endl;
                    }
                }
                else if (line_symbol == "a") {//read graph data
                    sin >> ID1 >> ID2;
                    if(temp_id != ID1){
                        if(temp_id != -1){
                            max_degree = max(max_degree,temp_degree);
                            ave_degree += temp_degree;
                        }
                        temp_degree = 1;
                        temp_id = ID1;
                    }else{
                        ++temp_degree;
                    }
                    temp_w.clear();
                    for(int i=0;i<num_cri;++i) {//get edge weights
                        sin >> weight;
                        max_weights[i] = max(max_weights[i],weight);
                        min_weights[i] = min(min_weights[i],weight);
                        if(i<num_cri)
                            temp_w.push_back(weight);
                    }
//                    MCEdges[ID1-1].emplace_back(MCEdge(ID2-1,temp_w));
                    MCEdges[ID1].emplace_back(MCEdge(ID2,temp_w));
                    if (num_line < num_show) {//show some edges
                        cout << ID1 << "\t" << ID2 ;
                        for(auto it=temp_w.begin();it!=temp_w.end();++it){
                            cout << "\t" << *it;
                        }
                        cout << endl;
                        ++num_line;
                    }
                    else if (num_line == num_show) {
                        cout << "..." << endl;
                        ++num_line;
                    }

                }
            }
        }
        ave_degree /= node_num;
        tt.stop();
        cout << "The time used for data loading:" << tt.GetRuntime() << " s." << endl;
        inFile.close();
    }



    //function for reading MCEdges from disk, store as MCEdgesMap. compute the basic statics information too.
    void Graph_pre::Read_MCEdgesMap(const string& filename){//for MCEdgesMap
        Timer tt;
        tt.start();
        string lineStr;//data of line
        string line_symbol;
        int num_line = 0;
        int num_show = 6;//number of lines to be printed
        string temp_str;
        int ID1, ID2, weight;
        int num_cri;
        MCEdge temp_edge;
        vector<int> temp_w;
        int temp_id = -1;
        uint temp_degree = 1;

        max_weights.assign(NUM_OF_CRITERIA,INT32_MIN);
        min_weights.assign(NUM_OF_CRITERIA,INT32_MAX);

        //Open file
        ifstream inFile(filename, ios::in);
        if (!inFile) {
            cout << "File opening failed." << endl;
            assert(inFile);
            exit(1);
        }
        cout << "MCEdges Data loading..." << endl;
        while (getline(inFile, lineStr)) {//read each line to lineStr
            if (lineStr != "") {
                istringstream sin(lineStr); //read each string to sin
                sin >> line_symbol;
                if (line_symbol == "p") {//read graph basic information
                    sin >> temp_str;
                    if(temp_str == "sp"){
                        sin >> node_num >> edge_num;
                        MCEdgesMap.assign(node_num,map<NodeId,MCEdge>());//allocate space for MCEdges
                    }
                    else if(temp_str == "criteria"){
                        sin >> num_cri;
                        num_criteria = num_cri;
                        cout << "There are " << num_cri <<" criteria in this graph." << endl;
                        cout << "The criteria includes: ";
                        for(int i=0;i<num_cri;++i){
                            sin >> temp_str;
                            cout << temp_str << " " ;
                            if(i<num_cri)
                                mc_criteria.push_back(temp_str);
                        }
                        cout << endl;
                    }
                }
                else if (line_symbol == "a") {//read graph data
                    sin >> ID1 >> ID2;
                    if(temp_id != ID1){
                        if(temp_id != -1){
                            max_degree = max(max_degree,temp_degree);
                            ave_degree += temp_degree;
                        }
                        temp_degree = 1;
                        temp_id = ID1;
                    }else{
                        ++temp_degree;
                    }
                    temp_w.clear();
                    for(int i=0;i<num_cri;++i) {//get edge weights
                        sin >> weight;
                        max_weights[i] = max(max_weights[i],weight);
                        min_weights[i] = min(min_weights[i],weight);
                        if(i<num_cri)
                            temp_w.push_back(weight);
                    }
//                    MCEdgesMap[ID1-1][ID2-1]=MCEdge(ID2-1,temp_w);
                    MCEdgesMap[ID1][ID2]=MCEdge(ID2,temp_w);
                    if (num_line < num_show) {//show some edges
                        cout << ID1 << "\t" << ID2 ;
                        for(auto it=temp_w.begin();it!=temp_w.end();++it){
                            cout << "\t" << *it;
                        }
                        cout << endl;
                        ++num_line;
                    }
                    else if (num_line == num_show) {
                        cout << "..." << endl;
                        ++num_line;
                    }

                }
            }
        }
        ave_degree /= node_num;
        tt.stop();
        cout << "The time used for data loading:" << tt.GetRuntime() << " s." << endl;
        inFile.close();
    }

//function for reading MCEdges from disk, store as MCEdgesMap. compute the basic statics information too.
    void Graph_pre::Read_MCEdgesMapS(const string& filename){//for MCEdgesMap
        Timer tt;
        tt.start();
        string lineStr;//data of line
        string line_symbol;
        int num_line = 0;
        int num_show = 6;//number of lines to be printed
        string temp_str;
        int ID1, ID2, weight;
        int num_cri;
        MCEdge temp_edge;
        vector<int> temp_w;
        int temp_id = -1;
        uint temp_degree = 1;
        int i = 0;

        //Open file
        ifstream inFile(filename, ios::in);
        if (!inFile) {
            cout << "File opening failed." << endl;
            assert(inFile);
            exit(1);
        }
        cout << "MCEdgesS Data loading..." << endl;

        inFile >> node_num >> edge_num >> num_cri;
        cout<<"Node number: "<<node_num<<"\tEdge number: "<<edge_num<<endl;
        MCEdgesMap.assign(node_num,map<NodeId,MCEdge>());
        num_criteria = num_cri;
        max_weights.assign(num_criteria,INT32_MIN);
        min_weights.assign(num_criteria,INT32_MAX);
        mc_criteria.assign(num_criteria,"");
        for(i=0;i<num_criteria;++i){
            inFile >> mc_criteria[i];
        }
        while(inFile){
            inFile >> ID1 >> ID2;
            if(temp_id != ID1){
                if(temp_id != -1){
                    max_degree = max(max_degree,temp_degree);
                    ave_degree += temp_degree;
                }
                temp_degree = 1;
                temp_id = ID1;
            }else{
                ++temp_degree;
            }
            temp_w.clear();
            for(i=0;i<num_criteria;++i) {//get edge weights
                inFile >> weight;
                if(weight != INF100){
                    max_weights[i] = max(max_weights[i],weight);
                    min_weights[i] = min(min_weights[i],weight);
                }
                temp_w.push_back(weight);
            }
//            MCEdgesMap[ID1-1][ID2-1]=MCEdge(ID2-1,temp_w);
            MCEdgesMap[ID1][ID2]=MCEdge(ID2,temp_w);
            if (num_line < num_show) {//show some edges
                cout << ID1 << "\t" << ID2 ;
                for(auto it=temp_w.begin();it!=temp_w.end();++it){
                    cout << "\t" << *it;
                }
                cout << endl;
                ++num_line;
            }
            else if (num_line == num_show) {
                cout << "..." << endl;
                ++num_line;
            }

        }

        ave_degree /= node_num;
        tt.stop();
        cout << "The time used for data loading:" << tt.GetRuntime() << " s." << endl;

        inFile.close();
    }

    void Graph_pre::Read_Multiplex(const string& filename){
        Timer tt;
        tt.start();
        int layerID = -1;
        int num_line = 0;
        int num_show = 6;//number of lines to be printed
        int tempID, ID1, ID2, weight;
        vector<EdgeWeight> temp_w;
//        int temp_max = INT32_MIN;
//        int temp_min = INT32_MAX;
        vector<int> edgenum;

        //Open file
        ifstream inFile(filename, ios::in);
        if (!inFile) {
            cout << "File opening failed." << endl;
            assert(inFile);
            exit(1);
        }
        cout << dataset<<" multiplex graph edges Data loading..." << endl;
        inFile >> node_num >> num_criteria;
        cout << "Node number: " << node_num << "\tNumber of criteria: " << num_criteria << endl;
        max_weights.assign(num_criteria,INT32_MIN);
        min_weights.assign(num_criteria,INT32_MAX);
        MCEdgesMap.assign(node_num,map<NodeId,MCEdge>());
        edgenum.assign(num_criteria,0);
        mc_criteria.resize(num_criteria);
        for(int i=0;i<num_criteria;++i){
            inFile >> mc_criteria[i];
        }
        edge_num = 0;
        while (inFile) {//read each line to lineStr
            inFile >> tempID >> ID1 >> ID2 >> weight;
            if(ID1 == ID2 || weight<=0){
//                cout << "Wrong edge e("<<ID1<<", "<<ID2<<") = "<<weight << endl;
                continue;
            }
            if (tempID != layerID+1){
                if(layerID >= 0){
                    cout << "Edge number of layer " << layerID+1 << " is: " << edgenum[layerID]<<endl;
                    cout << "The maximal edge weight is " << max_weights[layerID] << " , while the minimal edge weight is "<< min_weights[layerID] << endl;
                }

                cout << "New layer " << tempID << "." << endl;
                layerID = tempID-1;
            }

            max_weights[layerID] = max(max_weights[layerID],weight);
            min_weights[layerID] = min(min_weights[layerID],weight);
//            edge_num++;

            if(MCEdgesMap[ID1].find(ID2) == MCEdgesMap[ID1].end()){//if not found
                temp_w.assign(num_criteria,INF100);
                temp_w[layerID] = weight;
                MCEdgesMap[ID1].insert({ID2,MCEdge(ID2, temp_w)});//-1
                MCEdgesMap[ID2].insert({ID1,MCEdge(ID1, temp_w)});//-1
                edgenum[layerID]++;
                edge_num+=2;
            }else{//if found
                MCEdgesMap[ID1][ID2].w[layerID] = weight;
                MCEdgesMap[ID2][ID1].w[layerID] = weight;
                edgenum[layerID]++;
            }

            if (num_line < num_show) {
                cout << ID1 << "\t" << ID2 << "\t" << weight << endl;
                ++num_line;
            }
            else if (num_line == num_show) {
                cout << "..." << endl;
                ++num_line;
            }

        }

        cout << "Edge number of layer " << layerID+1 << " is: " << edgenum[layerID]<<endl;
        cout << "The maximal edge weight is " << max_weights[layerID] << " , while the minimal edge weight is "<< min_weights[layerID] << endl;

        cout << "Data loaded. ";
        inFile.close();
        tt.stop();
        cout << "The number used for data loading: " << tt.GetRuntime() << " s." << endl;
        uint temp_deg;
        for(auto it=MCEdgesMap.begin();it!=MCEdgesMap.end();it++){
            temp_deg = it->size();
            ave_degree += temp_deg;
            max_degree = max(max_degree,temp_deg);
            min_degree = min(min_degree,temp_deg);
        }
        ave_degree/=node_num;

        cout << "Number of nodes: " << node_num << endl;
        cout << "Number of edges: " << edge_num << endl;
        int sum_edge=0;
        for(int i=0;i<edgenum.size();i++)
            sum_edge+=edgenum[i];

        cout << "Summation of the edges of all layer: " << sum_edge << endl;
        cout << "Maximum degree: " << max_degree <<"\t Minimum degree: " << min_degree <<"\t Average degree: " <<ave_degree<<endl;
        cout << "--------------------" << endl;

    }

    //function for writing the multiplex network into disk
    void Graph_pre::Write_Multiplex(const string& w_file){

        /// Write
        ofstream outFile(w_file, ios::out);
        if (!outFile) {
            cout << "File opening failed." << endl;
            assert(outFile);
        }
        cout << "Writing simplified MCEdges data..."<<endl;
        outFile<<this->node_num<<" "<<this->edge_num<<endl;
        outFile<<num_criteria;

        for (int i = 0; i < num_criteria; ++i) {
            if(i == num_criteria-1)
                outFile<<" "<<mc_criteria[i]<<endl;
            else
                outFile<<" "<<mc_criteria[i];
        }
        int ID1 = 0;

        for(auto it=MCEdgesMap.begin();it!=MCEdgesMap.end();++it){
            for(auto it2=it->begin();it2!=it->end();it2++){
                outFile << ID1<<" "<< it2->first;//+1
                for(auto it3=it2->second.w.begin();it3!=it2->second.w.end();it3++){
                    outFile <<" "<< *it3;
                }
                outFile << endl;
            }
            ++ID1;
        }
        outFile.close();
        assert(ID1 == node_num);
        cout << "Write done." << endl;
    }


    //function of writing the largest component to disk
    void Graph_pre::Component_WriteMultiplex(set<NodeId>& w_set,uint new_edge_num){
        string filename = string(DataPath)  + dataset + "/" + dataset + "_MCEdgesS_Original.txt";
        uint temp_degree;
        unsigned long long i=0;
        //get the map from old vertex id to the new one
//        for (auto it=w_set.begin(); it != w_set.end(); ++it) {
//            nodes_map[*it] = i;
//            ++i;
//        }
        //update node_num and edge_num
        this->node_num = w_set.size();
        this->edge_num = new_edge_num;

        vector<EdgeWeight> temp_w;
        //write edges to disk
        cout << "Writing Edges data..." << endl;

        ofstream outFile(filename, ios::out);
        if (!outFile) {
            cout << "File opening failed." << endl;
            assert(outFile);
            exit(1);
        }
        cout << "Writing simplified MCEdges data..."<<endl;
        outFile<<this->node_num<<" "<<new_edge_num<<endl;
        outFile<<num_criteria;

        for (int i = 0; i < num_criteria; ++i) {
            if(i == num_criteria-1)
                outFile<<" "<<mc_criteria[i]<<endl;
            else
                outFile<<" "<<mc_criteria[i];
        }

        i=0; int ID1=0;
        for (auto it=w_set.begin(); it != w_set.end(); ++it) {//get ID1
            for(auto it2=MCEdgesMap[*it].begin();it2!=MCEdgesMap[*it].end();it2++){//get ID2
                outFile << *it<<" "<< it2->first;//+1
                for(auto it3=it2->second.w.begin();it3!=it2->second.w.end();it3++){
                    outFile <<" "<< *it3;
                }
                outFile << endl;
                i += 1;
            }

            ++ID1;
            temp_degree = MCEdgesMap[*it].size();
//                Nodes[nodes_map[*it]].degree = temp_degree;
            max_degree = max(max_degree,temp_degree);
            ave_degree += temp_degree;
        }
        outFile.close();
        assert(i==new_edge_num);
        if(i!=new_edge_num){
            cout<<"!!!!!!! The edge number is different! "<<new_edge_num<<" "<<i<<endl;
        }
        cout << "Write done." << endl;
        cout<<"Edges size of new graph: "<<i<<endl;

        ave_degree /= node_num;
        cout << "The maximal vertex degree is " << max_degree << " , the average degree is " << ave_degree << endl;

    }

    //function for reading partitions
//    void Graph_pre::Partition_Read(){
//        Timer tt;
//        tt.start();
//        string filename = string(DataPath) + this->dataset + "/" + this->dataset + "_Partition_node_clusters.txt";
//        string file_w_map = string(DataPath) + this->dataset + "/" + this->dataset + "_Partition_map.txt";
//        string file_w_cluster = string(DataPath) + this->dataset + "/" + this->dataset + "_Partition_info.txt";
//        string lineStr;//data of line
//        string line_symbol;
//        string temp_str;
//        int temp_int;
//        int cluster_id;
//        int new_id = 0;
//
//        /*** Read Partitions ***/
//        //Open file
//        ifstream inFile(filename, ios::in);
//        if (!inFile) {
//            cout << "File opening failed." << endl;
//            assert(inFile);
//            exit(1);
//        }
//        cout << "Partition Data loading..." << endl;
//        while (getline(inFile, lineStr)) {//read each line to lineStr
//            if (lineStr != "") {
//                istringstream sin(lineStr); //read each string to sin
//                sin >> line_symbol;
//                if (line_symbol == "p") {//read graph basic information
//                    sin >> temp_str;
//                    if(temp_str == "all"){
//                        sin >> temp_int;
//                        partitions_info.assign(temp_int,pair<int,pair<int,int>>(0, make_pair(INT32_MAX,0)));//allocate space for MCEdges
//                    }
//                    else if(temp_str == "partition"){
//                        sin >> cluster_id;
//                    }
//                }
//                else if (line_symbol == "v") {//read graph data
//                    sin >> temp_int;
//                    old_to_new[temp_int-1] = new_id;
//                    if(new_id < partitions_info[cluster_id].second.first){
//                        partitions_info[cluster_id].second.first = new_id;
//                    }
//                    if(new_id > partitions_info[cluster_id].second.second){
//                        partitions_info[cluster_id].second.second = new_id;
//                    }
//                    ++new_id;
//                }
//            }
//        }
//        tt.stop();
//        cout << "The time used for data loading:" << tt.GetRuntime() << " s." << endl;
//        inFile.close();
//        cout << "Vertices number: " << new_id << " "<<old_to_new.size()<< endl;
//        /*** Write node map file ***/
//        ofstream outFile_m(file_w_map, ios::out);
//        if (!outFile_m) {
//            cout << "Write File opening failed." << endl;
//            assert(outFile_m);
//            exit(1);
//        }
//
//        new_to_old.assign(old_to_new.size(),-1);
//        outFile_m<<"new_id old_id"<<endl;
//        for(auto it=old_to_new.begin();it!=old_to_new.end();++it){
//            new_to_old[it->second] = it->first;
//        }
//        for(int i=0;i<new_to_old.size();++i){
//            outFile_m<<i+1<<" "<<new_to_old[i]+1<<endl;//The first element is old id, the second one is new id
//        }
//        outFile_m.close();
//        /*** Write partition information file ***/
//        ofstream outFile_c(file_w_cluster, ios::out);
//        if (!outFile_c) {
//            cout << "Write File opening failed." << endl;
//            assert(outFile_c);
//            exit(1);
//        }
//        outFile_c << partitions_info.size()<<endl;
//        for(int i=0;i<partitions_info.size();++i){
//            partitions_info[i].first = partitions_info[i].second.second - partitions_info[i].second.first + 1;
//            outFile_c<< i <<" "<<partitions_info[i].first<<" "<<partitions_info[i].second.first<<" "<<partitions_info[i].second.second<<endl;
//        }
//        outFile_c.close();
//    }
////function for writing partitions to files
//    void Graph_pre::Partition_Write(){
//        map<pair<int,int>,vector<int>>::iterator iter = edge_list.begin();
//        for(int i=0;i<partitions_info.size();++i){
//            string file_w_p = DataPath + this->dataset + "/Partitions/" + this->dataset + "_Partition_"+ to_string(i);
//            ofstream outFile(file_w_p, ios::out);
//            if (!outFile) {
//                cout << "Write File opening failed." << endl;
//                assert(outFile);
//                exit(1);
//            }
//            outFile << to_string(i)<<" "<<partitions_info[i].first<<" "<<partitions_info[i].second.first+1<<" "<<partitions_info[i].second.second+1<<endl;
//            for (int j = partitions_info[i].second.first; j <= partitions_info[i].second.second; ++j) {
//                outFile << j + 1 << " " << Nodes[j].degree;
//                for(int k=0;k<Nodes[j].degree;k++){
//                    outFile << " " << iter->first.second + 1;
//                    for (auto it2 = iter->second.begin(); it2 != iter->second.end(); it2++) {
//                        outFile << " " << *it2;
//                    }
//                    ++iter;
//                }
//                outFile << endl;
//            }
//            outFile.close();
//        }
//    }
//
//    void Graph_pre::Read_MC(bool ifCoordinate){
//        string file_r_mc = DataPath + this->dataset + "/" + this->dataset + "_MCEdges.txt";
//        string file_r_co = DataPath + this->dataset + "/" + this->dataset + "_Coordinates.txt";
//        Read_MCEdges(file_r_mc);
//
//    }
}

#endif //MCSPS_PREPROCESS_HPP
