/*
 * Filename:    emgraph.hpp
 * Description: functions for external-memory graph data processing
 * Created:     04 March 2021
 * Authors:     Xinjie ZHOU
 */

#ifndef MCSPS_EMGRAPH_HPP
#define MCSPS_EMGRAPH_HPP

#include "emgraph.h"

namespace gbxxl{
    //one-off edges reading for MC-graph
    /*void EMMCGraph::MCReadGraph_MCEdges(string filename) {
        Timer tt;
        tt.start();
        string line_symbol; //line symbol in the beginning of each line
        string temp_str;    //temp variable of string
        int num_line = 0;   //recording the number of lines
        int num_show = 6;   //line number of graph data to be printed
        int num_cri = 5;
        int ID1, ID2, weight;
        MCEdge temp_edge;
        vector<int> temp_w;

        //Open file
        ifstream inFile(filename, ios::in);
        if (!inFile) {
            cout << "File opening failed." << endl;
            exit(1);
        }
        cout << "MCEdges Data loading..." << endl;
        while (inFile) {
            inFile >> line_symbol;
            if (line_symbol == "p") {//read graph basic information
                inFile >> temp_str;
                if (temp_str == "sp") {
                    inFile >> node_num >> edge_num;
                    MCEdges.assign(node_num, vector<MCEdge>());//allocate space for MCEdges
                } else if (temp_str == "criteria") {
                    inFile >> num_cri;
                    cout << "There are " << num_cri << " criteria in this graph." << endl;
                    cout << "The criteria includes: ";
                    for (int i = 0; i < num_cri; ++i) {
                        inFile >> temp_str;
                        cout << temp_str << " ";
                        if (i < num_of_cri)
                            mc_criteria.emplace_back(temp_str);
                    }
                    cout << endl;
                }
            } else if (line_symbol == "a") {//read graph data
                inFile >> ID1 >> ID2;
                temp_w.clear();
                for (int i = 0; i < num_cri; ++i) {//get edge weights
                    inFile >> weight;
                    if (i < num_cri)
                        temp_w.emplace_back(weight);
                }
                MCEdges[ID1 - 1].emplace_back(MCEdge(ID2 - 1, temp_w));
                if (num_line < num_show) {//show some edges
                    cout << ID1 << "\t" << ID2;
                    for (auto it = temp_w.begin(); it != temp_w.end(); ++it) {
                        cout << "\t" << *it;
                    }
                    cout << endl;
                } else if (num_line == num_show) {
                    cout << "..." << endl;
                }
                ++num_line;
            } else {
                getline(inFile, temp_str);
            }
        }

        cout << "Data loaded. ";
        inFile.close();
        tt.stop();
        cout << "The number used for data loading: " << tt.GetRuntime() << " s." << endl;
        cout << "Number of nodes:" << node_num << endl;
        cout << "Number of edges:" << edge_num << endl;
        cout << "--------------------" << endl;
    }*/
    //function for reading partition information
    void EMMCGraph::ClusterInfoLoad(bool ifMap) {//load partition information, ifMap is used to indicate whether to used node-to-cluster mapping
        /*** Load partition information ***/
        string filename = string(DataPath) + dataset + "/" + partition_type+"/"+dataset + "_Partitions_" + to_string(PartitionSize) + "_Info.txt";
//        string filename = string(DataPath) + dataset + "/" + dataset + "_Partitions_Info.txt";
        int p_id;   //partition id
        int ID1,ID2;

        //Open file
        ifstream inFile(filename, ios::in);
        if (!inFile) {
            cout << "File opening failed." << endl;
            assert(inFile);
            exit(1);
        }
        inFile >> node_num >> edge_num;
        mc_criteria.assign(NUM_OF_CRITERIA,string());
        for(int i=0;i<NUM_OF_CRITERIA;++i){
            inFile >> mc_criteria[i];
        }
        inFile >> partition_number;
        cout << "Partition number: " << partition_number <<endl;
        partitions_info.assign(partition_number+1,0);
        cluster_to_IO.assign(partition_number,0);

        if(ifMap){//if it needs to map vertex id to partition id
            node_to_cluster.assign(node_num,0);
            for(int i=0;i<partition_number;++i){
                inFile >> p_id >> ID1 >> ID2;
                partitions_info[p_id] = ID1;
                if(i==partition_number-1)
                    partitions_info[p_id+1] = ID2+1;
                for(int j=ID1;j<=ID2;++j){
                    node_to_cluster[j] = p_id;
                }
            }
        }else{
            for(int i=0;i<partition_number;++i){
                inFile >> p_id >> ID1 >> ID2;
                partitions_info[p_id] = ID1;
            }
        }
        inFile.close();
    }
    //Initiate graph for each round of IO-efficient algorithms
    void EMMCGraph::CommonInitiation_IO()
    {
        //initialise metrics
        io_num = 0;//update in reading function
        io_time = 0;
        p_num = 0;
        visited_number = 0;//update in functions share_and_slack and slack
        cluster_to_IO.assign(partition_number,0);
    }
    //Disk-based MCSP Evaluation
    void EMMCGraph::MC_Evaluate_IO(const string& qtype)
    {
        /*--variables about time record--*/
        //std::chrono::high_resolution_clock::time_point t_s, t_e;//variables for time record
        double time_EMDijk = 0;
        double time_DijkstraIO = 0;
        double time_BiDijkstraIO = 0;
        double time_OneHopNoIO = 0;
        double time_MultiHopsNoIO = 0;
        double time_BiMultiHopsNoIO = 0;
        double time_OneHop = 0;
        double time_MultiHops = 0;
        double time_BiMultiHops = 0;

        assert(BlockPerPage*Block_SZ <= PartitionSize*1024);//page size should be smaller than partition size
        /*** Read graph data ***/
        ClusterInfoLoad(true);//load partition information, using node-to-cluster mapping
        ///for one-off read testing
//        string gr_Edge = string(DataPath)  + dataset + "/" + dataset + "_MCEdges.txt";
//        ReadGraphToExternalVector(gr_Edge);
//        MCReadGraph_MCEdges(gr_Edge);
        cout << "Node number: " << node_num << "\t\tEdge number: " << edge_num << endl;

        /*** OnePass algorithms ***/
        /// execution of MC_OneHop without IO optimization
        cout << "----- Below are results of MC_OneHop_NoIO algorithm -----" << endl;
        time_OneHopNoIO = MC_OnePass(qtype,OneHopNoIO);
        cout << "-----------------------------------------------\n" << endl;
        /// execution of MC_MultiHop without IO optimization
        cout << "----- Below are results of MC_MultiHop_NoIO algorithm -----" << endl;
        time_MultiHopsNoIO = MC_OnePass(qtype,MultiHopsNoIO);
        cout << "-----------------------------------------------\n" << endl;
        /// execution of MC_BiMultiHop without IO optimization
        cout << "----- Below are results of MC_BiMultiHop_NoIO algorithm -----" << endl;
        time_BiMultiHopsNoIO = MC_OnePass(qtype,BiMultiHopsNoIO);
        cout << "-----------------------------------------------\n" << endl;
        /// execution of MC_OneHop
        cout << "----- Below are results of MC_OneHop algorithm -----" << endl;
        time_OneHop = MC_OnePass(qtype,OneHop);
        cout << "-----------------------------------------------\n" << endl;
        /// execution of MC_MultiHop
        cout << "----- Below are results of MC_MultiHop algorithm -----" << endl;
        time_MultiHops = MC_OnePass(qtype,MultiHops);
        cout << "-----------------------------------------------\n" << endl;
        /// execution of MC_BiMultiHop
        cout << "----- Below are results of MC_BiMultiHop algorithm -----" << endl;
        time_BiMultiHops = MC_OnePass(qtype,BiMultiHops);
        cout << "-----------------------------------------------\n" << endl;

        cout << "The total run time of MC_OneHop_NoIO algorithm is: "<< time_OneHopNoIO /1000<<" s."<<endl;
        cout << "The total run time of MC_MultiHop_NoIO algorithm is: "<< time_MultiHopsNoIO /1000<<" s."<<endl;
        cout << "The total run time of MC_BiMultiHop_NoIO algorithm is: " << time_BiMultiHopsNoIO / 1000 << " s." << endl;
        cout << "The total run time of MC_OneHop algorithm is: "<< time_OneHop /1000<<" s."<<endl;
        cout << "The total run time of MC_MultiHop algorithm is: "<< time_MultiHops /1000<<" s."<<endl;
        cout << "The total run time of MC_BiMultiHop algorithm is: " << time_BiMultiHops / 1000 << " s." << endl;
        node_to_cluster.clear();
    }

    //function of reading graph data to HotPool4
    template <class T>
    void EMMCGraph::GraphReadCluster4(HotPool4<T> & mcHotPool,int target_p_id){
//        cout << "target_p_id: " << target_p_id << endl;
        int evict_p_id = mcHotPool.LRU_IMCluster.put(target_p_id,1);
        //// evict old partition for new one
//        cout<<mcHotPool.LRU_IMCluster.ls_.size()<<endl;
        if(evict_p_id != -1) {//if it is necessary to evict old partition
            int temp_id = -1;
            for(int i=partitions_info[evict_p_id];i<partitions_info[evict_p_id+1];++i){//store all elements to stxxl vector
                if(i != temp_id){
                    get<1>(EMEdgesIndex[i]) = mcHotPool.MCEdges_EM.size();
                    if(temp_id != -1){
                        get<2>(EMEdgesIndex[temp_id]) = mcHotPool.MCEdges_EM.size();
                    }
                    temp_id = i;
                }
                assert(mcHotPool.HotPools[evict_p_id].find(i)!=mcHotPool.HotPools[evict_p_id].end());
//                auto range=mcHotPool.HotPools[evict_p_id].equal_range(i);
//                for(auto it=range.first;it!=range.second;++it){
//                    mcHotPool.MCEdges_EM.push_back(it->second);
//                }
                for(auto it=mcHotPool.HotPools[evict_p_id][i].begin();it!=mcHotPool.HotPools[evict_p_id][i].end();++it){
                    mcHotPool.MCEdges_EM.push_back(*it);
                }
            }
            get<2>(EMEdgesIndex[temp_id]) = mcHotPool.MCEdges_EM.size();

            mcHotPool.clusterStatus[evict_p_id] = 'S';
            mcHotPool.HotPools[evict_p_id].clear();
        }
        //// read partition from original disk file to internal vector
        //check status of storage
        char filePath[200];
        int u, u_deg;
        int v, weight;

        strcpy(filePath, DataPath);
        strcat(filePath, dataset.c_str());
        strcat(filePath, "/");
        strcat(filePath, partition_type.c_str());
        strcat(filePath, "/Partitions/");
        strcat(filePath, dataset.c_str());
        strcat(filePath, "_Partition_");
        strcat(filePath, to_string(PartitionSize).c_str());
        strcat(filePath, "_");
        strcat(filePath, to_string(target_p_id).c_str());
        strcat(filePath, ".bin");
        Timer tt;
        tt.start();
        FILE* file = fopen(filePath, "rb");
        assert(file != NULL);
        ReadBuffer read_io(file);
        bool flag_sizeUpdate = false;
        int temp_id = 0;
        MCEdgeT temp_edge;
        vector<MCEdgeT> temp_vedge;
        //read partition
        while(!read_io.isend) {//if not read out
            read_io.read(&u);//get node id
            temp_edge.ID1 = u-1;
            read_io.read(&u_deg);//get node degree
            //vertex store in indexNode, edge store in e
            for (int i = 0; i < u_deg; ++i) {
                read_io.read(&v);//read end node

                temp_edge.ID2 = v-1;
                for (int j = 0; j < NUM_OF_CRITERIA; ++j) {
                    read_io.read(&weight);
                    if(j<num_of_cri)
                        temp_edge.putW(j,weight);
//                    if(j == sc_i){
//                        mcHotPool.HotPools[target_p_id].emplace(u-1,Edge(v-1,weight));
//                    }
                }
                temp_vedge.emplace_back(temp_edge);
//                mcHotPool.HotPools[target_p_id].insert({u-1,temp_edge});
            }
            mcHotPool.HotPools[target_p_id].insert({u-1,temp_vedge});
            temp_vedge.clear();
        }

        io_num += read_io.io_number;//update #IO
        fclose(file);
        tt.stop();
        io_time += tt.GetRuntime() * 1000;
        cluster_to_IO[target_p_id] = read_io.io_number;

        mcHotPool.clusterStatus[target_p_id] = 'I';
    }
    //Entry for one-pass algorithms
    double EMMCGraph::MC_OnePass(const string& qtype,int algo_choice){
        string LongDis,MediumDis,ShortDis,RandomDis;
        LongDis = string(DataPath) + dataset + "/" + partition_type+"/"+dataset + "_OD_LongDis_" + to_string(PartitionSize) + ".txt";
        MediumDis = string(DataPath) + dataset + "/" + partition_type+"/"+dataset + "_OD_MediumDis_"+ to_string(PartitionSize) +".txt";
        ShortDis = string(DataPath) + dataset + "/" + partition_type+"/"+dataset + "_OD_ShortDis_"+ to_string(PartitionSize) +".txt";
        RandomDis = string(DataPath) +  dataset + "/" + partition_type+"/"+dataset + "_OD_Random_"+ to_string(PartitionSize) +".txt";

        double total_time = 0;

        ///Shortest path querying
        if(qtype == "S"){
            QueryType = "S";
            total_time = MC_OnePass_one(ShortDis, "SDis",algo_choice);//Efficiency evaluation on short distance
        }else if(qtype == "M"){
            QueryType = "M";
            total_time = MC_OnePass_one(MediumDis, "MDis",algo_choice);//Efficiency evaluation on medium distance
        }else if(qtype == "L"){
            QueryType = "L";
            total_time = MC_OnePass_one(LongDis, "LDis",algo_choice);//Efficiency evaluation on long distance
        }else if(qtype == "R"){
            QueryType = "R";
            total_time = MC_OnePass_one(RandomDis, "RDis",algo_choice);//Efficiency evaluation on long distance
        }else if (qtype == "all") {
            QueryType = "S";
            total_time += MC_OnePass_one(ShortDis, "SDis", algo_choice);//Efficiency evaluation on short distance
            QueryType = "M";
            total_time += MC_OnePass_one(MediumDis, "MDis", algo_choice);//Efficiency evaluation on medium distance
            QueryType = "L";
            total_time += MC_OnePass_one(LongDis, "LDis", algo_choice);//Efficiency evaluation on long distance
        }
        PeakMemory();
        return total_time;
    }
    //Function for calling different types of one-pass algorithms
    double EMMCGraph::MC_OnePass_one(const string& filename, const string& qtype, int algo_choice){
        string fname;
        int num, ID1, ID2;
        vector<pair<int, int>> ODpair;
        uint ave_visited = 0;
        double ave_io_num = 0;
        uint ave_share = 0;
        double ave_time = 0;
        double ave_io_time = 0;
        double ave_p_num = 0;
        vector<Distance> ave_cost(num_of_cri,0);
        query_time = 0;

        //Open file
        ifstream inFile(filename, ios::in);
        if (!inFile) {
            cout << "File opening failed." << endl;
            exit(1);
        }
        inFile >> num;
        for (int i = 0; i < run_times; ++i) {
            inFile >> ID1 >> ID2;
            ODpair.emplace_back(make_pair(ID1, ID2));
        }
        inFile.close();
        uint basic_mem;
        
        switch (algo_choice){
            case OneHopNoIO:{
                basic_mem = node_num*(num_of_cri*5 + 11)/(1024*1024) + PQMemory*num_of_cri;
                cout << "*** The basic memory consumption of OneHopNoIO is: " << basic_mem << " MB. The graph data memory is "<< MemGraph_NoIO<<" MB. ***"<<  endl; break;
            }
            case MultiHopsNoIO:{
                basic_mem = node_num*(num_of_cri*5 + 11)/(1024*1024) + PQMemory*num_of_cri;
                cout << "*** The basic memory consumption of MultiHopsNoIO is: " << basic_mem << " MB. The graph data memory is "<< MemGraph_NoIO<<" MB. ***"<<  endl; break;
            }
            case BiMultiHopsNoIO:{
                cout << "PageNumber4_NoIO_Bi: " << PageNumber4_NoIO_Bi << endl;
                basic_mem = node_num*(num_of_cri*10 + 11)/(1024*1024) + 2*PQMemory*num_of_cri;
                cout << "*** The basic memory consumption of BiMultiHopsNoIO is: " << basic_mem << " MB. The graph data memory is "<< MemGraph_NoIO_Bi<<" MB. ***"<<  endl; break;
            }
            case OneHop:{
                basic_mem = node_num*(num_of_cri*5 + 12)/(1024*1024) + PQMemory*num_of_cri;
                cout << "*** The basic memory consumption of OneHop is: " << basic_mem << " MB. The graph data memory is "<< MemGraph_IO<<" MB. ***"<< endl;
                cout << "Alpha: " << alpha << "\tMu: " << double(MuForEM)/10 << endl;
                break;
            }
            case MultiHops:{
                basic_mem = node_num*(num_of_cri*5 + 12)/(1024*1024) + PQMemory*num_of_cri;
                cout << "*** The basic memory consumption of MultiHops is: " << basic_mem << " MB. The graph data memory is "<< MemGraph_IO<<" MB. ***"<<  endl;
                cout << "Alpha: " << alpha_multi <<"\tMu: " << double(MuForEM) / 10 << endl;
                break;
            }
            case BiMultiHops:{
                basic_mem = node_num*(num_of_cri*10 + 12)/(1024*1024) + 2*PQMemory*num_of_cri;
                cout << "*** The basic memory consumption of BiMultiHops is: " << basic_mem << " MB. The graph data memory is "<< MemGraph_IO_Bi<<" MB. ***"<<  endl;
                cout << "Alpha: " << alpha_bi <<"\tMu: " << double(MuForEM_Bi) / 10 << endl;
                break;
            }
            default:
                cout<<"\nwrong!"<<endl; break;
        }
        cout << "Query type: " << QueryType << " \tRun times: " << run_times << endl;
        EM_MC_PQueue.resize(num_of_cri);
        gbxxl::PriorityQueue pqueue0(gbxxl::PQ_Pool); gbxxl::PriorityQueue pqueue1(gbxxl::PQ_Pool);
        gbxxl::PriorityQueue pqueue2(gbxxl::PQ_Pool);
        gbxxl::PriorityQueue pqueue3(gbxxl::PQ_Pool); gbxxl::PriorityQueue pqueue4(gbxxl::PQ_Pool);
        if(num_of_cri>=2){
            EM_MC_PQueue[0] = &pqueue0; EM_MC_PQueue[1] = &pqueue1;
            if(num_of_cri>=3){
                EM_MC_PQueue[2] = &pqueue2;
                if(num_of_cri>=4){
                    EM_MC_PQueue[3] = &pqueue3;
                    if(num_of_cri==5){
                        EM_MC_PQueue[4] = &pqueue4;
                    }
                }
            }
        }

        if(algo_choice == BiMultiHops || algo_choice == BiMultiHopsNoIO){
            EM_MC_PQueue_r.resize(num_of_cri);
            gbxxl::PriorityQueue2 pqueue0_r(gbxxl::PQ_Pool2); gbxxl::PriorityQueue2 pqueue1_r(gbxxl::PQ_Pool2);
            gbxxl::PriorityQueue2 pqueue2_r(gbxxl::PQ_Pool2); gbxxl::PriorityQueue2 pqueue3_r(gbxxl::PQ_Pool2);
            gbxxl::PriorityQueue2 pqueue4_r(gbxxl::PQ_Pool2);
            if(num_of_cri>=2){
                EM_MC_PQueue_r[0] = &pqueue0_r; EM_MC_PQueue_r[1] = &pqueue1_r;
                if(num_of_cri>=3){
                    EM_MC_PQueue_r[2] = &pqueue2_r;
                    if(num_of_cri>=4){
                        EM_MC_PQueue_r[3] = &pqueue3_r;
                        if(num_of_cri==5){
                            EM_MC_PQueue_r[4] = &pqueue4_r;
                        }
                    }
                }
            }
            for (int i = 0; i < run_times; ++i) {//0
                ID1 = ODpair[i].first;
                ID2 = ODpair[i].second;
//                cout << "Query "<<i<<" : "<<ID1-1 << " " << ID2-1 << endl;
                CommonInitiation_IO();//common initiation of each round
                OnePassInitiation_IO(algo_choice);

                if(algo_choice == BiMultiHops){
                    EM_MC_BiMultiHop(ID1-1,ID2-1);
                }else if(algo_choice == BiMultiHopsNoIO){
//                    EM_MC_BiMultiHop_NoIO(ID1-1,ID2-1);
                    EM_MC_BiMultiHop_NoIO_new(ID1-1,ID2-1);
                }

                //post-processing
//                for(int a=0;a<num_of_cri;++a){
//                    MCDij_getPath(ID1-1, ID2-1, a);
//                }
                //clear priority queue
                while(!pqueue0.empty()) pqueue0.pop();
                while(!pqueue1.empty()) pqueue1.pop();
                while(!pqueue2.empty()) pqueue2.pop();
                while(!pqueue3.empty()) pqueue3.pop();
                while(!pqueue4.empty()) pqueue4.pop();
                OnePassClear_IO(algo_choice);

                while (!pqueue0_r.empty()) pqueue0_r.pop();
                while (!pqueue1_r.empty()) pqueue1_r.pop();
                while (!pqueue2_r.empty()) pqueue2_r.pop();
                while (!pqueue3_r.empty()) pqueue3_r.pop();
                while (!pqueue4_r.empty()) pqueue4_r.pop();

                for(int j=0;j<num_of_cri;++j){
                    if(mc_min_cost[j]<INF100)
                        ave_cost[j] += mc_min_cost[j];
                }
                ave_visited += visited_number;
                ave_share += share_number;
                ave_io_time += io_time;
                ave_io_num += io_num;
                ave_p_num += p_num;
            }
        }else{
            for (int i = 0; i < run_times; ++i) {
                ID1 = ODpair[i].first;
                ID2 = ODpair[i].second;
//                cout << "Query "<<i<<" : "<<ID1-1 << " " << ID2-1 << endl;

                CommonInitiation_IO();//common initiation of each round
                OnePassInitiation_IO(algo_choice);
                switch (algo_choice){
                    case OneHopNoIO:{
//                        EM_MC_OneHop_NoIO(ID1-1, ID2-1); break;
                        EM_MC_OneHop_NoIO_new(ID1-1, ID2-1); break;
                    }
                    case MultiHopsNoIO:{
//                        EM_MC_MultiHop_NoIO(ID1-1, ID2-1); break;
                        EM_MC_MultiHop_NoIO_new(ID1-1, ID2-1); break;
                    }
                    case OneHop:{
                        EM_MC_OneHop(ID1-1, ID2-1); break;
                    }
                    case MultiHops:{
                        EM_MC_MultiHop(ID1-1, ID2-1); break;
                    }
                    default:
                        cout<<"\nwrong!"<<endl; break;
                }
                //post-processing
//                for(int a=0;a<num_of_cri;++a){
//                    MCDij_getPath(ID1-1, ID2-1, a);
//                }
                //clear priority queue
                while(!pqueue0.empty()) pqueue0.pop();
                while(!pqueue1.empty()) pqueue1.pop();
                while(!pqueue2.empty()) pqueue2.pop();
                while(!pqueue3.empty()) pqueue3.pop();
                while(!pqueue4.empty()) pqueue4.pop();
                OnePassClear_IO(algo_choice);

                for(int j=0;j<num_of_cri;++j){
                    if(mc_min_cost[j]<INF100)
                        ave_cost[j] += mc_min_cost[j];
                }
                ave_visited += visited_number;
                ave_share += share_number;
                ave_io_time += io_time;
                ave_io_num += io_num;
                ave_p_num += p_num;
            }
        }

        switch (algo_choice){
            case OneHopNoIO:{
                cout << "Average performance of MC_OneHop_NoIO algorithm (on " << run_times << " " << QueryType << " OD pairs) is shown below " << endl;
                break;
            }
            case MultiHopsNoIO:{
                cout << "Average performance of MC_MultiHop_NoIO algorithm (on " << run_times << " " << QueryType << " OD pairs) is shown below " << endl;
                break;
            }
            case BiMultiHopsNoIO:{
                cout << "Average performance of MC_BiMultiHop_NoIO algorithm (on " << run_times << " " << QueryType << " OD pairs) is shown below " << endl;
                break;
            }
            case OneHop:{
                cout << "Average performance of MC_OneHop algorithm (on " << run_times << " " << QueryType << " OD pairs) is shown below " << endl;
                break;
            }
            case MultiHops:{
                cout << "Average performance of MC_MultiHop algorithm (on " << run_times << " " << QueryType << " OD pairs) is shown below " << endl;
                break;
            }
            case BiMultiHops:{
                cout << "Average performance of MC_BiMultiHop algorithm (on " << run_times << " " << QueryType << " OD pairs) is shown below " << endl;
                break;
            }
            default:
                cout<<"\nThe algorithm choice is wrong!"<<endl;
                break;
        }
        ave_time = query_time / run_times;
        ave_visited /= run_times;
        ave_share /= run_times;
        ave_io_time /= run_times;
        ave_io_num /= run_times;
        ave_p_num /= run_times;
        for(int i=0;i<num_of_cri;++i) {
            ave_cost[i] /= run_times;
            cout << "Average distance of criteria " << mc_criteria[i] <<": " << ave_cost[i] << endl;
        }
        cout << "! Average query time: " << ave_time/1000 << " s." << endl;
        cout << "! Average io time: " << ave_io_time/1000 << " s." << endl;
        cout << "! Average number of IO: " << ave_io_num << endl;
        cout << "! Average number of nodes visited: " << ave_visited << endl;
//        cout << "! Average share visited nodes: " << ave_share << endl;
        cout << "! Average number of read partitions: " << ave_p_num << endl;
        cout << "Average proportion of io time: " << 100*ave_io_time/ave_time << " %" << endl;
        return query_time;
    }
    //function of initiating One-Pass algorithms
    void EMMCGraph::OnePassInitiation_IO(int algo_choice) {
        //initialise variables for one-pass algorithms
        mc_finished.assign(num_of_cri, false);
        mc_min_cost.assign(num_of_cri, INF);
        share_number = 0;//update in functions share_and_slack

        if(algo_choice == BiMultiHops){
            EMEdgesIndex_Bi.assign(node_num, make_tuple(false,false,-1,-1));
        }else{
            EMEdgesIndex.assign(node_num,make_tuple(false,-1,-1));
        }

        //for algorithms with io optimization
        if(algo_choice == OneHop || algo_choice == MultiHops || algo_choice == BiMultiHops){
//            mc_HotPool.init(partition_number);
            if(algo_choice == OneHop || algo_choice == MultiHops){
//                mc_HotPool.set_capacity(Partition_N);
                vertex_cri.assign(node_num,-1);
            }else if(algo_choice == BiMultiHops){
//                mc_HotPool.set_capacity(Partition_N_Bi);
                vertex_cri_Bi.assign(node_num, make_pair(-1,-1));
            }
        }

    }
    //function of initiating One-Pass algorithms
    void EMMCGraph::OnePassClear_IO(int algo_choice) {
        set_cri.clear();
        cluster_to_IO.clear();
        EMEdgesIndex.clear();
        vertex_cri.clear();
        //for algorithms with io optimization
        if(algo_choice == OneHop || algo_choice == MultiHops || algo_choice == BiMultiHops){
            //clear
//            mc_HotPool.clear();
            if(algo_choice == BiMultiHops){
                EMEdgesIndex_Bi.clear();
                vertex_cri_Bi.clear();
            }
        }
    }

    //One-hop algorithm without io optimization (LRU+map+stxxl)
    void EMMCGraph::EM_MC_OneHop_NoIO_new(int node_start, int node_end){
        if (node_start == node_end) {
            cout << "Source vertex and Target vertex are the same!!!" << endl;
            return;
        }
        if(ifshow){
            cout<<"LRU+unordered_map+stxxl."<<endl;
            ifshow= false;
        }
        //statistics of IO
        Stats stats_b, stats_e;
        stats_b.start();
        //variables
        NodeId item_id;
        NodeId temp_id;
        Distance item_dis;
        Distance temp_dis;
        int cri_i;
        int temp_number;
        int partition_id;

//        VectorMCEdgesEMTuple_NoIO EMMCEdges;//the external vectors of different criteria
        vector<vector<Distance>> mc_cost(num_of_cri, vector<Distance>(node_num, INF));   //vector of cost
        //vector<vector<NodeId>> mc_pre(num_of_cri, vector<int>(node_num, -1));      //vector of predecessor id
        vector<vector<bool>> mc_closed(num_of_cri, vector<bool>(node_num));    //flag of vertex closing
        vector<int> cluster_to_cri(partition_number,0);//used to map which criterion read this partition
        set<int> set_remain;//set for storing the id of unshared criteria

        HotPool4<VectorMCEdgesEMTuple4_NoIO> myHotPool(partition_number,PartitionNumber_NoIO);

        Timer tt;
        tt.start();

        //Initiation of start node
        for (int i = 0; i < num_of_cri; ++i) {
            set_cri.insert(i);
            mc_cost[i][node_start] = 0;//cost of start node
            EM_MC_PQueue[i]->push(VertexCost(node_start, 0));
        }
        //Iteration
        while (!EM_JudgeEmpty(EM_MC_PQueue)) {//for every node in pqueue !mc_finished[0] || !mc_finished[1] || !mc_finished[2]
            set_remain = set_cri;
            //EM_OneHop_strategy_NoIO(set_remain,node_end,EM_MC_PQueue);//divide_and_slack
            while (!set_remain.empty()) {
                cri_i = *set_remain.begin();
//                while (mc_closed[cri_i][EM_MC_PQueue[cri_i]->top().id]) {//if found in visited set
//                    EM_MC_PQueue[cri_i]->pop();
//                }
                while (mc_closed[cri_i][EM_MC_PQueue[cri_i]->top().id] && !EM_MC_PQueue[cri_i]->empty()) {//if found in visited set
                    EM_MC_PQueue[cri_i]->pop();
                }
                if(EM_MC_PQueue[cri_i]->empty()){//if empty
                    mc_finished[cri_i] = true;
                    set_cri.erase(cri_i);
                    mc_min_cost[cri_i] = INF;
                    set_remain.erase(cri_i);
                    continue;
                }
                item_id = EM_MC_PQueue[cri_i]->top().id;// top min item
                item_dis = EM_MC_PQueue[cri_i]->top().cost;// top min item
                assert(item_id < node_num);
                EM_MC_PQueue[cri_i]->pop();
                if (item_id == node_end) {//if reach target node, this criteria end.
                    mc_finished[cri_i] = true;
                    set_cri.erase(cri_i);
                    mc_min_cost[cri_i] = mc_cost[cri_i][item_id];
                    ++visited_number;
//                    cout<<"Minimal cost of criterion " << cri_i << ": "<<mc_min_cost[cri_i]<<endl;
                }
                else {
                    partition_id = node_to_cluster[item_id];
                    if(myHotPool.clusterStatus[partition_id] == 'D'){//if in disk
                        GraphReadCluster4(myHotPool,partition_id);
                        ++p_num;
                    }
                    if(myHotPool.clusterStatus[partition_id] == 'I'){
                        ++visited_number;

                        /// read from EMMCEdges
//                        auto adjlist = myHotPool.HotPools[partition_id].equal_range(item_id);
//                        for(auto it = adjlist.first; it != adjlist.second; ++it){
                        for(auto it = myHotPool.HotPools[partition_id][item_id].begin(); it != myHotPool.HotPools[partition_id][item_id].end(); ++it){
                            temp_id = it->ID2;
                            if (it->getW(cri_i)==INF100)
                                continue;
                            if (mc_closed[cri_i][temp_id])//if closed
                                continue;
                            temp_dis = item_dis + it->getW(cri_i);
                              if (mc_cost[cri_i][temp_id] > temp_dis) {//slack operation
                                mc_cost[cri_i][temp_id] = temp_dis;
                                //mc_pre[cri_i][temp_id] = item_id;
                                EM_MC_PQueue[cri_i]->push(VertexCost(temp_id, temp_dis));
                            }
                        }
                        //set closed
                        mc_closed[cri_i][item_id] = true;
                    }else if(myHotPool.clusterStatus[partition_id] == 'S'){
                        ++visited_number;
                        for (int i = get<1>(EMEdgesIndex[item_id]); i < get<2>(EMEdgesIndex[item_id]); ++i) {
                            assert(myHotPool.MCEdges_EM[i].ID1 == item_id);
                            temp_id = myHotPool.MCEdges_EM[i].ID2;
                            if (myHotPool.MCEdges_EM[i].getW(cri_i)==INF100)
                                continue;
                            if (mc_closed[cri_i][temp_id])//if closed
                                continue;
                            temp_dis = item_dis + myHotPool.MCEdges_EM[i].getW(cri_i);
                            if (mc_cost[cri_i][temp_id] > temp_dis) {//slack operation
                                mc_cost[cri_i][temp_id] = temp_dis;
                                //mc_pre[cri_i][temp_id] = item_id;
                                EM_MC_PQueue[cri_i]->push(VertexCost(temp_id, temp_dis));
                            }
                        }
                        //set closed
                        mc_closed[cri_i][item_id] = true;
                    }
                }
                set_remain.erase(cri_i);
            }
        }
        tt.stop();
        query_time += tt.GetRuntime() * 1000;
        //io record
        stats_e = stats_b.get_stats();
        if (stats_e.m_noOfWrites > 0 || stats_e.m_noOfReads > 0) {
            cout << "There are " << stats_e.m_noOfWrites + stats_e.m_noOfReads << " io for stxxl." << endl;
            //cout<<stats_e<<endl;
            EM_IORecord(stats_e);
        }
    }
    //Multi-hops algorithm without io optimization (LRU+map+stxxl)
    void EMMCGraph::EM_MC_MultiHop_NoIO_new(int node_start, int node_end){
        if (node_start == node_end) {
            cout << "Source vertex and Target vertex are the same!!!" << endl;
            return;
        }
        if(ifshow){
            cout<<"LRU+unordered_map+stxxl."<<endl;
            ifshow= false;
        }
        //statstics of IO
        Stats stats_b, stats_e;
        stats_b.start();
        //variables

//        VectorMCEdgesEMTuple_NoIO EMMCEdges;
        vector<vector<Distance>> mc_cost(num_of_cri, vector<Distance>(node_num, INF));   //vector of cost
        //vector<vector<NodeId>> mc_pre(num_of_cri, vector<int>(node_num, -1));      //vector of predecessor id
        vector<vector<bool>> mc_closed(num_of_cri, vector<bool>(node_num));    //flag of vertex closing
        vector<int> cluster_to_cri(partition_number,0);//used to map which criterion read this partition
        bool flag_finish = false;
        bool flag_mcFinished = false;
        int partition_id;
        int cri_i;
        NodeId item_id;
        Distance item_dis;
        NodeId temp_id;
        Distance temp_dis;
        int temp_number;
        int temp_cluster_id;
        benchmark::heap<2, int, int> cri_pqueue(num_of_cri);
        vector<int> cri_hops(num_of_cri, 0);
        vector<int> cri_to_cluster(num_of_cri, 0);
        set<int> set_remain;//set for storing the id of unshared criteria

        //hot pool
        HotPool4<VectorMCEdgesEMTuple4_NoIO> myHotPool(partition_number,PartitionNumber_NoIO);

        Timer tt;
        tt.start();

        //Initiation of start node
        for (int i = 0; i < num_of_cri; ++i) {
            set_cri.insert(i);
            mc_cost[i][node_start] = 0;//cost of start node
            EM_MC_PQueue[i]->push(VertexCost(node_start,0));
        }
        //first read
//        ReadClusterToExternalVectorMC(node_to_cluster[node_start],EMMCEdges);
//        GraphMCReadCluster_NoIO(myHotPool,node_to_cluster[node_start]);
        GraphReadCluster4(myHotPool,node_to_cluster[node_start]);
        ++p_num;
        cluster_to_cri[node_to_cluster[node_start]] = 0;
        //Iteration
        while (!EM_JudgeEmpty(EM_MC_PQueue)) {//for every node in pqueue !mc_finished[0] || !mc_finished[1] || !mc_finished[2]
            set_remain = set_cri;
            //EM_MultiHop_strategy_NoIO(set_remain,node_end, EM_MC_PQueue,EMMCEdges);//divide_and_slack
            while (!set_remain.empty()) {
                cri_i = *set_remain.begin();

                while (mc_closed[cri_i][EM_MC_PQueue[cri_i]->top().id]) {//if found in visited set
                    EM_MC_PQueue[cri_i]->pop();
                    if(EM_MC_PQueue[cri_i]->empty()){
                        mc_finished[cri_i] = true;
                        set_cri.erase(cri_i);
                        mc_min_cost[cri_i] = INF;
                        cri_hops[cri_i] = INF;
                        cri_pqueue.update(cri_i, INF);
                        set_remain.erase(cri_i);
                        goto jump1;
                    }
                }
                item_id = EM_MC_PQueue[cri_i]->top().id;// top min item
                item_dis = EM_MC_PQueue[cri_i]->top().cost;

                partition_id = node_to_cluster[item_id];

                while (myHotPool.clusterStatus[partition_id]!='D') {//if found get<1>(EMEdgesIndex[item_id]) != -1 myHotPool.FlagIM[partition_id]
                    if (item_id == node_end) {//if reach target node, this criterion ends.
                        //                mc_pqueue[cri_i].clear();
                        mc_finished[cri_i] = true;
                        set_cri.erase(cri_i);
                        mc_min_cost[cri_i] = mc_cost[cri_i][item_id];
//                        cout<<"Minimal cost of criterion " << cri_i << ": "<<mc_min_cost[cri_i]<<endl;
                        flag_finish = true;
                        cri_hops[cri_i] = INF;
                        cri_pqueue.update(cri_i, INF);
                        ++visited_number;
                        flag_mcFinished = true;
                        for (auto it = mc_finished.begin(); it != mc_finished.end(); ++it) {
                            if (!*it) {//if none of them is false, then flag_mcFinished = true
                                flag_mcFinished = false;
                                break;
                            }
                        }
                        break;
                    }
                    else {
                        ++cri_hops[cri_i];
                        //slack
                        /// io-efficient graph data reading
                        if (cluster_to_cri[node_to_cluster[item_id]] != cri_i) {//if the cluster is not read by this criterion
                            ++share_number;
                        }
                        //set closed
                        mc_closed[cri_i][item_id] = true;
                        EM_MC_PQueue[cri_i]->pop();//pop min item
                        ++visited_number;//update #visited_number

                        if(myHotPool.clusterStatus[partition_id] == 'I'){
                            /// read from EMMCEdges
//                            auto adjlist = myHotPool.HotPools[partition_id].equal_range(item_id);
//                            for(auto it = adjlist.first; it != adjlist.second; ++it){
                            for(auto it = myHotPool.HotPools[partition_id][item_id].begin(); it != myHotPool.HotPools[partition_id][item_id].end(); ++it){
                                temp_id = it->ID2;
                                if (it->getW(cri_i)==INF100)
                                    continue;
                                if (mc_closed[cri_i][temp_id])//if closed
                                    continue;
                                temp_dis = item_dis + it->getW(cri_i);
                                if (mc_cost[cri_i][temp_id] > temp_dis) {//slack operation
                                    mc_cost[cri_i][temp_id] = temp_dis;
                                    //mc_pre[cri_i][temp_id] = item_id;
                                    EM_MC_PQueue[cri_i]->push(VertexCost(temp_id, temp_dis));
                                }
                            }
                            //set closed
                            mc_closed[cri_i][item_id] = true;
                        }else if(myHotPool.clusterStatus[partition_id] == 'S'){
                            for (int i = get<1>(EMEdgesIndex[item_id]); i < get<2>(EMEdgesIndex[item_id]); ++i) {
                                assert(myHotPool.MCEdges_EM[i].ID1 == item_id);
                                temp_id = myHotPool.MCEdges_EM[i].ID2;
                                if (myHotPool.MCEdges_EM[i].getW(cri_i)==INF100)
                                    continue;
                                if (mc_closed[cri_i][temp_id])//if closed
                                    continue;
                                temp_dis = item_dis + myHotPool.MCEdges_EM[i].getW(cri_i);
                                if (mc_cost[cri_i][temp_id] > temp_dis) {//slack operation
                                    mc_cost[cri_i][temp_id] = temp_dis;
                                    //mc_pre[cri_i][temp_id] = item_id;
                                    EM_MC_PQueue[cri_i]->push(VertexCost(temp_id, temp_dis));
                                }
                            }
                            //set closed
                            mc_closed[cri_i][item_id] = true;
                        }

                        while (mc_closed[cri_i][EM_MC_PQueue[cri_i]->top().id] && !EM_MC_PQueue[cri_i]->empty()) {//if already closed
                            EM_MC_PQueue[cri_i]->pop();
                        }
                        if(EM_MC_PQueue[cri_i]->empty()){
                            mc_finished[cri_i] = true;
                            set_cri.erase(cri_i);
                            mc_min_cost[cri_i] = INF;
                            cri_hops[cri_i] = INF;
                            cri_pqueue.update(cri_i, INF);
                            set_remain.erase(cri_i);
                            goto jump1;
                        }
                        item_id = EM_MC_PQueue[cri_i]->top().id;// top min item
                        item_dis = EM_MC_PQueue[cri_i]->top().cost;

                        partition_id = node_to_cluster[item_id];
                    }
                }
                set_remain.erase(cri_i);
                if (!flag_finish) {
                    cluster_to_cri[node_to_cluster[item_id]] = cri_i;
                    cri_to_cluster[cri_i] = node_to_cluster[item_id];
                    cri_pqueue.update(cri_i, cri_hops[cri_i]);
                }
                else {
                    flag_finish = false;
                }
                if (set_remain.empty() && !flag_mcFinished) {//if all criteria are processed &&!EM_JudgeEmpty(em_mc_pqueue)
                    /*** read by synchronizing the hops ***/
                    temp_cluster_id = cri_to_cluster[cri_pqueue.top_value()];
                    assert(get<1>(EMEdgesIndex[partitions_info[temp_cluster_id]]) == -1);
//                    ReadClusterToExternalVectorMC(temp_cluster_id, EMMCEdges);
//                    GraphMCReadCluster_NoIO(myHotPool,temp_cluster_id);
                    GraphReadCluster4(myHotPool,temp_cluster_id);
                    ++p_num;
                }
                jump1: "there is a jump.";
            }
        }
        tt.stop();
        query_time += tt.GetRuntime() * 1000;
        //io record
        stats_e = stats_b.get_stats();
        if (stats_e.m_noOfWrites > 0 || stats_e.m_noOfReads > 0) {
            cout << "There are " << stats_e.m_noOfWrites + stats_e.m_noOfReads << " io for stxxl." << endl;
            //cout<<stats_e<<endl;
            EM_IORecord(stats_e);
        }
        //clear
//        EMMCEdges.clear();
    }
    //BiMulti-hops algorithm without io optimization (LRU+map+stxxl)
    void EMMCGraph::EM_MC_BiMultiHop_NoIO_new(int node_start, int node_end) {
        if (node_start == node_end) {
            cout << "Source vertex and Target vertex are the same!!!" << endl;
            return;
        }
        if(ifshow){
            cout<<"LRU+unordered_map+stxxl."<<endl;
            ifshow= false;
        }
        //statstics of IO
        Stats stats_b, stats_e;
        stats_b.start();
        //variables
        set<int> set_remain;//set for storing the id of unshared criteria

//        VectorMCEdgesEMTuple_NoIO_Bi EMMCEdges;
        vector<vector<Distance>> mc_cost(num_of_cri, vector<Distance>(node_num, INF));   //vector of cost
        //vector<vector<NodeId>> mc_pre(num_of_cri, vector<int>(node_num, -1));      //vector of predecessor id
        vector<vector<bool>> mc_closed(num_of_cri, vector<bool>(node_num));    //flag of vertex closing
        vector<vector<Distance>> mc_cost_r(num_of_cri, vector<Distance>(node_num, INF));   //vector of cost
        //vector<vector<NodeId>> mc_pre_r(num_of_cri, vector<int>(node_num, -1));      //vector of predecessor id
        vector<vector<bool>> mc_closed_r(num_of_cri, vector<bool>(node_num));    //flag of vertex closing
        vector<int> cluster_to_cri(partition_number,0);//used to map which criterion read this partition
        vector<int8_t> cluster_to_bi(partition_number,0);   //used to indicate the partition is read by forward search or reverse search
        vector<NodeId> mc_terminate_id(num_of_cri,-1);;//termination id of bi-dijkstra
        NodeId item_id, item_id_r;
        NodeId temp_id, temp_id_r;
        Distance temp_dis, temp_dis_r;
        Distance item_dis, item_dis_r;
        int partition_id, partition_id_r;
        int cri_i;
        int temp_cluster_id, temp_cluster_id_r;
        bool temp_bool;
        int temp_degree;
        bool flag_mcFinished = false;
        benchmark::heap<2, int, int> cri_pqueue(num_of_cri);
        vector<int> cri_hops(num_of_cri, 0);
        vector<int> cri_to_cluster(num_of_cri, 0);
        benchmark::heap<2, int, int> cri_pqueue_r(num_of_cri);
        vector<int> cri_hops_r(num_of_cri, 0);
        vector<int> cri_to_cluster_r(num_of_cri, 0);
        //hot pool
        HotPool4<VectorMCEdgesEMTuple4_NoIO_Bi> myHotPool(partition_number,PartitionNumber_NoIO_Bi);

        Timer tt;
        tt.start();

        //Initiation of start node
        for (int i = 0; i < num_of_cri; i++) {
            set_cri.insert(i);
            mc_cost[i][node_start] = 0;//cost of start node
            EM_MC_PQueue[i]->push(VertexCost(node_start,0));
            mc_cost_r[i][node_end] = 0;
            EM_MC_PQueue_r[i]->push(VertexCost(node_end,0));
        }

        //first read
//        ReadClusterToExternalVectorMC(node_to_cluster[node_start],EMMCEdges);
//        GraphMCReadCluster_NoIO(myHotPool,node_to_cluster[node_start]);
        GraphReadCluster4(myHotPool,node_to_cluster[node_start]);
        ++p_num;
        cluster_to_cri[node_to_cluster[node_start]] = 0;
        cluster_to_bi[node_to_cluster[node_start]] = FORWARD;
        if(node_to_cluster[node_end]!=node_to_cluster[node_start]){
//            ReadClusterToExternalVectorMC(node_to_cluster[node_end],EMMCEdges);
//            GraphMCReadCluster_NoIO(myHotPool,node_to_cluster[node_end]);
            GraphReadCluster4(myHotPool,node_to_cluster[node_end]);
            ++p_num;
            cluster_to_cri[node_to_cluster[node_end]] = 0;
            cluster_to_bi[node_to_cluster[node_end]] = REVERSE;
        }

        //Iteration
        while (!EM_JudgeEmptyBi(EM_MC_PQueue, EM_MC_PQueue_r)) {//for every node in pqueue !mc_finished[0] || !mc_finished[1] || !mc_finished[2]
            set_remain = set_cri;

            while (!set_remain.empty()) {
                cri_i = *set_remain.begin();//pick one criterion to process
                //deal with invalid top elements
                while (mc_closed[cri_i][EM_MC_PQueue[cri_i]->top().id]) {
                    EM_MC_PQueue[cri_i]->pop();
                    if(EM_MC_PQueue[cri_i]->empty()){
                        mc_finished[cri_i] = true;
                        set_cri.erase(cri_i);
                        cri_hops[cri_i] = INF;
                        cri_pqueue.update(cri_i, INF);
                        cri_hops_r[cri_i] = INF;
                        cri_pqueue_r.update(cri_i, INF);
                        goto jump1;
                    }
                }
                while (mc_closed_r[cri_i][EM_MC_PQueue_r[cri_i]->top().id]) {
                    EM_MC_PQueue_r[cri_i]->pop();
                    if(EM_MC_PQueue[cri_i]->empty()){
                        mc_finished[cri_i] = true;
                        set_cri.erase(cri_i);
                        cri_hops[cri_i] = INF;
                        cri_pqueue.update(cri_i, INF);
                        cri_hops_r[cri_i] = INF;
                        cri_pqueue_r.update(cri_i, INF);
                        goto jump1;
                    }
                }

                //// Forward searching
                item_id = EM_MC_PQueue[cri_i]->top().id;// top min item
                item_dis = EM_MC_PQueue[cri_i]->top().cost;// top min item

                partition_id = node_to_cluster[item_id];
                while (myHotPool.clusterStatus[partition_id]!='D') {//if found get<1>(EMEdgesIndex[item_id]) != -1  && !mc_closed_r[cri_i][item_id]
                    if(mc_closed_r[cri_i][item_id] || item_id == node_end){
                        mc_finished[cri_i] = true;
                        set_cri.erase(cri_i);
                        if(item_id == node_end) {
                            mc_min_cost[cri_i] = mc_cost[cri_i][item_id];
                        }
//                            cout<<"Minimal cost of criterion " << cri_i << ": "<<mc_min_cost[cri_i]<<endl;
                        cri_hops[cri_i] = INF;
                        cri_pqueue.update(cri_i, INF);
                        cri_hops_r[cri_i] = INF;
                        cri_pqueue_r.update(cri_i, INF);
                        ++visited_number;
                        set_remain.erase(cri_i);
                        goto jump1;
                    }

                    ++cri_hops[cri_i];
                    if (cluster_to_cri[node_to_cluster[item_id]] != cri_i) {//if the cluster is not read by this criterion
                        ++share_number;
                    }
                    //set closed
                    mc_closed[cri_i][item_id] = true;
                    EM_MC_PQueue[cri_i]->pop();//pop min item
                    ++visited_number;//update #visited_number

                    if(myHotPool.clusterStatus[partition_id] == 'I'){
                        /// read from EMMCEdges
//                            auto adjlist = myHotPool.HotPools[partition_id].equal_range(item_id);
//                            for(auto it = adjlist.first; it != adjlist.second; ++it){
                        for(auto it = myHotPool.HotPools[partition_id][item_id].begin(); it != myHotPool.HotPools[partition_id][item_id].end(); ++it){
                            temp_id = it->ID2;
                            if (it->getW(cri_i)==INF100)
                                continue;
                            temp_dis = item_dis + it->getW(cri_i);
                            if (!mc_closed[cri_i][temp_id]) {//if not closed
                                if (mc_cost[cri_i][temp_id] > temp_dis) {//slack operation
                                    mc_cost[cri_i][temp_id] = temp_dis;
                                    //mc_pre[cri_i][temp_id] = item_id;
                                    EM_MC_PQueue[cri_i]->push(VertexCost(temp_id, temp_dis));
                                }
                            }
                            //update min cost
                            temp_bool = mc_closed_r[cri_i][temp_id];
                            if (temp_bool && temp_dis + mc_cost_r[cri_i][temp_id] < mc_min_cost[cri_i]) {
                                mc_min_cost[cri_i] = temp_dis + mc_cost_r[cri_i][temp_id];
                                mc_terminate_id[cri_i] = temp_id;
                            }
                        }
                        //set closed
                        mc_closed[cri_i][item_id] = true;
                    }else if(myHotPool.clusterStatus[partition_id] == 'S'){
                        for (int i = get<1>(EMEdgesIndex[item_id]); i < get<2>(EMEdgesIndex[item_id]); ++i) {
                            assert(myHotPool.MCEdges_EM[i].ID1 == item_id);
                            temp_id = myHotPool.MCEdges_EM[i].ID2;
                            if (myHotPool.MCEdges_EM[i].getW(cri_i)==INF100)
                                continue;
                            temp_dis = item_dis + myHotPool.MCEdges_EM[i].getW(cri_i);
                            if (!mc_closed[cri_i][temp_id]) {//if not closed
                                if (mc_cost[cri_i][temp_id] > temp_dis) {//slack operation
                                    mc_cost[cri_i][temp_id] = temp_dis;
                                    //mc_pre[cri_i][temp_id] = item_id;
                                    EM_MC_PQueue[cri_i]->push(VertexCost(temp_id, temp_dis));
                                }
                            }
                            //update min cost
                            temp_bool = mc_closed_r[cri_i][temp_id];
                            if (temp_bool && temp_dis + mc_cost_r[cri_i][temp_id] < mc_min_cost[cri_i]) {
                                mc_min_cost[cri_i] = temp_dis + mc_cost_r[cri_i][temp_id];
                                mc_terminate_id[cri_i] = temp_id;
                            }
                        }
                        //set closed
                        mc_closed[cri_i][item_id] = true;
                    }

                    while (mc_closed[cri_i][EM_MC_PQueue[cri_i]->top().id] && !EM_MC_PQueue[cri_i]->empty()) {//if already closed
                        EM_MC_PQueue[cri_i]->pop();
                    }
                    if(EM_MC_PQueue[cri_i]->empty()){
                        mc_finished[cri_i] = true;
                        set_cri.erase(cri_i);
                        mc_min_cost[cri_i] = INF;
                        cri_hops[cri_i] = INF;
                        cri_pqueue.update(cri_i, INF);
                        cri_hops_r[cri_i] = INF;
                        cri_pqueue_r.update(cri_i, INF);
                        set_remain.erase(cri_i);
                        goto jump1;
                    }

                    item_id = EM_MC_PQueue[cri_i]->top().id;// top min item
                    item_dis = EM_MC_PQueue[cri_i]->top().cost;// top min item

                    partition_id = node_to_cluster[item_id];
                }
                cluster_to_cri[node_to_cluster[item_id]] = cri_i;
                cri_to_cluster[cri_i] = node_to_cluster[item_id];
                cri_pqueue.update(cri_i, cri_hops[cri_i]);
                //// Reverse searching
                item_id_r = EM_MC_PQueue_r[cri_i]->top().id;// top min item
                item_dis_r = EM_MC_PQueue_r[cri_i]->top().cost;// top min item

                partition_id_r = node_to_cluster[item_id_r];
                while (myHotPool.clusterStatus[partition_id_r]!='D') {//if found get<1>(EMEdgesIndex[item_id_r]) != -1  && !mc_closed[cri_i][item_id_r]
                    if(mc_closed[cri_i][item_id_r] || item_id_r == node_start){
                        mc_finished[cri_i] = true;
                        set_cri.erase(cri_i);
                        if(item_id_r == node_start){
                            mc_min_cost[cri_i] = mc_cost_r[cri_i][item_id_r];
                        }
//                            cout<<"Minimal cost of criterion " << cri_i << ": "<<mc_min_cost[cri_i]<<endl;
                        cri_hops[cri_i] = INF;
                        cri_pqueue.update(cri_i, INF);
                        cri_hops_r[cri_i] = INF;
                        cri_pqueue_r.update(cri_i, INF);
                        ++visited_number;
                        set_remain.erase(cri_i);
                        goto jump1;
                    }
                    ++cri_hops_r[cri_i];
                    if (cluster_to_cri[node_to_cluster[item_id_r]] != cri_i) {//if the cluster is not read by this criterion
                        ++share_number;
                    }
                    //set closed
                    mc_closed_r[cri_i][item_id_r] = true;
                    EM_MC_PQueue_r[cri_i]->pop();//pop min item
                    ++visited_number;//update #visited_number
                    //slack

                    if(myHotPool.clusterStatus[partition_id_r] == 'I'){
                        /// read from EMMCEdges
//                            auto adjlist_r = myHotPool.HotPools[partition_id_r].equal_range(item_id_r);
//                            for(auto it = adjlist_r.first; it != adjlist_r.second; ++it){
                        for(auto it = myHotPool.HotPools[partition_id_r][item_id_r].begin(); it != myHotPool.HotPools[partition_id_r][item_id_r].end(); ++it){
                            temp_id_r = it->ID2;
                            if (it->getW(cri_i)==INF100)
                                continue;
                            temp_dis_r = item_dis_r + it->getW(cri_i);
                            if (!mc_closed_r[cri_i][temp_id_r]){//if closed
                                if (mc_cost_r[cri_i][temp_id_r] > temp_dis_r) {//slack operation
                                    mc_cost_r[cri_i][temp_id_r] = temp_dis_r;
                                    //mc_pre_r[cri_i][temp_id_r] = item_id_r;
                                    EM_MC_PQueue_r[cri_i]->push(VertexCost(temp_id_r, temp_dis_r));
                                }
                            }
                            //update min cost
                            temp_bool = mc_closed[cri_i][temp_id_r];
                            if (temp_bool && temp_dis_r + mc_cost[cri_i][temp_id_r] < mc_min_cost[cri_i]) {
                                mc_min_cost[cri_i] = temp_dis_r + mc_cost[cri_i][temp_id_r];
                                mc_terminate_id[cri_i] = temp_id_r;
                            }
                        }
                        //set closed
                        mc_closed_r[cri_i][item_id_r] = true;
                    }else if(myHotPool.clusterStatus[partition_id_r] == 'S'){
                        for (int i = get<1>(EMEdgesIndex[item_id_r]); i < get<2>(EMEdgesIndex[item_id_r]); ++i) {
                            assert(myHotPool.MCEdges_EM[i].ID1 == item_id_r);
                            temp_id_r = myHotPool.MCEdges_EM[i].ID2;
                            if (myHotPool.MCEdges_EM[i].getW(cri_i)==INF100)
                                continue;
                            temp_dis_r = item_dis_r + myHotPool.MCEdges_EM[i].getW(cri_i);
                            if (!mc_closed_r[cri_i][temp_id_r]){//if closed
                                if (mc_cost_r[cri_i][temp_id_r] > temp_dis_r) {//slack operation
                                    mc_cost_r[cri_i][temp_id_r] = temp_dis_r;
                                    //mc_pre_r[cri_i][temp_id_r] = item_id_r;
                                    EM_MC_PQueue_r[cri_i]->push(VertexCost(temp_id_r, temp_dis_r));
                                }
                            }
                            //update min cost
                            temp_bool = mc_closed[cri_i][temp_id_r];
                            if (temp_bool && temp_dis_r + mc_cost[cri_i][temp_id_r] < mc_min_cost[cri_i]) {
                                mc_min_cost[cri_i] = temp_dis_r + mc_cost[cri_i][temp_id_r];
                                mc_terminate_id[cri_i] = temp_id_r;
                            }
                        }
                        //set closed
                        mc_closed_r[cri_i][item_id_r] = true;
                    }
                    while (mc_closed_r[cri_i][EM_MC_PQueue_r[cri_i]->top().id] && !EM_MC_PQueue_r[cri_i]->empty()) {
                        EM_MC_PQueue_r[cri_i]->pop();
                    }
                    if(EM_MC_PQueue_r[cri_i]->empty()){
                        mc_finished[cri_i] = true;
                        set_cri.erase(cri_i);
                        mc_min_cost[cri_i] = INF;
                        cri_hops[cri_i] = INF;
                        cri_pqueue.update(cri_i, INF);
                        cri_hops_r[cri_i] = INF;
                        cri_pqueue_r.update(cri_i, INF);
                        set_remain.erase(cri_i);
                        goto jump1;
                    }
                    item_id_r = EM_MC_PQueue_r[cri_i]->top().id;// top min item
                    item_dis_r = EM_MC_PQueue_r[cri_i]->top().cost;// top min item

                    partition_id_r = node_to_cluster[item_id_r];
                }
                cri_to_cluster_r[cri_i] = node_to_cluster[item_id_r];
                cluster_to_cri[cri_to_cluster_r[cri_i]] = cri_i;//node_to_cluster[item_id_r]
                cri_pqueue_r.update(cri_i, cri_hops_r[cri_i]);

                //// Termination judging
                if (EM_MC_PQueue[cri_i]->top().cost + EM_MC_PQueue_r[cri_i]->top().cost >= mc_min_cost[cri_i]) {//condition of termination
                    mc_finished[cri_i] = true;
                    set_cri.erase(cri_i);
                    cri_hops[cri_i] = INF;
                    cri_pqueue.update(cri_i, INF);
                    cri_hops_r[cri_i] = INF;
                    cri_pqueue_r.update(cri_i, INF);
//                    cout<<"Minimal cost of criterion "<<cri_i<<": "<<mc_min_cost[cri_i]<<endl;
                    flag_mcFinished = true;
                    for (auto it = mc_finished.begin(); it != mc_finished.end(); ++it) {
                        if (!*it) {//if none of them is false, then flag_mcFinished = true
                            flag_mcFinished = false;
                            break;
                        }
                    }
                }
                set_remain.erase(cri_i);
                // judge reading
                if (set_remain.empty() && !flag_mcFinished) {//if all criteria are processed and not all criteria are finished
                    /// read by synchronizing the hops
                    temp_cluster_id = cri_to_cluster[cri_pqueue.top_value()];
                    temp_cluster_id_r = cri_to_cluster_r[cri_pqueue_r.top_value()];

                    if(cri_hops[cri_pqueue.top_value()] <= cri_hops_r[cri_pqueue_r.top_value()]){
                        if (cluster_to_bi[temp_cluster_id] == 0 || get<1>(EMEdgesIndex[partitions_info[temp_cluster_id]]) == -1) {
//                        ReadClusterToExternalVectorMC(temp_cluster_id, EMMCEdges);
//                        GraphMCReadCluster_NoIO(myHotPool,temp_cluster_id);
                            GraphReadCluster4(myHotPool,temp_cluster_id);
                            cluster_to_bi[temp_cluster_id] = FORWARD;
                            ++p_num;
                        }
                    }else{
                        if (cluster_to_bi[temp_cluster_id_r] == 0 || get<1>(EMEdgesIndex[partitions_info[temp_cluster_id_r]]) == -1) {
//                            if (temp_cluster_id != temp_cluster_id_r) {//if different
//                            ReadClusterToExternalVectorMC(temp_cluster_id_r, EMMCEdges);
//                            GraphMCReadCluster_NoIO(myHotPool,temp_cluster_id_r);
                                GraphReadCluster4(myHotPool,temp_cluster_id_r);
                                cluster_to_bi[temp_cluster_id_r] = REVERSE;
                                ++p_num;
//                            }
                        }
                    }

                }
                jump1: "there is a jump.";
            }
        }
        tt.stop();
        query_time += tt.GetRuntime() * 1000;
        //io record
        stats_e = stats_b.get_stats();
        if (stats_e.m_noOfWrites > 0 || stats_e.m_noOfReads > 0) {
            cout << "There are " << stats_e.m_noOfWrites + stats_e.m_noOfReads << " io for stxxl." << endl;
            //cout<<stats_e<<endl;
            EM_IORecord(stats_e);
        }
    }
    //Function for one-stop algorithm, i.e. BASE
    void EMMCGraph::EM_MC_OneHop(int node_start, int node_end){
        if (node_start == node_end) {
            cout << "Source vertex and Target vertex are the same!!!" << endl;
            return;
        }
        //statistics of IO
        Stats stats_b, stats_e;
        stats_b.start();
        //variables
        vector<vector<Distance>> mc_cost(num_of_cri, vector<Distance>(node_num, INF));   //vector of cost
        //vector<vector<NodeId>> mc_pre(num_of_cri, vector<int>(node_num, -1));      //vector of predecessor id
        vector<vector<bool>> mc_closed(num_of_cri, vector<bool>(node_num));    //flag of vertex closing
        vector<int> cluster_to_cri(partition_number,0);//used to map which criterion read this partition
        NodeId item_id;
        NodeId temp_id;
        Distance temp_dis;
        Distance item_dis;
        int cri_i;
        int temp_number;
        int temp_degree = 0;
        int partition_id;
        set<int> set_remain;//set for storing the id of unshared criteria
        int cri_empty = 0;
        partition_left = 0;
        //hot pool
        MCHotPool<VectorMCEdgesEMTuple_IO> mc_HotPool(Partition_N);       //Hot pool for one-pass algorithms
        mc_HotPool.init(partition_number);
        mc_HotPool.set_capacity(Partition_N);

        Timer tt;
        tt.start();

        //Initiation of start node
        for (int i = 0; i < num_of_cri; ++i) {
            set_cri.insert(i);
            mc_cost[i][node_start] = 0;//cost of start node
            EM_MC_PQueue[i]->push(VertexCost(node_start, 0));
        }
        //Iteration
        while (!EM_JudgeEmpty(EM_MC_PQueue)) {//for every node in pqueue !mc_finished[0] || !mc_finished[1] || !mc_finished[2]
            set_remain = set_cri;
            //EM_OneHop_strategy_New(mc_HotPool,set_remain,node_end,EM_MC_PQueue);//divide_and_slack
            while (!set_remain.empty()) {
                cri_i = *set_remain.begin();
                while (mc_closed[cri_i][EM_MC_PQueue[cri_i]->top().id] && !EM_MC_PQueue[cri_i]->empty()) {//if found in visited set
                    EM_MC_PQueue[cri_i]->pop();
                }
                if(EM_MC_PQueue[cri_i]->empty()){
                    mc_finished[cri_i] = true;
                    set_cri.erase(cri_i);
                    mc_min_cost[cri_i] = INF;
                    set_remain.erase(cri_i);
                    ++cri_empty;
                    continue;
                }
                item_id = EM_MC_PQueue[cri_i]->top().id;// top min item
                item_dis = EM_MC_PQueue[cri_i]->top().cost;// top min item
                EM_MC_PQueue[cri_i]->pop();

                assert(item_id < node_num);
                if (item_id == node_end) {//if reach target node, this criteria end.
                    mc_finished[cri_i] = true;
                    set_cri.erase(cri_i);
                    mc_min_cost[cri_i] = mc_cost[cri_i][item_id];
                    ++visited_number;
//                    cout<<"Minimal cost of criterion " << cri_i << ": "<<mc_min_cost[cri_i]<<endl;
                }
                else {
                    //slack
                    /// read from mc_HotPool
                    partition_id = node_to_cluster[item_id];
                    switch (get<0>(mc_HotPool.clusterStatus[partition_id])) {// deal according to partition storage status
                    case 'D': //if in disk, read corresponding partition
                        GraphMCReadCluster_New(mc_HotPool, item_id, partition_id, OneHop);
                    case 'I': {//if in memory
                        if (cluster_to_cri[partition_id] != cri_i) {//if the cluster is not read by this criterion
                            ++share_number;
                        }
                        //set closed
                        mc_closed[cri_i][item_id] = true;
                        ++visited_number;//update #visited_number
                        temp_degree = 0;
                        for (int i = get<1>(EMEdgesIndex[item_id]); i < get<2>(EMEdgesIndex[item_id]); ++i) {
                            assert(mc_HotPool.HotPools[partition_id][i].ID1 == item_id);
                            ++temp_degree;
                            temp_id = mc_HotPool.HotPools[partition_id][i].ID2;
                            if (mc_HotPool.HotPools[partition_id][i].getW(cri_i)==INF100)
                                continue;
                            if (mc_closed[cri_i][temp_id])//if closed
                                continue;
                            temp_dis = item_dis + mc_HotPool.HotPools[partition_id][i].getW(cri_i);
                            if (mc_cost[cri_i][temp_id] > temp_dis) {//slack operation
                                mc_cost[cri_i][temp_id] = temp_dis;
                                //mc_pre[cri_i][temp_id] = item_id;
                                EM_MC_PQueue[cri_i]->push(VertexCost(temp_id, temp_dis));
                            }
                        }
                        break;
                    }
                    case 'S': {//if in external vector
                        assert(get<2>(EMEdgesIndex[item_id]) > 0);
                        if (cluster_to_cri[partition_id] != cri_i) {//if the cluster is not read by this criterion
                            ++share_number;
                        }
                        //set closed
                        mc_closed[cri_i][item_id] = true;
                        ++visited_number;//update #visited_number
                        temp_degree = 0;
                        for (int i = get<1>(EMEdgesIndex[item_id]); i < get<2>(EMEdgesIndex[item_id]); ++i) {
                            assert(mc_HotPool.MCEdges_EM[i].ID1 == item_id);
                            ++temp_degree;
                            temp_id = mc_HotPool.MCEdges_EM[i].ID2;
                            if (mc_HotPool.MCEdges_EM[i].getW(cri_i)==INF100)
                                continue;
                            if (mc_closed[cri_i][temp_id])//if closed
                                continue;
                            temp_dis = item_dis + mc_HotPool.MCEdges_EM[i].getW(cri_i);
                            if (mc_cost[cri_i][temp_id] > temp_dis) {//slack operation
                                mc_cost[cri_i][temp_id] = temp_dis;
                                //mc_pre[cri_i][temp_id] = item_id;
                                EM_MC_PQueue[cri_i]->push(VertexCost(temp_id, temp_dis));
                            }
                        }
                        break;
                    }
                    default: {
                        cout << "Wrong partition status!!!" << endl; break;
                    }
                    }
                    --vertex_cri[item_id];
                    assert(vertex_cri[item_id] >= 0);
                    //// evict vertex
                    if (vertex_cri[item_id]-cri_empty <= 0) {//evict valid in-memory vertices immediately ==
                        if (get<0>(mc_HotPool.clusterStatus[partition_id]) == 'I') {
                            get<0>(EMEdgesIndex[item_id]) = true;
                            get<2>(mc_HotPool.clusterStatus[partition_id]) -= temp_degree;
                            if (get<2>(mc_HotPool.clusterStatus[partition_id]) <= 0) {//if empty, clear it immediately
                                assert(get<2>(mc_HotPool.clusterStatus[partition_id]) == 0);
                                for(int i=partitions_info[partition_id];i<partitions_info[partition_id+1];++i){
                                    get<1>(EMEdgesIndex[i]) = -1;//hard delete
                                    get<2>(EMEdgesIndex[i]) = -1;
                                }
//                                for(auto it=mc_HotPool.HotPools[partition_id].begin();it!=mc_HotPool.HotPools[partition_id].end();++it){
//                                    get<1>(EMEdgesIndex[it->ID1]) = -1;//hard delete
//                                    get<2>(EMEdgesIndex[it->ID1]) = -1;
//                                }
                                mc_HotPool.HotPools[partition_id].clear();
                                mc_HotPool.HotPools[partition_id].shrink_to_fit();//free memory
//                                mc_HotPool.im_partitions.erase(partition_id);
                                mc_HotPool.LRU_IMCluster.erase(partition_id);
                                get<0>(mc_HotPool.clusterStatus[partition_id]) = 'D';
                            }
                        }
                    }
                }
                set_remain.erase(cri_i);
            }
        }
        tt.stop();
        query_time += tt.GetRuntime() * 1000;
//        cout<<"Query time: "<<tt.GetRuntime()<<" s."<<endl;
        //io record
        stats_e = stats_b.get_stats();
        if (stats_e.m_noOfWrites > 0 || stats_e.m_noOfReads > 0) {
            cout << "There are " << stats_e.m_noOfWrites + stats_e.m_noOfReads << " io for stxxl." << endl;
            //cout<<stats_e<<endl;
            EM_IORecord(stats_e);
        }
        mc_HotPool.clear();
    }
    //function of reading partition for one-pass algorithms with io optimization: version 2, deal with second read
    void EMMCGraph::GraphMCReadCluster_New(MCHotPool<VectorMCEdgesEMTuple_IO> & mcHotPool,NodeId ID1,int target_p_id,int algo_choice){
        int evict_p_id = mcHotPool.LRU_IMCluster.put(target_p_id,1);
        //// evict old partition for new one
        if(evict_p_id != -1){//if it is necessary to evict old partition
            assert(get<0>(mcHotPool.clusterStatus[evict_p_id]) == 'I');
            unsigned int evict_p_size = get<2>(mcHotPool.clusterStatus[evict_p_id]);
            uint temp_sz = 0;
            if(algo_choice == OneHop){
                temp_sz = get<1>(mcHotPool.clusterStatus[evict_p_id])*alpha;
            }else if(algo_choice == MultiHops){
                temp_sz = get<1>(mcHotPool.clusterStatus[evict_p_id])*alpha_multi;
            }
            assert(temp_sz>0);
            assert(evict_p_size>0);

            if(evict_p_size > temp_sz){//if the remaining size is too large or no node remain, evict all
                get<2>(mcHotPool.clusterStatus[evict_p_id]) = 0;
                get<0>(mcHotPool.clusterStatus[evict_p_id]) = 'D';
                for(int i=partitions_info[evict_p_id];i<partitions_info[evict_p_id+1];++i){
                    get<1>(EMEdgesIndex[i]) = -1;//hard delete
                    get<2>(EMEdgesIndex[i]) = -1;
                }
//                for(auto it=mcHotPool.HotPools[evict_p_id].begin();it!=mcHotPool.HotPools[evict_p_id].end();++it){
//                    get<1>(EMEdgesIndex[it->ID1]) = -1;//hard delete
//                    get<2>(EMEdgesIndex[it->ID1]) = -1;
////                    get<0>(EMEdgesIndex[it->ID1]) = false;
//                }
                mcHotPool.HotPools[evict_p_id].clear();//clear
                mcHotPool.HotPools[evict_p_id].shrink_to_fit();//free memory
//                mcHotPool.im_partitions.erase(evict_p_id);
            }else if(evict_p_size > MCEDGE_PER_BLK){//if the remaining size is moderate, store vertices in temp file; if the remaining size is smaller than block size, remain the vertices in hop pool
                assert(mcHotPool.HotPools[evict_p_id].size()>0);
                /// write to stxxl vector
                int temp_id = -1;
                Timer tt;
                tt.start();

                for(auto it=mcHotPool.HotPools[evict_p_id].begin();it!=mcHotPool.HotPools[evict_p_id].end();++it){
                    if(!get<0>(EMEdgesIndex[it->ID1])){//if false, i.e. if the node needs to be pushed into external vector
                        if(it->ID1 != temp_id){
                            get<1>(EMEdgesIndex[it->ID1]) = mcHotPool.MCEdges_EM.size();
                            if(temp_id != -1){
                                get<2>(EMEdgesIndex[temp_id]) = mcHotPool.MCEdges_EM.size();
                            }
                            temp_id = it->ID1;
                        }
                        mcHotPool.MCEdges_EM.push_back(*it);
                    }
                }
                get<2>(EMEdgesIndex[temp_id]) = mcHotPool.MCEdges_EM.size();
                tt.stop();
                io_time += tt.GetRuntime() * 1000;

                get<0>(mcHotPool.clusterStatus[evict_p_id]) = 'S';
                mcHotPool.HotPools[evict_p_id].clear();//clear
                mcHotPool.HotPools[evict_p_id].shrink_to_fit();//free memory
//                mcHotPool.im_partitions.erase(evict_p_id);
            }
//            else {
//                ++partition_left;
//                if (partition_left > (partition_number/2)) {
//                    if (partition_left % 10 == 0) {
//                        cout << "The number of partitions which have less than 4KB memory consumption: " << partition_left << endl;
//                    }
//                }
//            }
        }
        
        //// read partition from original disk file to internal vector
        //check status of storage
        char filePath[200];
        int u, u_deg;
        int v, weight;
        MCEdgeT temp_edge;

        strcpy(filePath, DataPath);
        strcat(filePath, dataset.c_str());
        strcat(filePath, "/");
        strcat(filePath, partition_type.c_str());
        strcat(filePath, "/Partitions/");
        strcat(filePath, dataset.c_str());
        strcat(filePath, "_Partition_");
        strcat(filePath, to_string(PartitionSize).c_str());
        strcat(filePath, "_");
        strcat(filePath, to_string(target_p_id).c_str());
        strcat(filePath, ".bin");
        Timer tt;
        tt.start();
        FILE* file = fopen(filePath, "rb");
        assert(file != NULL);
        ReadBuffer read_io(file);
        bool flag_sizeUpdate = false;
        int temp_id = 0;

        if(get<1>(mcHotPool.clusterStatus[target_p_id]) == 0){//if this is the first read
            flag_sizeUpdate = true;
            //read partition
            while(!read_io.isend) {//if not read out
                read_io.read(&u);//get node id
                temp_edge.ID1 = u-1;
                if(u != temp_id){
                    get<1>(EMEdgesIndex[u-1]) = mcHotPool.HotPools[target_p_id].size();
                    if(temp_id != 0){
                        get<2>(EMEdgesIndex[temp_id-1]) = mcHotPool.HotPools[target_p_id].size();
                    }
                    temp_id = u;
                }

//                if(vertex_cri[u-1] == -1){//if not found
                assert(vertex_cri[u-1] == -1);
                    vertex_cri[u-1] = num_of_cri;
//                }
                read_io.read(&u_deg);//get node degree
                //vertex store in indexNode, edge store in e
                for (int i = 0; i < u_deg; ++i) {
                    read_io.read(&v);//read end node
                    temp_edge.ID2 = v-1;
                    for (int j = 0; j < NUM_OF_CRITERIA; ++j) {
                        read_io.read(&weight);
                        if(j < num_of_cri){
                            temp_edge.putW(j,weight);
                        }
                    }
                    mcHotPool.HotPools[target_p_id].emplace_back(temp_edge);

                }
                get<2>(mcHotPool.clusterStatus[target_p_id]) += u_deg;

//                if(flag_sizeUpdate)
                    get<1>(mcHotPool.clusterStatus[target_p_id]) += u_deg;
            }
            io_num += read_io.io_number;//update #IO
        }else{//if this is not the first read
            //read partition
            while(!read_io.isend) {//if not read out
                read_io.read(&u);//get node id
                temp_edge.ID1 = u-1;
                if(u != temp_id){
                    get<1>(EMEdgesIndex[u-1]) = mcHotPool.HotPools[target_p_id].size();
                    if(temp_id != 0){
                        get<2>(EMEdgesIndex[temp_id-1]) = mcHotPool.HotPools[target_p_id].size();
                    }
                    temp_id = u;
                }

                if(vertex_cri[u-1] == -1){//if not found
                    vertex_cri[u-1] = num_of_cri;
                }
                read_io.read(&u_deg);//get node degree
                //vertex store in indexNode, edge store in e
                for (int i = 0; i < u_deg; ++i) {
                    read_io.read(&v);//read end node
                    temp_edge.ID2 = v-1;
                    for (int j = 0; j < NUM_OF_CRITERIA; ++j) {
                        read_io.read(&weight);
                        if(j < num_of_cri){
                            temp_edge.putW(j,weight);
                        }
                    }
                    mcHotPool.HotPools[target_p_id].emplace_back(temp_edge);

                }
                if (!get<0>(EMEdgesIndex[u - 1]))//if false, if it has not been evicted
                    get<2>(mcHotPool.clusterStatus[target_p_id]) += u_deg;
            }
            io_num += 1;//update #IO
        }


        get<2>(EMEdgesIndex[u-1]) = mcHotPool.HotPools[target_p_id].size();
//        cout<<mcHotPool.HotPools[target_p_id].size()<<endl;
//        mcHotPool.im_partitions.insert(target_p_id);
        //io_num += read_io.io_number;//update #IO
        fclose(file);
        tt.stop();
        io_time += tt.GetRuntime() * 1000;
        cluster_to_IO[target_p_id] = read_io.io_number;

        get<0>(mcHotPool.clusterStatus[target_p_id]) = 'I';
    }
    //Function of multi-stop algorithm with io optimization
    void EMMCGraph::EM_MC_MultiHop(int node_start, int node_end){
        if (node_start == node_end) {
            cout << "Source vertex and Target vertex are the same!!!" << endl;
            return;
        }
        //statstics of IO
        Stats stats_b, stats_e;
        stats_b.start();
        //variables
        vector<vector<Distance>> mc_cost(num_of_cri, vector<Distance>(node_num, INF));   //vector of cost
        //vector<vector<NodeId>> mc_pre(num_of_cri, vector<int>(node_num, -1));      //vector of predecessor id
        vector<vector<bool>> mc_closed(num_of_cri, vector<bool>(node_num));    //flag of vertex closing
        vector<int> cluster_to_cri(partition_number,0);//used to map which criterion read this partition
        set<int> set_remain;//set for storing the id of unshared criteria
        NodeId item_id;
        bool flag_finish = false;
        bool flag_mcFinished = false;
        int partition_id;
        int cri_i;
        NodeId temp_id;
        Distance temp_dis;
        Distance item_dis;
        int temp_number;
        int temp_degree = 0;
        benchmark::heap<2, int, int> cri_pqueue(num_of_cri);
        vector<int> cri_hops(num_of_cri, 0);
        vector<int> cri_to_cluster(num_of_cri, 0);
        int cri_empty = 0;
        partition_left = 0;
        //hot pool
        MCHotPool<VectorMCEdgesEMTuple_IO> mc_HotPool(Partition_N);       //Hot pool for one-pass algorithms
        mc_HotPool.init(partition_number);
        mc_HotPool.set_capacity(Partition_N);

        Timer tt;
        tt.start();

        //Initiation of start node
        for (int i = 0; i < num_of_cri; ++i) {
            set_cri.insert(i);
            mc_cost[i][node_start] = 0;//cost of start node
            EM_MC_PQueue[i]->push(VertexCost(node_start,0));
            cri_to_cluster[i] = node_to_cluster[node_start];
            cri_pqueue.update(i, 0);
        }
        //first read
        //GraphMCReadCluster_New(mc_HotPool, node_start, node_to_cluster[node_start], MultiHops);//if in disk, read corresponding partition
        //++p_num;
        cluster_to_cri[node_to_cluster[node_start]] = 0;
        //Iteration
        while (!EM_JudgeEmpty(EM_MC_PQueue)) {//for every node in pqueue !mc_finished[0] || !mc_finished[1] || !mc_finished[2]
            set_remain = set_cri;
            /// read by synchronizing the hops
//                temp_cluster_id = cri_to_cluster[cri_pqueue.top_value()];
            GraphMCReadCluster_New(mc_HotPool, item_id, cri_to_cluster[cri_pqueue.top_value()], MultiHops);//if in disk, read corresponding partition
            ++p_num;
            while (!set_remain.empty()) {
                cri_i = *set_remain.begin();

                while (mc_closed[cri_i][EM_MC_PQueue[cri_i]->top().id]) {//if found in visited set  && !EM_MC_PQueue[cri_i]->empty()
                    EM_MC_PQueue[cri_i]->pop();
                    if(EM_MC_PQueue[cri_i]->empty()){
                        mc_finished[cri_i] = true;
                        set_cri.erase(cri_i);
                        mc_min_cost[cri_i] = INF;
                        cri_hops[cri_i] = INF;
                        cri_pqueue.update(cri_i, INF);
                        set_remain.erase(cri_i);
                        ++cri_empty;
                        goto jump1;
                    }
                }
                item_id = EM_MC_PQueue[cri_i]->top().id;// top min item
                item_dis = EM_MC_PQueue[cri_i]->top().cost;// top min item

                partition_id = node_to_cluster[item_id];

                while (get<0>(mc_HotPool.clusterStatus[partition_id]) == 'I' || get<0>(mc_HotPool.clusterStatus[partition_id]) == 'S') {//if found
                    if (item_id == node_end) {//if reach target node, this criterion ends.
    //                mc_pqueue[cri_i].clear();
                        mc_finished[cri_i] = true;
                        set_cri.erase(cri_i);
                        mc_min_cost[cri_i] = mc_cost[cri_i][item_id];
//                        cout<<"Minimal cost of criterion " << cri_i << ": "<<mc_min_cost[cri_i]<<endl;
                        flag_finish = true;
                        cri_hops[cri_i] = INF;
                        cri_pqueue.update(cri_i, INF);
                        ++visited_number;
                        flag_mcFinished = true;
                        for (auto it = mc_finished.begin(); it != mc_finished.end(); ++it) {
                            if (!*it) {//if none of them is false, then flag_mcFinished = true
                                flag_mcFinished = false;
                                break;
                            }
                        }
                        break;
                    }
                    else {
                        ++cri_hops[cri_i];
                        //slack
                        /// read from mc_HotPool
                        partition_id = node_to_cluster[item_id];
                        if (cluster_to_cri[partition_id] != cri_i) {//if the cluster is not read by this criterion
                            ++share_number;
                        }
                        //set closed
                        mc_closed[cri_i][item_id] = true;
                        EM_MC_PQueue[cri_i]->pop();//pop min item
                        ++visited_number;//update #visited_number
                        switch (get<0>(mc_HotPool.clusterStatus[partition_id])) {// deal according to partition storage status
                        case 'I': {//if in memory
                            temp_degree = 0;
                            for (int i = get<1>(EMEdgesIndex[item_id]); i < get<2>(EMEdgesIndex[item_id]); ++i) {
                                assert(mc_HotPool.HotPools[partition_id][i].ID1 == item_id);
                                ++temp_degree;
                                temp_id = mc_HotPool.HotPools[partition_id][i].ID2;
                                if (mc_HotPool.HotPools[partition_id][i].getW(cri_i)==INF100)
                                    continue;
                                if (mc_closed[cri_i][temp_id])//if closed
                                    continue;
                                temp_dis = item_dis + mc_HotPool.HotPools[partition_id][i].getW(cri_i);
                                if (mc_cost[cri_i][temp_id] > temp_dis) {//slack operation
                                    mc_cost[cri_i][temp_id] = temp_dis;
                                    //mc_pre[cri_i][temp_id] = item_id;
                                    EM_MC_PQueue[cri_i]->push(VertexCost(temp_id, temp_dis));
                                }
                            }
                            break;
                        }
                        case 'S': {//if in external vector
                            assert(get<2>(EMEdgesIndex[item_id]) > 0);
                            temp_degree = 0;
                            for (int i = get<1>(EMEdgesIndex[item_id]); i < get<2>(EMEdgesIndex[item_id]); ++i) {
                                assert(mc_HotPool.MCEdges_EM[i].ID1 == item_id);
                                ++temp_degree;
                                temp_id = mc_HotPool.MCEdges_EM[i].ID2;
                                if (mc_HotPool.MCEdges_EM[i].getW(cri_i)==INF100)
                                    continue;
                                if (mc_closed[cri_i][temp_id])//if closed
                                    continue;
                                temp_dis = item_dis + mc_HotPool.MCEdges_EM[i].getW(cri_i);
                                if (mc_cost[cri_i][temp_id] > temp_dis) {//slack operation
                                    mc_cost[cri_i][temp_id] = temp_dis;
                                    //mc_pre[cri_i][temp_id] = item_id;
                                    EM_MC_PQueue[cri_i]->push(VertexCost(temp_id, temp_dis));
                                }
                            }
                            break;
                        }
                        default: {
                            cout << "Wrong partition status!!!" << endl;
                            break;
                        }
                        }
                        --vertex_cri[item_id];
                        assert(vertex_cri[item_id] >= 0);
                        //// evict vertex
                        if (vertex_cri[item_id] - cri_empty <= 0) {//evict valid in-memory vertices immediately
                            if (get<0>(mc_HotPool.clusterStatus[partition_id]) == 'I') {
                                get<0>(EMEdgesIndex[item_id]) = true;
                                get<2>(mc_HotPool.clusterStatus[partition_id]) -= temp_degree;
                                if (get<2>(mc_HotPool.clusterStatus[partition_id]) <= 0) {//if empty, clear it immediately
                                    assert(get<2>(mc_HotPool.clusterStatus[partition_id]) == 0);
                                    for(int i=partitions_info[partition_id];i<partitions_info[partition_id+1];++i){
                                        get<1>(EMEdgesIndex[i]) = -1;//hard delete
                                        get<2>(EMEdgesIndex[i]) = -1;
                                    }
//                                    for(auto it=mc_HotPool.HotPools[partition_id].begin();it!=mc_HotPool.HotPools[partition_id].end();++it){
//                                        get<1>(EMEdgesIndex[it->ID1]) = -1;//hard delete
//                                        get<2>(EMEdgesIndex[it->ID1]) = -1;
//                                    }
                                    mc_HotPool.HotPools[partition_id].clear();
                                    mc_HotPool.HotPools[partition_id].shrink_to_fit();//free memory
//                                    mc_HotPool.im_partitions.erase(partition_id);
                                    mc_HotPool.LRU_IMCluster.erase(partition_id);
                                    get<0>(mc_HotPool.clusterStatus[partition_id]) = 'D';
                                }
                            }
                        }

                        while (mc_closed[cri_i][EM_MC_PQueue[cri_i]->top().id] && !EM_MC_PQueue[cri_i]->empty()) {//if already closed
                            EM_MC_PQueue[cri_i]->pop();
                        }
                        if(EM_MC_PQueue[cri_i]->empty()){
                            mc_finished[cri_i] = true;
                            set_cri.erase(cri_i);
                            mc_min_cost[cri_i] = INF;
                            cri_hops[cri_i] = INF;
                            cri_pqueue.update(cri_i, INF);
                            set_remain.erase(cri_i);
                            ++cri_empty;
                            goto jump1;
                        }
                        item_id = EM_MC_PQueue[cri_i]->top().id;// top min item
                        item_dis = EM_MC_PQueue[cri_i]->top().cost;// top min item

                        partition_id = node_to_cluster[item_id];
                    }
                }
                set_remain.erase(cri_i);
                if (!flag_finish) {
                    cluster_to_cri[node_to_cluster[item_id]] = cri_i;
                    cri_to_cluster[cri_i] = node_to_cluster[item_id];
                    cri_pqueue.update(cri_i, cri_hops[cri_i]);
                }else{
                    flag_finish = false;
                }
    //            if (set_remain.empty() && !flag_mcFinished) {//if all criteria are processed &&!EM_JudgeEmpty(em_mc_pqueue)
    //                /// read by synchronizing the hops
    ////                temp_cluster_id = cri_to_cluster[cri_pqueue.top_value()];
    //                GraphMCReadCluster_New(mc_HotPool, item_id, cri_to_cluster[cri_pqueue.top_value()],MultiHops);//if in disk, read corresponding partition
    //                ++p_num;
    //            }
                jump1: "there is a jump.";
            }
        }
        tt.stop();
        query_time += tt.GetRuntime() * 1000;
//        cout<<"Query time: "<<tt.GetRuntime()<<" s."<<endl;
        //io record
        stats_e = stats_b.get_stats();
        if (stats_e.m_noOfWrites > 0 || stats_e.m_noOfReads > 0) {
            cout << "There are " << stats_e.m_noOfWrites + stats_e.m_noOfReads << " io for stxxl." << endl;
            //cout<<stats_e<<endl;
            EM_IORecord(stats_e);
        }
        mc_HotPool.clear();
    }
    //Function of Bidirectional Multi-Hop algorithm, version 1
    void EMMCGraph::EM_MC_BiMultiHop(int node_start, int node_end) {
        if (node_start == node_end) {
            cout << "Source vertex and Target vertex are the same!!!" << endl;
            return;
        }
        //statstics of IO
        Stats stats_b, stats_e;
        stats_b.start();
        //variables
        set<int> set_remain;//set for storing the id of unshared criteria
        vector<vector<Distance>> mc_cost(num_of_cri, vector<Distance>(node_num, INF));   //vector of cost
        //vector<vector<NodeId>> mc_pre(num_of_cri, vector<int>(node_num, -1));      //vector of predecessor id
        vector<vector<bool>> mc_closed(num_of_cri, vector<bool>(node_num));    //flag of vertex closing
        vector<vector<Distance>> mc_cost_r(num_of_cri, vector<Distance>(node_num, INF));   //vector of cost
        //vector<vector<NodeId>> mc_pre_r(num_of_cri, vector<int>(node_num, -1));      //vector of predecessor id
        vector<vector<bool>> mc_closed_r(num_of_cri, vector<bool>(node_num));    //flag of vertex closing
        vector<int> cluster_to_cri(partition_number, 0);//used to map which criterion read this partition
        vector<int8_t> cluster_to_bi(partition_number, 0);   //used to indicate the partition is read by forward search or reverse search
        vector<NodeId> mc_terminate_id(num_of_cri, -1);;//termination id of bi-dijkstra
        NodeId item_id, item_id_r;
        NodeId temp_id, temp_id_r;
        Distance temp_dis, temp_dis_r;
        Distance item_dis, item_dis_r;
        int partition_id, partition_id_r;
        int cri_i;
        int temp_cluster_id, temp_cluster_id_r;
        bool temp_bool;
        int temp_degree;
        bool flag_mcFinished = false;
        benchmark::heap<2, int, int> cri_pqueue(num_of_cri);
        vector<int> cri_hops(num_of_cri, 0);
        vector<int> cri_to_cluster(num_of_cri, 0);
        benchmark::heap<2, int, int> cri_pqueue_r(num_of_cri);
        vector<int> cri_hops_r(num_of_cri, 0);
        vector<int> cri_to_cluster_r(num_of_cri, 0);
        vector<bool> clusterRead(partition_number, false);
        int cri_empty = 0;
        partition_left = 0;
        //hot pool
        MCHotPool<VectorMCEdgesEMTuple_IO_Bi> mc_HotPool(Partition_N_Bi);       //Hot pool for one-pass algorithms
        mc_HotPool.init(partition_number);
        mc_HotPool.set_capacity(Partition_N_Bi);

        Timer tt;
        tt.start();

        //Initiation of start node
        for (int i = 0; i < num_of_cri; i++) {
            set_cri.insert(i);
            mc_cost[i][node_start] = 0;//cost of start node
            EM_MC_PQueue[i]->push(VertexCost(node_start, 0));
            mc_cost_r[i][node_end] = 0;
            EM_MC_PQueue_r[i]->push(VertexCost(node_end, 0));
            cri_to_cluster[i] = node_to_cluster[node_start];
            cri_to_cluster_r[i] = node_to_cluster[node_end];
            cri_pqueue.update(i,0);
            cri_pqueue_r.update(i,0);
        }

        //first read
        //GraphMCReadCluster_New_Bi(mc_HotPool, node_start, node_to_cluster[node_start], cluster_to_bi, FORWARD, clusterRead);//if in disk, read corresponding partition
        //++p_num;
        //cluster_to_cri[node_to_cluster[node_start]] = 0;
        //cluster_to_bi[node_to_cluster[node_start]] = FORWARD;
        //if (node_to_cluster[node_end] != node_to_cluster[node_start]) {
        //    //            GraphMCReadCluster_New_Bi(mc_HotPool, node_end, node_to_cluster[node_end],cluster_to_bi, false);//if in disk, read corresponding partition
        //    GraphMCReadCluster_New_Bi(mc_HotPool, node_end, node_to_cluster[node_end], cluster_to_bi, REVERSE, clusterRead);//if in disk, read corresponding partition
        //    ++p_num;
        //    cluster_to_cri[node_to_cluster[node_end]] = 0;
        //    cluster_to_bi[node_to_cluster[node_end]] = REVERSE;
        //}
        cluster_to_cri[node_to_cluster[node_start]] = 0;
        cluster_to_cri[node_to_cluster[node_end]] = 0;
        //Iteration
        while (!EM_JudgeEmptyBi(EM_MC_PQueue, EM_MC_PQueue_r)) {//for every node in pqueue !mc_finished[0] || !mc_finished[1] || !mc_finished[2]
            set_remain = set_cri;
            /// read by synchronizing the hops
            temp_cluster_id = cri_to_cluster[cri_pqueue.top_value()];
            temp_cluster_id_r = cri_to_cluster_r[cri_pqueue_r.top_value()];
            if(cri_hops[cri_pqueue.top_value()] <= cri_hops_r[cri_pqueue_r.top_value()]){
                GraphMCReadCluster_New_Bi(mc_HotPool, temp_cluster_id, cluster_to_bi, FORWARD, clusterRead);//if in disk, read corresponding partition item_id,
                cluster_to_bi[temp_cluster_id] = FORWARD;
                ++p_num;
            }else{
                GraphMCReadCluster_New_Bi(mc_HotPool, temp_cluster_id_r, cluster_to_bi, REVERSE, clusterRead);//if in disk, read corresponding partition item_id_r,
                cluster_to_bi[temp_cluster_id_r] = REVERSE;
                ++p_num;
            }

            while (!set_remain.empty()) {
                cri_i = *set_remain.begin();//pick one criterion to process
                //deal with invalid top elements
                while (mc_closed[cri_i][EM_MC_PQueue[cri_i]->top().id]) {
                    EM_MC_PQueue[cri_i]->pop();
                    if(EM_MC_PQueue[cri_i]->empty())
                    {
                        mc_finished[cri_i] = true;
                        set_cri.erase(cri_i);
                        cri_hops[cri_i] = INF;
                        cri_pqueue.update(cri_i, INF);
                        cri_hops_r[cri_i] = INF;
                        cri_pqueue_r.update(cri_i, INF);
                        ++cri_empty;
                        goto jump1;
                    }
                }
                while (mc_closed_r[cri_i][EM_MC_PQueue_r[cri_i]->top().id]) {
                    EM_MC_PQueue_r[cri_i]->pop();
                    if(EM_MC_PQueue_r[cri_i]->empty())
                    {
                        mc_finished[cri_i] = true;
                        set_cri.erase(cri_i);
                        cri_hops[cri_i] = INF;
                        cri_pqueue.update(cri_i, INF);
                        cri_hops_r[cri_i] = INF;
                        cri_pqueue_r.update(cri_i, INF);
                        ++cri_empty;
                        goto jump1;
                    }
                }

                //// Forward searching
                item_id = EM_MC_PQueue[cri_i]->top().id;// top min item
                item_dis = EM_MC_PQueue[cri_i]->top().cost;// top min item

                partition_id = node_to_cluster[item_id];
                while (get<2>(EMEdgesIndex_Bi[item_id]) != -1  && !mc_closed_r[cri_i][item_id]) {//if the vertex is in memory or stxxl vector.
                    if(mc_closed_r[cri_i][item_id] || item_id == node_end){//// Termination judging  || item_id == node_end
                        mc_finished[cri_i] = true;
                        set_cri.erase(cri_i);
                        cri_hops[cri_i] = INF;
                        cri_pqueue.update(cri_i, INF);
                        cri_hops_r[cri_i] = INF;
                        cri_pqueue_r.update(cri_i, INF);
                        if(item_id == node_end){
                            cout<<"Unidirectional search!"<<endl;
                            mc_min_cost[cri_i] = mc_cost[cri_i][item_id];
                        }
//                            cout<<node_start<<" "<<node_end<<endl;
//                        cout<<"T2 Minimal cost of criterion "<<cri_i<<": "<<mc_min_cost[cri_i]<<endl;
                        set_remain.erase(cri_i);
                        goto jump1;
                    }
                    ++cri_hops[cri_i];

                    // relaxation
                    if (cluster_to_cri[node_to_cluster[item_id]] != cri_i) {//if the cluster is not read by this criterion
                        ++share_number;
                    }
                    //set closed
                    mc_closed[cri_i][item_id] = true;
                    EM_MC_PQueue[cri_i]->pop();//pop min item
                    ++visited_number;//update #visited_number
                    //slack
                    switch (get<0>(mc_HotPool.clusterStatus[partition_id])) {// deal according to partition storage status
                    case 'I': {//if in memory
                        temp_degree = 0;
                        for (int i = get<2>(EMEdgesIndex_Bi[item_id]); i < get<3>(EMEdgesIndex_Bi[item_id]); ++i) {
                            assert(mc_HotPool.HotPools[partition_id][i].ID1 == item_id);
                            ++temp_degree;
                            temp_id = mc_HotPool.HotPools[partition_id][i].ID2;
                            if (mc_HotPool.HotPools[partition_id][i].getW(cri_i)==INF100)
                                continue;
                            temp_dis = item_dis + mc_HotPool.HotPools[partition_id][i].getW(cri_i);
                            if (!mc_closed[cri_i][temp_id]) {//if not closed
                                if (mc_cost[cri_i][temp_id] > temp_dis) {//slack operation
                                    mc_cost[cri_i][temp_id] = temp_dis;
                                    //mc_pre[cri_i][temp_id] = item_id;
                                    EM_MC_PQueue[cri_i]->push(VertexCost(temp_id, temp_dis));
                                }
                            }
                            //update min cost
                            temp_bool = mc_closed_r[cri_i][temp_id];
                            if (temp_bool && temp_dis + mc_cost_r[cri_i][temp_id] < mc_min_cost[cri_i]) {
                                mc_min_cost[cri_i] = temp_dis + mc_cost_r[cri_i][temp_id];
                                mc_terminate_id[cri_i] = temp_id;
                            }
                        }
                        break;
                    }
                    case 'S': {//if in external vector
                        assert(get<3>(EMEdgesIndex_Bi[item_id]) > 0);
                        temp_degree = 0;
                        for (int i = get<2>(EMEdgesIndex_Bi[item_id]); i < get<3>(EMEdgesIndex_Bi[item_id]); ++i) {
                            assert(mc_HotPool.MCEdges_EM[i].ID1 == item_id);
                            ++temp_degree;
                            temp_id = mc_HotPool.MCEdges_EM[i].ID2;
                            if (mc_HotPool.MCEdges_EM[i].getW(cri_i)==INF100)
                                continue;
                            temp_dis = item_dis + mc_HotPool.MCEdges_EM[i].getW(cri_i);
                            if (!mc_closed[cri_i][temp_id]) {//if not closed
                                if (mc_cost[cri_i][temp_id] > temp_dis) {//slack operation
                                    mc_cost[cri_i][temp_id] = temp_dis;
                                    //mc_pre[cri_i][temp_id] = item_id;
                                    EM_MC_PQueue[cri_i]->push(VertexCost(temp_id, temp_dis));
                                }
                            }
                            //update min cost
                            temp_bool = mc_closed_r[cri_i][temp_id];
                            if (temp_bool && temp_dis + mc_cost_r[cri_i][temp_id] < mc_min_cost[cri_i]) {
                                mc_min_cost[cri_i] = temp_dis + mc_cost_r[cri_i][temp_id];
                                mc_terminate_id[cri_i] = temp_id;
                            }
                        }
                        break;
                    }
                    default: {
                        cout << "Wrong vertex status!!!" << endl; break;
                    }
                    }
                    if (cluster_to_bi[partition_id] == FORWARD) {//if the partition is read by forward search, reduce the counter; otherwise, it means the partition is shared by both forward search and reverse search
                        --vertex_cri_Bi[item_id].first;
                        assert(vertex_cri_Bi[item_id].first >= 0);
                        //// evict vertex
                        if (vertex_cri_Bi[item_id].first - cri_empty <= 0) {//set the flag of evict to true
//                                assert(get<0>(mc_HotPool.clusterStatus[partition_id]) == 'I');
                            if (get<0>(mc_HotPool.clusterStatus[partition_id]) == 'I') {
                                get<0>(EMEdgesIndex_Bi[item_id]) = true;//soft delete
//                                    get<2>(EMEdgesIndex_Bi[item_id]) = -1; get<3>(EMEdgesIndex_Bi[item_id]) = -1;
                                get<2>(mc_HotPool.clusterStatus[partition_id]) -= temp_degree;
                                if (get<2>(mc_HotPool.clusterStatus[partition_id]) <= 0) {//if empty, evict the partition immediately
//                                    cout << "Toss immediate: Forward." << endl;
                                    assert(get<2>(mc_HotPool.clusterStatus[partition_id]) == 0);
                                    for(int i=partitions_info[partition_id];i<partitions_info[partition_id+1];++i){
                                        get<2>(EMEdgesIndex_Bi[i]) = -1;//hard delete
                                        get<3>(EMEdgesIndex_Bi[i]) = -1;
                                    }
//                                        for (auto it = mc_HotPool.HotPools[partition_id].begin(); it != mc_HotPool.HotPools[partition_id].end(); ++it) {
//                                            get<2>(EMEdgesIndex_Bi[it->ID1]) = -1;//hard delete
//                                            get<3>(EMEdgesIndex_Bi[it->ID1]) = -1;
//                                        }
                                    mc_HotPool.HotPools[partition_id].clear();
                                    mc_HotPool.HotPools[partition_id].shrink_to_fit();//free memory
//                                        mc_HotPool.im_partitions.erase(partition_id);
                                    mc_HotPool.LRU_IMCluster.erase(partition_id);
                                    get<0>(mc_HotPool.clusterStatus[partition_id]) = 'D';
                                }
                            }
                        }
                    }

                    while (mc_closed[cri_i][EM_MC_PQueue[cri_i]->top().id] && !EM_MC_PQueue[cri_i]->empty()) {//if already closed
                        EM_MC_PQueue[cri_i]->pop();
                    }
                    if(EM_MC_PQueue[cri_i]->empty()){
                        mc_finished[cri_i] = true;
                        set_cri.erase(cri_i);
                        mc_min_cost[cri_i] = INF;
                        cri_hops[cri_i] = INF;
                        cri_pqueue.update(cri_i, INF);
                        cri_hops_r[cri_i] = INF;
                        cri_pqueue_r.update(cri_i, INF);
                        set_remain.erase(cri_i);
                        ++cri_empty;
                        goto jump1;
                    }
                    item_id = EM_MC_PQueue[cri_i]->top().id;// top min item
                    item_dis = EM_MC_PQueue[cri_i]->top().cost;// top min item
//                        if(item_id == 22683035)
//                            cout<<item_id<<endl;
                    partition_id = node_to_cluster[item_id];
                }
                //if(mc_finished[cri_i])
                //    break;
                //cluster_to_cri[node_to_cluster[item_id]] = cri_i;
                //cri_to_cluster[cri_i] = node_to_cluster[item_id];
                //cri_pqueue.update(cri_i, cri_hops[cri_i]);
                //// Reverse searching
                item_id_r = EM_MC_PQueue_r[cri_i]->top().id;// top min item
                item_dis_r = EM_MC_PQueue_r[cri_i]->top().cost;// top min item

                partition_id_r = node_to_cluster[item_id_r];

                while (get<2>(EMEdgesIndex_Bi[item_id_r]) != -1 && !mc_closed[cri_i][item_id_r]) {//if the vertex is in memory or stxxl vector.
                    if(mc_closed[cri_i][item_id_r] || item_id_r == node_start){//// Termination judging
                        mc_finished[cri_i] = true;
                        set_cri.erase(cri_i);
                        cri_hops[cri_i] = INF;
                        cri_pqueue.update(cri_i, INF);
                        cri_hops_r[cri_i] = INF;
                        cri_pqueue_r.update(cri_i, INF);
                        if(item_id_r == node_start){
                            cout<<"Unidirectional search!"<<endl;
                            mc_min_cost[cri_i] = mc_cost_r[cri_i][item_id_r];
                        }
//                    cout<<node_start<<" "<<node_end<<endl;
//                        cout<<"T2 Minimal cost of criterion "<<cri_i<<": "<<mc_min_cost[cri_i]<<endl;
                        set_remain.erase(cri_i);
                        goto jump1;
                    }

                    ++cri_hops_r[cri_i];
                    // relaxation
                    if (cluster_to_cri[node_to_cluster[item_id_r]] != cri_i) {//if the cluster is not read by this criterion
                        ++share_number;
                    }
                    //set closed
                    mc_closed_r[cri_i][item_id_r] = true;
                    EM_MC_PQueue_r[cri_i]->pop();//pop min item
                    ++visited_number;//update #visited_number
                    //slack
                    switch (get<0>(mc_HotPool.clusterStatus[partition_id_r])) {// deal according to partition storage status
                    case 'I': {//if in memory
                        temp_degree = 0;
                        for (int i = get<2>(EMEdgesIndex_Bi[item_id_r]); i < get<3>(EMEdgesIndex_Bi[item_id_r]); ++i) {
                            assert(mc_HotPool.HotPools[partition_id_r][i].ID1 == item_id_r);
                            ++temp_degree;
                            temp_id_r = mc_HotPool.HotPools[partition_id_r][i].ID2;
                            if (mc_HotPool.HotPools[partition_id_r][i].getW(cri_i)==INF100)
                                continue;
                            temp_dis_r = item_dis_r + mc_HotPool.HotPools[partition_id_r][i].getW(cri_i);
                            if (!mc_closed_r[cri_i][temp_id_r]) {//if not closed
                                if (mc_cost_r[cri_i][temp_id_r] > temp_dis_r) {//slack operation
                                    mc_cost_r[cri_i][temp_id_r] = temp_dis_r;
                                    //mc_pre_r[cri_i][temp_id_r] = item_id_r;
                                    EM_MC_PQueue_r[cri_i]->push(VertexCost(temp_id_r, temp_dis_r));
                                }
                            }
                            //update min cost
                            temp_bool = mc_closed[cri_i][temp_id_r];
                            if (temp_bool && temp_dis_r + mc_cost[cri_i][temp_id_r] < mc_min_cost[cri_i]) {
                                mc_min_cost[cri_i] = temp_dis_r + mc_cost[cri_i][temp_id_r];
                                mc_terminate_id[cri_i] = temp_id_r;
                            }
                        }
                        break;
                    }
                    case 'S': {//if in external vector
                        assert(get<3>(EMEdgesIndex_Bi[item_id_r]) > 0);
                        temp_degree = 0;
                        for (int i = get<2>(EMEdgesIndex_Bi[item_id_r]); i < get<3>(EMEdgesIndex_Bi[item_id_r]); ++i) {
//                                if (mc_HotPool.MCEdges_EM[i].ID1 != item_id_r) {
//                                    cout << item_id_r << " " << mc_HotPool.MCEdges_EM[i].ID1 << endl;
//                                }
                            assert(mc_HotPool.MCEdges_EM[i].ID1 == item_id_r);
                            ++temp_degree;
                            temp_id_r = mc_HotPool.MCEdges_EM[i].ID2;
                            if (mc_HotPool.MCEdges_EM[i].getW(cri_i)==INF100)
                                continue;
                            temp_dis_r = item_dis_r + mc_HotPool.MCEdges_EM[i].getW(cri_i);
                            if (!mc_closed_r[cri_i][temp_id_r]) {//if not closed
                                if (mc_cost_r[cri_i][temp_id_r] > temp_dis_r) {//slack operation
                                    mc_cost_r[cri_i][temp_id_r] = temp_dis_r;
                                    //mc_pre_r[cri_i][temp_id_r] = item_id_r;
                                    EM_MC_PQueue_r[cri_i]->push(VertexCost(temp_id_r, temp_dis_r));
                                }
                            }
                            //update min cost
                            temp_bool = mc_closed[cri_i][temp_id_r];
                            if (temp_bool && temp_dis_r + mc_cost[cri_i][temp_id_r] < mc_min_cost[cri_i]) {
                                mc_min_cost[cri_i] = temp_dis_r + mc_cost[cri_i][temp_id_r];
                                mc_terminate_id[cri_i] = temp_id_r;
                            }
                        }
                        break;
                    }
                    default: {
                        cout << "Wrong vertex status!!!" << endl; break;
                    }
                    }
                    if (cluster_to_bi[partition_id_r] == REVERSE) {//if the partition is read by reverse search, reduce the counter; otherwise, if means the partition is shared by both forward search and reverse search
                        --vertex_cri_Bi[item_id_r].second;
                        assert(vertex_cri_Bi[item_id_r].second >= 0);
                        //// evict vertex
                        if (vertex_cri_Bi[item_id_r].second - cri_empty <= 0) {//evict valid in-memory vertices immediately
//                                assert(get<0>(mc_HotPool.clusterStatus[partition_id_r]) == 'I');
                            if (get<0>(mc_HotPool.clusterStatus[partition_id_r]) == 'I') {
                                get<1>(EMEdgesIndex_Bi[item_id_r]) = true;//soft delete
//                                    get<2>(EMEdgesIndex_Bi[item_id_r]) = -1; get<3>(EMEdgesIndex_Bi[item_id_r]) = -1;
                                get<2>(mc_HotPool.clusterStatus[partition_id_r]) -= temp_degree;
                                if (get<2>(mc_HotPool.clusterStatus[partition_id_r]) <= 0) {//if empty, clear it immediately
//                                    cout << "Toss immediate: Reverse." << endl;
                                    assert(get<2>(mc_HotPool.clusterStatus[partition_id_r]) == 0);
                                    for(int i=partitions_info[partition_id_r];i<partitions_info[partition_id_r+1];++i){
                                        get<2>(EMEdgesIndex_Bi[i]) = -1;//hard delete
                                        get<3>(EMEdgesIndex_Bi[i]) = -1;
                                    }
//                                        for (auto it = mc_HotPool.HotPools[partition_id_r].begin(); it != mc_HotPool.HotPools[partition_id_r].end(); ++it) {
//                                            get<2>(EMEdgesIndex_Bi[it->ID1]) = -1;//hard delete
//                                            get<3>(EMEdgesIndex_Bi[it->ID1]) = -1;
//                                        }
                                    mc_HotPool.HotPools[partition_id_r].clear();
                                    mc_HotPool.HotPools[partition_id_r].shrink_to_fit();//free memory
//                                        mc_HotPool.im_partitions.erase(partition_id_r);
                                    mc_HotPool.LRU_IMCluster.erase(partition_id_r);
                                    get<0>(mc_HotPool.clusterStatus[partition_id_r]) = 'D';
                                }
                            }
                        }
                    }

                    while (mc_closed_r[cri_i][EM_MC_PQueue_r[cri_i]->top().id] && !EM_MC_PQueue_r[cri_i]->empty()) {
                        EM_MC_PQueue_r[cri_i]->pop();
                    }
                    if(EM_MC_PQueue_r[cri_i]->empty()){
                        mc_finished[cri_i] = true;
                        set_cri.erase(cri_i);
                        mc_min_cost[cri_i] = INF;
                        cri_hops[cri_i] = INF;
                        cri_pqueue.update(cri_i, INF);
                        cri_hops_r[cri_i] = INF;
                        cri_pqueue_r.update(cri_i, INF);
                        set_remain.erase(cri_i);
                        ++cri_empty;
                        goto jump1;
                    }
                    item_id_r = EM_MC_PQueue_r[cri_i]->top().id;// top min item
                    item_dis_r = EM_MC_PQueue_r[cri_i]->top().cost;// top min item

                    partition_id_r = node_to_cluster[item_id_r];
                }
                //cri_to_cluster_r[cri_i] = node_to_cluster[item_id_r];
                //cluster_to_cri[cri_to_cluster_r[cri_i]] = cri_i;//node_to_cluster[item_id_r]
                //cri_pqueue_r.update(cri_i, cri_hops_r[cri_i]);
                //// Termination Judge
                if (EM_MC_PQueue[cri_i]->top().cost + EM_MC_PQueue_r[cri_i]->top().cost >= mc_min_cost[cri_i]) {//condition of termination
                    mc_finished[cri_i] = true;
                    set_cri.erase(cri_i);
                    cri_hops[cri_i] = INF;
                    cri_pqueue.update(cri_i, INF);
                    cri_hops_r[cri_i] = INF;
                    cri_pqueue_r.update(cri_i, INF);
//                    cout<<node_start<<" "<<node_end<<endl;
//                    cout<<"Minimal cost of criterion "<<cri_i<<": "<<mc_min_cost[cri_i]<<endl;
                    flag_mcFinished = true;
                    for (auto it = mc_finished.begin(); it != mc_finished.end(); ++it) {
                        if (!*it) {//if none of them is false, then flag_mcFinished = true
                            flag_mcFinished = false;
                            break;
                        }
                    }
                }
                set_remain.erase(cri_i);
                if (!mc_finished[cri_i]) {
                    cluster_to_cri[node_to_cluster[item_id]] = cri_i;
                    cri_to_cluster[cri_i] = node_to_cluster[item_id];
                    cri_pqueue.update(cri_i, cri_hops[cri_i]);
                    cri_to_cluster_r[cri_i] = node_to_cluster[item_id_r];
                    cluster_to_cri[cri_to_cluster_r[cri_i]] = cri_i;//node_to_cluster[item_id_r]
                    cri_pqueue_r.update(cri_i, cri_hops_r[cri_i]);
                }

                jump1: "there is a jump.";
            }
        }
        tt.stop();
        query_time += tt.GetRuntime() * 1000;
//        cout<<"Query time: " <<tt.GetRuntime()<<"s."<<endl;
        //io record
        stats_e = stats_b.get_stats();
        if (stats_e.m_noOfWrites > 0 || stats_e.m_noOfReads > 0) {
            cout << "There are " << stats_e.m_noOfWrites + stats_e.m_noOfReads << " io for stxxl." << endl;
            //cout<<stats_e<<endl;
            EM_IORecord(stats_e);
        }
        mc_HotPool.clear();
    }
    //function of reading partition for one-pass algorithms with io optimization
    void EMMCGraph::GraphMCReadCluster_New_Bi(MCHotPool<VectorMCEdgesEMTuple_IO_Bi> & mcHotPool,int target_p_id,vector<int8_t>& cluster_to_bi,int8_t direction,vector<bool>& clusterRead){ //NodeId& ID1,
//        cout << "target_p_id: " << target_p_id << endl;
        int evict_p_id = mcHotPool.LRU_IMCluster.put(target_p_id,1);
        //// evict old partition for new one
        if(evict_p_id != -1){//if it is necessary to evict old partition
//            cout<<"There is an eviction!"<<endl;
            assert(get<0>(mcHotPool.clusterStatus[evict_p_id]) == 'I');
            unsigned int evict_p_size = get<2>(mcHotPool.clusterStatus[evict_p_id]);
            uint temp_sz = get<1>(mcHotPool.clusterStatus[evict_p_id])*alpha_bi;//
            assert(temp_sz>0);
            assert(evict_p_size>0);

            if(evict_p_size > temp_sz){//if the remaining size is too large or no node remain, evict all
//                get<2>(mcHotPool.clusterStatus[evict_p_id]) = 0;
                get<0>(mcHotPool.clusterStatus[evict_p_id]) = 'D';
//                cluster_to_bi[evict_p_id] = 0;
                for(int i=partitions_info[evict_p_id];i<partitions_info[evict_p_id+1];++i){
                    get<2>(EMEdgesIndex_Bi[i]) = -1;//hard delete
                    get<3>(EMEdgesIndex_Bi[i]) = -1;
                }
//                for(auto it=mcHotPool.HotPools[evict_p_id].begin();it!=mcHotPool.HotPools[evict_p_id].end();++it){
//                    get<2>(EMEdgesIndex_Bi[it->ID1]) = -1;//hard delete
//                    get<3>(EMEdgesIndex_Bi[it->ID1]) = -1;
//                }
                mcHotPool.HotPools[evict_p_id].clear();//clear
                mcHotPool.HotPools[evict_p_id].shrink_to_fit();//free memory
//                mcHotPool.im_partitions.erase(evict_p_id);
            }else if(evict_p_size > MCEDGE_PER_BLK){//if the remaining size is moderate, store vertices in temp file; if the remaining size is smaller than block size, remain the vertices in hop pool
                assert(mcHotPool.HotPools[evict_p_id].size()>0);
                /// write to stxxl vector
                int temp_id = -1;
                Timer tt;
                tt.start();
                //The partition that is stored in stxxl vector has limited information for the other search, therefore, it may need to be re-read from the disk. There are two cases that the partition need to be written into stxxl vector: 1) the first time for one directional search to write, i.e. get<0>(mcHotPool.clusterStatus[evict_p_id]) = 'I'; 2) the second time for the other directional search to write, in this case, all the other useful vertices that are not in stxxl vector will be written into stxxl vector.
                if(!clusterRead[evict_p_id]){//the first time to write
                    if(cluster_to_bi[evict_p_id] == FORWARD){//if the evicting partition was read by forward search
                        for(auto it=mcHotPool.HotPools[evict_p_id].begin();it!=mcHotPool.HotPools[evict_p_id].end();++it){
                            if(!get<0>(EMEdgesIndex_Bi[it->ID1])){//if false, i.e. if the node needs to be pushed into external vector
                                if(it->ID1 != temp_id){
                                    get<2>(EMEdgesIndex_Bi[it->ID1]) = mcHotPool.MCEdges_EM.size();
                                    if(temp_id != -1){
                                        get<3>(EMEdgesIndex_Bi[temp_id]) = mcHotPool.MCEdges_EM.size();
                                    }
                                    temp_id = it->ID1;
                                }
                                mcHotPool.MCEdges_EM.push_back(*it);
                            }else{//if the vertex should be evicted
                                get<2>(EMEdgesIndex_Bi[it->ID1]) = -1;//hard delete
                                get<3>(EMEdgesIndex_Bi[it->ID1]) = -1;
                            }
                        }
                    }else if(cluster_to_bi[evict_p_id] == REVERSE){//if the evicting partition was read by reverse search
                        for(auto it=mcHotPool.HotPools[evict_p_id].begin();it!=mcHotPool.HotPools[evict_p_id].end();++it){
                            if(!get<1>(EMEdgesIndex_Bi[it->ID1])){//if false, i.e. if the node needs to be pushed into external vector
                                if(it->ID1 != temp_id){
                                    get<2>(EMEdgesIndex_Bi[it->ID1]) = mcHotPool.MCEdges_EM.size();
                                    if(temp_id != -1){
                                        get<3>(EMEdgesIndex_Bi[temp_id]) = mcHotPool.MCEdges_EM.size();
                                    }
                                    temp_id = it->ID1;
                                }
                                mcHotPool.MCEdges_EM.push_back(*it);
                            }else{//if the vertex should be evicted
                                get<2>(EMEdgesIndex_Bi[it->ID1]) = -1;//hard delete
                                get<3>(EMEdgesIndex_Bi[it->ID1]) = -1;
                            }
                        }
                    }
                    get<3>(EMEdgesIndex_Bi[temp_id]) = mcHotPool.MCEdges_EM.size();
                }else{//if this is the second time to write
//                    cout << "The second write of partition " <<evict_p_id<<endl;
                    if(cluster_to_bi[evict_p_id] == FORWARD){//if the evicting partition was read by forward search
                        for(auto it=mcHotPool.HotPools[evict_p_id].begin();it!=mcHotPool.HotPools[evict_p_id].end();++it){
                            if(!get<1>(EMEdgesIndex_Bi[it->ID1])){//if it has never been evicted by reverse search, i.e. it already exists in stxxl vector
                                continue;
                            }
                            if(!get<0>(EMEdgesIndex_Bi[it->ID1])){//if false for forward search, i.e. if the node needs to be pushed into external vector
                                if(it->ID1 != temp_id){
                                    get<2>(EMEdgesIndex_Bi[it->ID1]) = mcHotPool.MCEdges_EM.size();
                                    if(temp_id != -1){
                                        get<3>(EMEdgesIndex_Bi[temp_id]) = mcHotPool.MCEdges_EM.size();
                                    }
                                    temp_id = it->ID1;
                                }
                                mcHotPool.MCEdges_EM.push_back(*it);
                            }else{//if the vertex should be evicted
                                get<2>(EMEdgesIndex_Bi[it->ID1]) = -1;//hard delete
                                get<3>(EMEdgesIndex_Bi[it->ID1]) = -1;
                            }
                        }
                    }else if(cluster_to_bi[evict_p_id] == REVERSE){//if the evicting partition was read by reverse search
                        for(auto it=mcHotPool.HotPools[evict_p_id].begin();it!=mcHotPool.HotPools[evict_p_id].end();++it){
                            if(!get<0>(EMEdgesIndex_Bi[it->ID1])){//if it has never been evicted by reverse search, i.e. it already exists in stxxl vector
                                continue;
                            }
                            if(!get<1>(EMEdgesIndex_Bi[it->ID1])){//if false, i.e. if the node needs to be pushed into external vector
                                if(it->ID1 != temp_id){
                                    get<2>(EMEdgesIndex_Bi[it->ID1]) = mcHotPool.MCEdges_EM.size();
                                    if(temp_id != -1){
                                        get<3>(EMEdgesIndex_Bi[temp_id]) = mcHotPool.MCEdges_EM.size();
                                    }
                                    temp_id = it->ID1;
                                }
                                mcHotPool.MCEdges_EM.push_back(*it);
                            }else{//if the vertex should be evicted
                                get<2>(EMEdgesIndex_Bi[it->ID1]) = -1;//hard delete
                                get<3>(EMEdgesIndex_Bi[it->ID1]) = -1;
                            }
                        }
                    }
                    get<3>(EMEdgesIndex_Bi[temp_id]) = mcHotPool.MCEdges_EM.size();
                }

                tt.stop();
                io_time += tt.GetRuntime() * 1000;

                clusterRead[evict_p_id] = true;
                get<0>(mcHotPool.clusterStatus[evict_p_id]) = 'S';
                get<2>(mcHotPool.clusterStatus[evict_p_id]) = 0;
                mcHotPool.HotPools[evict_p_id].clear();//clear
                mcHotPool.HotPools[evict_p_id].shrink_to_fit();//free memory
//                mcHotPool.im_partitions.erase(evict_p_id);
            }
        }
        /// read partition from original disk file to internal vector
        /// there is three cases: 1) the partition is firstly read, which can be judged by get<1>(mcHotPool.clusterStatus[target_p_id]) == 0; 2) the partition has been read before and it is the same search demands for it, which can be judged by get<1>(mcHotPool.clusterStatus[target_p_id]) != 0 and flag_reverse == false; 3) the partition has been read before and it is the other search demands for it, which can be judged by get<1>(mcHotPool.clusterStatus[target_p_id]) != 0 and flag_reverse == true
        //check status of storage
        char filePath[200];
        int u = 0, u_deg;
        int v, weight;
        MCEdgeT temp_edge;

        strcpy(filePath, DataPath);
        strcat(filePath, dataset.c_str());
        strcat(filePath, "/");
        strcat(filePath, partition_type.c_str());
        strcat(filePath, "/Partitions/");
        strcat(filePath, dataset.c_str());
        strcat(filePath, "_Partition_");
        strcat(filePath, to_string(PartitionSize).c_str());
        strcat(filePath, "_");
        strcat(filePath, to_string(target_p_id).c_str());
        strcat(filePath, ".bin");

        Timer tt;
        tt.start();
        FILE* file = fopen(filePath, "rb");
        assert(file != NULL);
        ReadBuffer read_io(file);
        bool flag_sizeUpdate = false;
        int temp_id = 0;

        if(get<1>(mcHotPool.clusterStatus[target_p_id]) == 0){//if this is the first read, judged by the original size of the cluster
            flag_sizeUpdate = true;
            while(!read_io.isend) {//if not read out
                read_io.read(&u);//get node id
                temp_edge.ID1 = u-1;
                if(u != temp_id){
                    get<2>(EMEdgesIndex_Bi[u-1]) = mcHotPool.HotPools[target_p_id].size();
                    if(temp_id != 0){
                        get<3>(EMEdgesIndex_Bi[temp_id-1]) = mcHotPool.HotPools[target_p_id].size();
                    }
                    temp_id = u;
                }

                assert(vertex_cri_Bi[u-1].first == -1);
                assert(vertex_cri_Bi[u-1].second == -1);
                vertex_cri_Bi[u-1].first = num_of_cri;
                vertex_cri_Bi[u-1].second = num_of_cri;

                read_io.read(&u_deg);//get node degree

                for (int i = 0; i < u_deg; ++i) {
                    read_io.read(&v);//read end node
                    temp_edge.ID2 = v-1;
                    for (int j = 0; j < NUM_OF_CRITERIA; ++j) {
                        read_io.read(&weight);
                        if(j < num_of_cri){
                            temp_edge.putW(j,weight);
                        }
                    }
                    mcHotPool.HotPools[target_p_id].emplace_back(temp_edge);
                }
                get<2>(mcHotPool.clusterStatus[target_p_id]) += u_deg;
                get<1>(mcHotPool.clusterStatus[target_p_id]) += u_deg;
            }
            io_num += read_io.io_number;//update #IO
        }else{//if it is not the first read, read all
            cout<<"The second read for partition "<<target_p_id<<endl;
            mcHotPool.HotPools[target_p_id].clear();
            if(direction == FORWARD) {//if the partition now is read by forward search
                if (cluster_to_bi[target_p_id] == REVERSE) {//if it has been read by the other search, update the edge number in memory
                    get<2>(mcHotPool.clusterStatus[target_p_id]) = 0;
                }
                while (!read_io.isend) {//if not read out
                    read_io.read(&u);//get node id
                    temp_edge.ID1 = u - 1;

                    if (u != temp_id) {
                        get<2>(EMEdgesIndex_Bi[u - 1]) = mcHotPool.HotPools[target_p_id].size();

                        if (temp_id != 0) {
                            get<3>(EMEdgesIndex_Bi[temp_id - 1]) = mcHotPool.HotPools[target_p_id].size();
                        }
                        temp_id = u;
                    }

                    read_io.read(&u_deg);//get node degree

                    for (int i = 0; i < u_deg; ++i) {
                        read_io.read(&v);//read end node
                        temp_edge.ID2 = v - 1;
                        for (int j = 0; j < NUM_OF_CRITERIA; ++j) {
                            read_io.read(&weight);
                            if (j < num_of_cri) {
                                temp_edge.putW(j, weight);
                            }
                        }
                        mcHotPool.HotPools[target_p_id].emplace_back(temp_edge);
                    }
                    if (cluster_to_bi[target_p_id] == REVERSE) {//if it has been read by the other search
                        if (!get<0>(EMEdgesIndex_Bi[u - 1]))//if false, if it has not been evicted
                            get<2>(mcHotPool.clusterStatus[target_p_id]) += u_deg;
                    }
                }
            }else if(direction == REVERSE) {//if the partition now is read by reverse search
                if(cluster_to_bi[target_p_id] == FORWARD) {//if it has been read by the other search
                    get<2>(mcHotPool.clusterStatus[target_p_id]) = 0;
                }
                while(!read_io.isend) {//if not read out
                    read_io.read(&u);//get node id
                    temp_edge.ID1 = u-1;

                    if(u != temp_id){
                        get<2>(EMEdgesIndex_Bi[u-1]) = mcHotPool.HotPools[target_p_id].size();

                        if(temp_id != 0){
                            get<3>(EMEdgesIndex_Bi[temp_id-1]) = mcHotPool.HotPools[target_p_id].size();
                        }
                        temp_id = u;
                    }

                    read_io.read(&u_deg);//get node degree

                    for (int i = 0; i < u_deg; ++i) {
                        read_io.read(&v);//read end node
                        temp_edge.ID2 = v-1;
                        for (int j = 0; j < NUM_OF_CRITERIA; ++j) {
                            read_io.read(&weight);
                            if(j < num_of_cri){
                                temp_edge.putW(j,weight);
                            }
                        }
                        mcHotPool.HotPools[target_p_id].emplace_back(temp_edge);
                    }
                    if(cluster_to_bi[target_p_id] == FORWARD){//if it has been read by the other search
                        if(!get<1>(EMEdgesIndex_Bi[u-1]))//if false, if it has not been evicted
                            get<2>(mcHotPool.clusterStatus[target_p_id]) += u_deg;
                    }
                }
            }
            io_num += 1;//update #IO
        }

        get<3>(EMEdgesIndex_Bi[u-1]) = mcHotPool.HotPools[target_p_id].size();
//        mcHotPool.im_partitions.insert(target_p_id);
        //io_num += read_io.io_number;//update #IO
        fclose(file);
        tt.stop();
        io_time += tt.GetRuntime() * 1000;
        cluster_to_IO[target_p_id] = read_io.io_number;
//        if(get<0>(mcHotPool.clusterStatus[target_p_id]) != 'D'){
//            cout << "The second read of partition " <<target_p_id<<endl;
//        }
        get<0>(mcHotPool.clusterStatus[target_p_id]) = 'I';
    }
    //Function of judge whether the criteria have finished
    bool EMMCGraph::EM_JudgeEmpty(vector<PriorityQueue*> &em_mc_pqueue){//IM
        for(int i=0;i<num_of_cri;++i){
            if(!mc_finished[i] && !em_mc_pqueue[i]->empty())//if exist anyone not empty and anyone unfinished, return false
                return false;
        }
        return true;
    }
    //Function of judge whether the criteria have finished
    bool EMMCGraph::EM_JudgeEmptyBi(vector<PriorityQueue*> &em_mc_pqueue, vector<PriorityQueue2*> &em_mc_pqueue_r){
        for(int i=0;i<num_of_cri;i++){
            if(!mc_finished[i] && (!em_mc_pqueue[i]->empty() || !em_mc_pqueue_r[i]->empty()))//if exist anyone not empty and anyone unfinished, return false
                return false;
        }
        return true;
    }
    //function of printing basic information
    void EMMCGraph::InfoPrint() const {
        cout << "Dataset: " << dataset << endl;
        cout << "Number of criteria: " << num_of_cri << endl;
        cout << "Query type: " << query_type << endl;
        cout << "Run times: " << run_times << endl;
        cout << "Total memory size limitation: " << memSize << " MB." <<endl;
        cout << "Minimal memory size for graph date: " << MinMemThreshold << " MB." <<endl;
        cout << "Block size: " << Block_SZ/1024 << " KB." <<endl;
        cout << "Partition size: " << PartitionSize << " KB." <<endl;
        cout << "Partition strategy: " << partition_type << "." <<endl;
        cout << "PageNumber4_Dijk: " << PageNumber4_Dijk << endl;
    }
    //function of checking if the memory setting is feasible
    void EMMCGraph::MemoryCheck(int algo_choice, uint mem, bool& flag_exit) {
        uint basic_mem;

        ClusterInfoLoad(false);
        cout << "Node number: "<<node_num<<"\tEdge number: "<<edge_num<<endl;
        /// basic algorithm memory for in-memory algorithms
        basic_mem = node_num * 17 / (1024 * 1024);
        cout << "*** The basic memory consumption of in-memory Dijkstra is: " << basic_mem << " MB. ***" << endl;

        switch (algo_choice){
            case EMDijk:{
                basic_mem = node_num*12/(1024*1024) + WeightPowMax*QMemory/(1024*1024) + PQMemory; //
//                basic_mem = node_num*12/(1024*1024) + WeightPowMax/2 + PQMemory; //
                cout << "*** The basic memory consumption of EM_Dijk is: " << basic_mem << " MB. ***"<<  endl;
                if(mem <= basic_mem){
                    cout << "!!!The memory setting is too small!\n" << endl;
                }else{
                    memGraph = (mem - basic_mem > MinMemThreshold ? mem - basic_mem : MinMemThreshold);
                    cout << "*** The memory for graph is: " << memGraph << " MB. ***" << endl;
                    if(memGraph != MemGraph_EMDijk){
                        cout << "!!!Please change the memory size for graph data to " << memGraph << " in head.h manually!\n" << endl;
                        flag_exit = true;
                    }
                }
                break;
            }
            case DijkstraIO:{
                basic_mem = node_num*18/(1024*1024) + PQMemory;
                cout << "*** The basic memory consumption of Dijkstra_IO is: " << basic_mem << " MB. ***"<<  endl;
                if(mem <= basic_mem){
                    cout << "!!!The memory setting is too small!\n" << endl;
                }else{
                    memGraph = (mem - basic_mem > MinMemThreshold ? mem - basic_mem : MinMemThreshold);
                    cout << "*** The memory for graph is: " << memGraph << " MB. ***" << endl;
                    if(memGraph != MemGraph_DijkIO){
                        cout << "!!!Please change the memory size for graph data to " << memGraph << " in head.h manually!\n" << endl;
                        flag_exit = true;
                    }
                }
                break;
            }
            case BiDijkstraIO:{
                basic_mem = node_num*23/(1024*1024) + PQMemory;
                cout << "*** The basic memory consumption of BiDijkstra_IO is: " << basic_mem << " MB. ***"<<  endl;
                if(mem <= basic_mem){
                    cout << "!!!The memory setting is too small!\n" << endl;
                }else{
                    memGraph = (mem - basic_mem > MinMemThreshold ? mem - basic_mem : MinMemThreshold);
                    cout << "*** The memory for graph is: " << memGraph << " MB. ***" << endl;
                    if(memGraph != MemGraph_BiDijkIO){
                        cout << "!!!Please change the memory size for graph data to " << memGraph << " in head.h manually!\n" << endl;
                        flag_exit = true;
                    }
                }
                break;
            }
            case OneHopNoIO:{
                basic_mem = node_num*(num_of_cri*5 + 11)/(1024*1024) + PQMemory*num_of_cri;
                cout << "*** The basic memory consumption of OneHopNoIO is: " << basic_mem << " MB. ***"<<  endl;
                if(mem <= basic_mem){
                    cout << "!!!The memory setting is too small!\n" << endl;
                }else{
                    memGraph = (mem - basic_mem > MinMemThreshold ? mem - basic_mem : MinMemThreshold);
                    cout << "*** The memory for graph is: " << memGraph << " MB. ***" << endl;
                    if(memGraph != MemGraph_NoIO){
                        cout << "!!!Please change the memory size for graph data to " << memGraph << " in head.h manually!\n" << endl;
                        flag_exit = true;
                    }
                }
                break;
            }
            case MultiHopsNoIO:{
                basic_mem = node_num*(num_of_cri*5 + 11)/(1024*1024) + PQMemory*num_of_cri;
                cout << "*** The basic memory consumption of MultiHopsNoIO is: " << basic_mem << " MB. ***"<<  endl;
                if(mem <= basic_mem){
                    cout << "!!!The memory setting is too small!\n" << endl;
                }else{
                    memGraph = (mem - basic_mem > MinMemThreshold ? mem - basic_mem : MinMemThreshold);
                    cout << "*** The memory for graph is: " << memGraph << " MB. ***" << endl;
                    if(memGraph != MemGraph_NoIO){
                        cout << "!!!Please change the memory size for graph data to " << memGraph << " in head.h manually!\n" << endl;
                        flag_exit = true;
                    }
                }
                break;
            }
            case BiMultiHopsNoIO:{
                basic_mem = node_num*(num_of_cri*10 + 11)/(1024*1024) + 2*PQMemory*num_of_cri;
                cout << "*** The basic memory consumption of BiMultiHopsNoIO is: " << basic_mem << " MB. ***"<<  endl;
                if(mem <= basic_mem){
                    cout << "!!!The memory setting is too small!\n" << endl;
                }else{
                    memGraph = (mem - basic_mem > MinMemThreshold ? mem - basic_mem : MinMemThreshold);
                    cout << "*** The memory for graph is: " << memGraph << " MB. ***" << endl;
                    if(memGraph != MemGraph_NoIO_Bi){
                        cout << "!!!Please change the memory size for graph data to " << memGraph << " in head.h manually!\n" << endl;
                        flag_exit = true;
                    }
                }
                break;
            }
            case OneHop:{
                basic_mem = node_num*(num_of_cri*5 + 12)/(1024*1024) + PQMemory*num_of_cri;
                cout << "*** The basic memory consumption of OneHop is: " << basic_mem << " MB. ***"<<  endl;
                if(mem <= basic_mem){
                    cout << "!!!The memory setting is too small!\n" << endl;
                }else{
                    memGraph = (mem - basic_mem > MinMemThreshold ? mem - basic_mem : MinMemThreshold);
                    cout << "*** The memory for graph is: " << memGraph << " MB. ***" << endl;
                    if(memGraph != MemGraph_IO){
                        cout << "!!!Please change the memory size for graph data to " << memGraph << " in head.h manually!\n" << endl;
                        flag_exit = true;
                    }
                }
                break;
            }
            case MultiHops:{
                basic_mem = node_num*(num_of_cri*5 + 12)/(1024*1024) + PQMemory*num_of_cri;
                cout << "*** The basic memory consumption of MultiHops is: " << basic_mem << " MB. ***"<<  endl;
                if(mem <= basic_mem){
                    cout << "!!!The memory setting is too small!\n" << endl;
                }else{
                    memGraph = (mem - basic_mem > MinMemThreshold ? mem-basic_mem : MinMemThreshold);
                    cout << "*** The memory for graph is: " << memGraph << " MB. ***" << endl;
                    if(memGraph != MemGraph_IO){
                        cout << "!!!Please change the memory size for graph data to " << memGraph << " in head.h manually!\n" << endl;
                        flag_exit = true;
                    }
                }
                break;
            }
            case BiMultiHops:{
                basic_mem = node_num*(num_of_cri*10 + 12)/(1024*1024) + 2*PQMemory*num_of_cri;
                cout << "*** The basic memory consumption of BiMultiHops is: " << basic_mem << " MB. ***"<<  endl;
                if(mem <= basic_mem){
                    cout << "!!!The memory setting is too small!\n" << endl;
                }else{
                    memGraph = (mem - basic_mem > MinMemThreshold ? mem - basic_mem : MinMemThreshold);
                    cout << "*** The memory for graph is: " << memGraph << " MB. ***" << endl;
                    if(memGraph != MemGraph_IO_Bi){
                        cout << "!!!Please change the memory size for graph data to " << memGraph << " in head.h manually!\n" << endl;
                        flag_exit = true;
                    }
                }
                break;
            }
            default:
                cout<<"\nwrong!"<<endl; break;
        }
    }
    //function of converting edge weight to category
    int EMMCGraph::Weight_to_category(int w){
        return floor(log2(w));
    }
    //function for recording the io performance of EM_Dijk
    void EMMCGraph::EM_IORecord(Stats& a) {
        double ratio_r = (double)a.m_aveBlockSizeR / bufferSize;
        double ratio_w = (double)a.m_aveBlockSizeW / bufferSize;
        io_num = io_num + int(ratio_r * a.m_noOfReads) + int(ratio_w * a.m_noOfWrites);
        io_time += a.m_pioTime * 1000;
    }
    //function for computing the minimal io
    int EMMCGraph::Minimal_IO(){
        int minio = 0;
        for(auto it=set_readCluster.begin();it!=set_readCluster.end();++it){
            minio += cluster_to_IO[*it];
        }
        return minio;
    }

    //Function for traverse parental vertices to get the shortest path
    list<int> EMMCGraph::Dij_getPath(vector<NodeId> &pre, NodeId ID1, NodeId ID2) {
        list<int> list_sp;      //list to record nodes in the shortest path
        list<int> temp_list;
        int pre_temp = ID2;

        while (pre_temp != ID1) {//&& pre_temp != -1
            temp_list.push_back(pre_temp);
            pre_temp = pre[pre_temp];
        }
        if (pre_temp == -1) {
            cout << "Do not find the shortest path!" << endl;
        } else {
            temp_list.push_back(pre_temp);
            list_sp.clear();
            while (!temp_list.empty()) {//reverse the elements
                list_sp.push_back(temp_list.back());
                temp_list.pop_back();
            }
        }
        for(auto it=list_sp.begin();it!=list_sp.end();++it){
            cout<<*it<<"->";
        }
        cout<<endl;
        return list_sp;
    }

    //Function for traverse parental vertices to get the shortest path
    list<int> EMMCGraph::BiDij_getPath(vector<NodeId> &pre, vector<NodeId> &pre_b, NodeId ID1, NodeId ID2, NodeId terminate_id) {
        list<int> list_sp;      //list to record nodes in the shortest path
        list<int> temp_list;
        int pre_temp = terminate_id;
        while (pre_temp != ID1 && pre_temp != -1) {
            temp_list.push_back(pre_temp);
            pre_temp = pre[pre_temp];
        }
        if (pre_temp == -1) {
            cout << "Do not find the shortest path!" << endl;
        }
        //reverse the result
        temp_list.push_back(pre_temp);
        while (!temp_list.empty()) {
            list_sp.push_back(temp_list.back());
            temp_list.pop_back();
        }
        pre_temp = terminate_id;
        while (pre_temp != ID2) {
            pre_temp = pre_b[pre_temp];
            list_sp.push_back(pre_temp);
        }
        if (pre_temp == -1) {
            cout << "Do not find the shortest path!" << endl;
        }
        for(auto it=list_sp.begin();it!=list_sp.end();++it){
            cout<<*it<<"->";
        }
        cout<<endl;
        return list_sp;
    }
}
#endif //MCSPS_EMGRAPH_HPP
