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

    //function for reading partition information
    void EMMCGraph::ClusterInfoLoad(bool ifMap) {//load partition information, ifMap is used to indicate whether to used node-to-cluster mapping
        /*** Load partition information ***/
        string filename = string(DataPath) + dataset + "/" + dataset + "_Partitions_" + to_string(PartitionSize) + "_Info.txt";
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
        partitions_info.assign(partition_number,0);
        cluster_to_IO.assign(partition_number,0);

        if(ifMap){//if it needs to map vertex id to partition id
            node_to_cluster.assign(node_num,0);
            for(int i=0;i<partition_number;++i){
                inFile >> p_id >> ID1 >> ID2;
                partitions_info[p_id] = ID1;
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
    void EMMCGraph::MC_Evaluate_IO(const string& qtype, int algo_choice)
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
        cout << "Node number: " << node_num << "\t\tEdge number: " << edge_num << endl;

        switch (algo_choice) {
            case OneHop:{
                /// execution of MC_OneHop
                cout << "----- Below are results of MC_OneHop algorithm -----" << endl;
                time_OneHop = MC_OnePass(qtype,OneHop);
                cout << "-----------------------------------------------\n" << endl;
                cout << "The total run time of MC_OneHop algorithm is: "<< time_OneHop /1000<<" s."<<endl;
                break;
            }
            case MultiHops:{
                /// execution of MC_MultiHop
                cout << "----- Below are results of MC_MultiHop algorithm -----" << endl;
                time_MultiHops = MC_OnePass(qtype,MultiHops);
                cout << "-----------------------------------------------\n" << endl;
                cout << "The total run time of MC_MultiHop algorithm is: "<< time_MultiHops /1000<<" s."<<endl;
                break;
            }
            case BiMultiHops:{
                /// execution of MC_BiMultiHop
                cout << "----- Below are results of MC_BiMultiHop algorithm -----" << endl;
                time_BiMultiHops = MC_OnePass(qtype,BiMultiHops);
                cout << "-----------------------------------------------\n" << endl;
                break;
            }
            case OneHopNoIO:{
                /// execution of MC_OneHop without IO optimization
                cout << "----- Below are results of MC_OneHop_NoIO algorithm -----" << endl;
                time_OneHopNoIO = MC_OnePass(qtype,OneHopNoIO);
                cout << "-----------------------------------------------\n" << endl;
                cout << "The total run time of MC_OneHop_NoIO algorithm is: "<< time_OneHopNoIO /1000<<" s."<<endl;
                break;
            }
            case MultiHopsNoIO:{
                /// execution of MC_MultiHop without IO optimization
                cout << "----- Below are results of MC_MultiHop_NoIO algorithm -----" << endl;
                time_MultiHopsNoIO = MC_OnePass(qtype,MultiHopsNoIO);
                cout << "-----------------------------------------------\n" << endl;
                cout << "The total run time of MC_MultiHop_NoIO algorithm is: "<< time_MultiHopsNoIO /1000<<" s."<<endl;
                break;
            }
            case BiMultiHopsNoIO:{
                /// execution of MC_BiMultiHop without IO optimization
                cout << "----- Below are results of MC_BiMultiHop_NoIO algorithm -----" << endl;
                time_BiMultiHopsNoIO = MC_OnePass(qtype,BiMultiHopsNoIO);
                cout << "-----------------------------------------------\n" << endl;
                cout << "The total run time of MC_BiMultiHop_NoIO algorithm is: " << time_BiMultiHopsNoIO / 1000 << " s." << endl;
                break;
            }
            default:
                break;
        }
        node_to_cluster.clear();
    }
    //Disk-based MCSP Evaluation: Exp1-Parameters setting
    void EMMCGraph::MC_Evaluate_EXP1(const string& qtype)
    {
        /*--variables about time record--*/
        double time_OneHop = 0;
        double time_MultiHops = 0;
        double time_BiMultiHops = 0;

        assert(BlockPerPage*Block_SZ <= PartitionSize*1024);//page size should be smaller than partition size
        /*** Read graph data ***/
        ClusterInfoLoad(true);//load partition information, using node-to-cluster mapping
        cout << "Node number: " << node_num << "\t\tEdge number: " << edge_num << endl;

        /*** OnePass algorithms ***/
        ///// execution of MC_OneHop
        cout << "----- Below are results of MC_OneHop algorithm -----" << endl;
        cout << "The graph data memory size is " << MemGraph_IO << " MB." << endl;
        time_OneHop = MC_OnePass(qtype,OneHop);
        cout << "-----------------------------------------------\n" << endl;
        /// execution of MC_MultiHop
        cout << "----- Below are results of MC_MultiHop algorithm -----" << endl;
        cout << "The graph data memory size is " << MemGraph_IO << " MB." << endl;
        time_MultiHops = MC_OnePass(qtype,MultiHops);
        cout << "-----------------------------------------------\n" << endl;
        /// execution of MC_BiMultiHop
        cout << "----- Below are results of MC_BiMultiHop algorithm -----" << endl;
        cout << "The graph data memory size is " << MemGraph_IO_Bi << " MB." << endl;
        time_BiMultiHops = MC_OnePass(qtype,BiMultiHops);
        cout << "-----------------------------------------------\n" << endl;
        cout << "The total run time of MC_OneHop algorithm is: "<< time_OneHop /1000<<" s."<<endl;
        cout << "The total run time of MC_MultiHop algorithm is: "<< time_MultiHops /1000<<" s."<<endl;
        cout << "The total run time of MC_BiMultiHop algorithm is: " << time_BiMultiHops / 1000 << " s." << endl;
        node_to_cluster.clear();
    }
    //Disk-based MCSP Evaluation: Exp2-Effectiveness of optimizations
    void EMMCGraph::MC_Evaluate_EXP2(const string& qtype)
    {
        /*--variables about time record--*/
        double time_OneHopNoIO = 0;
        double time_MultiHopsNoIO = 0;
        double time_BiMultiHopsNoIO = 0;
        double time_OneHop = 0;
        double time_MultiHops = 0;
        double time_BiMultiHops = 0;

        assert(BlockPerPage*Block_SZ <= PartitionSize*1024);//page size should be smaller than partition size
        /*** Read graph data ***/
        ClusterInfoLoad(true);//load partition information, using node-to-cluster mapping
        cout << "Node number: " << node_num << "\t\tEdge number: " << edge_num << endl;

        /*** OnePass algorithms ***/
        /// execution of MC_OneHop without IO optimization
        cout << "----- Below are results of MC_OneHop_NoIO algorithm -----" << endl;
        cout << "The graph data memory size is " << MemGraph_NoIO << " MB." << endl;
        time_OneHopNoIO = MC_OnePass(qtype,OneHopNoIO);
        cout << "-----------------------------------------------\n" << endl;
        /// execution of MC_MultiHop without IO optimization
        cout << "----- Below are results of MC_MultiHop_NoIO algorithm -----" << endl;
        cout << "The graph data memory size is " << MemGraph_NoIO << " MB." << endl;
        time_MultiHopsNoIO = MC_OnePass(qtype, MultiHopsNoIO);
        cout << "-----------------------------------------------\n" << endl;
        ///// execution of MC_BiMultiHop without IO optimization
        cout << "----- Below are results of MC_BiMultiHop_NoIO algorithm -----" << endl;
        cout << "The graph data memory size is " << MemGraph_NoIO_Bi << " MB." << endl;
        time_BiMultiHopsNoIO = MC_OnePass(qtype,BiMultiHopsNoIO);
        cout << "-----------------------------------------------\n" << endl;
        

        /// execution of MC_OneHop
        cout << "----- Below are results of MC_OneHop algorithm -----" << endl;
        cout << "The graph data memory size is " << MemGraph_IO << " MB." << endl;
        time_OneHop = MC_OnePass(qtype, OneHop);
        cout << "-----------------------------------------------\n" << endl;
        /// execution of MC_MultiHop
        cout << "----- Below are results of MC_MultiHop algorithm -----" << endl;
        cout << "The graph data memory size is " << MemGraph_IO << " MB." << endl;
        time_MultiHops = MC_OnePass(qtype,MultiHops);
        cout << "-----------------------------------------------\n" << endl;
        /// execution of MC_BiMultiHop
        cout << "The graph data memory size is " << MemGraph_IO_Bi << " MB." << endl;
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
    //Disk-based MCSP Evaluation: Exp3-Comparison with baselines
    void EMMCGraph::MC_Evaluate_EXP3(const string& qtype)
    {
        /*--variables about time record--*/
        //std::chrono::high_resolution_clock::time_point t_s, t_e;//variables for time record
        double time_EMDijk = 0;
        double time_DijkstraIO = 0;
        double time_BiDijkstraIO = 0;
        double time_MultiHops = 0;
        double time_BiMultiHops = 0;
        double time_OneHopNoIO = 0;

        assert(BlockPerPage * Block_SZ <= PartitionSize * 1024);//page size should be smaller than partition size
        /*** Read graph data ***/
        ClusterInfoLoad(true);//load partition information, using node-to-cluster mapping

        cout << "Node number: " << node_num << "\t\tEdge number: " << edge_num << endl;

        /*** MultiPass algorithm ***/
        cout << "----- Below are results of EM_Dijk algorithm -----" << endl;
        cout << "The graph data memory size is " << MemGraph_IO << " MB." << endl;
        time_EMDijk = MC_Multipass(qtype,EMDijk);
        cout << "-----------------------------------------------\n" << endl;
        /// execution of Dijkstra with external vector
        cout << "----- Below are results of Dijkstra_IO algorithm -----" << endl;
        cout << "The graph data memory size is " << MemGraph_DijkIO << " MB." << endl;
        time_DijkstraIO = MC_Multipass(qtype, DijkstraIO);
        cout << "-----------------------------------------------\n" << endl;
        /// execution of BiDijkstra with external vector
        cout << "----- Below are results of BiDijkstra_IO algorithm -----" << endl;
        time_BiDijkstraIO = MC_Multipass(qtype,BiDijkstraIO);
        cout << "-----------------------------------------------\n" << endl;
        /*** OnePass algorithms ***/
        /// execution of MC_MultiHop
        cout << "----- Below are results of MC_MultiHop algorithm -----" << endl;
        cout << "The graph data memory size is " << MemGraph_IO << " MB." << endl;
        time_MultiHops = MC_OnePass(qtype, MultiHops);
        cout << "-----------------------------------------------\n" << endl;
        ///// execution of MC_BiMultiHop
        cout << "The graph data memory size is " << MemGraph_IO_Bi << " MB." << endl;
        cout << "----- Below are results of MC_BiMultiHop algorithm -----" << endl;
        time_BiMultiHops = MC_OnePass(qtype, BiMultiHops);
        cout << "-----------------------------------------------\n" << endl;

        cout << "The total run time of EM-Dijk algorithm is: " << time_EMDijk / 1000 << " s." << endl;
        cout << "The total run time of Dijktra_IO algorithm is: " << time_DijkstraIO / 1000 << " s." << endl;
        cout << "The total run time of BiDijktra_IO algorithm is: "<<time_BiDijkstraIO/1000<<" s."<<endl;
        cout << "The total run time of MC_MultiHop algorithm is: " << time_MultiHops / 1000 << " s." << endl;
        cout << "The total run time of MC_BiMultiHop algorithm is: " << time_BiMultiHops / 1000 << " s." << endl;
        node_to_cluster.clear();
    }
    //Disk-based MCSP Evaluation: Exp4-Effect of the number of criteria; Exp5-Effect of memory size; Exp6-Effect of partition size; Exp7-Effect of query types
    void EMMCGraph::MC_Evaluate_EXP456(const string& qtype)
    {
        /*--variables about time record--*/
        double time_DijkstraIO = 0;
        double time_BiDijkstraIO = 0;
        double time_MultiHops = 0;
        double time_BiMultiHops = 0;

        assert(BlockPerPage*Block_SZ <= PartitionSize*1024);//page size should be smaller than partition size
        /*** Read graph data ***/
        ClusterInfoLoad(true);//load partition information, using node-to-cluster mapping
        cout << "Node number: " << node_num << "\t\tEdge number: " << edge_num << endl;

        /*** MultiPass algorithm ***/
        /// execution of Dijkstra with external vector
        cout << "----- Below are results of Dijkstra_IO algorithm -----" << endl;
        cout << "The graph data memory size is " << MemGraph_DijkIO << " MB." << endl;
        time_DijkstraIO = MC_Multipass(qtype,DijkstraIO);
        cout << "-----------------------------------------------\n" << endl;
        /// execution of BiDijkstra with external vector
        cout << "----- Below are results of BiDijkstra_IO algorithm -----" << endl;
        time_BiDijkstraIO = MC_Multipass(qtype,BiDijkstraIO);
        cout << "-----------------------------------------------\n" << endl;
        /*** OnePass algorithms ***/
        /// execution of MC_MultiHop
        cout << "----- Below are results of MC_MultiHop algorithm -----" << endl;
        cout << "The graph data memory size is " << MemGraph_IO << " MB." << endl;
        time_MultiHops = MC_OnePass(qtype,MultiHops);
        cout << "-----------------------------------------------\n" << endl;
        //// execution of MC_BiMultiHop
        cout << "----- Below are results of MC_BiMultiHop algorithm -----" << endl;
        cout << "The graph data memory size is " << MemGraph_IO_Bi << " MB." << endl;
        time_BiMultiHops = MC_OnePass(qtype,BiMultiHops);
        cout << "-----------------------------------------------\n" << endl;

        cout << "The total run time of Dijktra_IO algorithm is: "<<time_DijkstraIO/1000<<" s."<<endl;
        cout << "The total run time of BiDijktra_IO algorithm is: "<<time_BiDijkstraIO/1000<<" s."<<endl;
        cout << "The total run time of MC_MultiHop algorithm is: "<< time_MultiHops /1000<<" s."<<endl;
        cout << "The total run time of MC_BiMultiHop algorithm is: " << time_BiMultiHops / 1000 << " s." << endl;
        node_to_cluster.clear();
    }
    //Function for multipass algorithms
    double EMMCGraph::MC_Multipass(const string & qtype, int algo_choice) {
        string LongDis,MediumDis,ShortDis,RandomDis;
        double total_time = 0;
        ///File path of OD pairs
        LongDis = string(DataPath) + dataset + "/" + dataset + "_OD_LongDis_" + to_string(PartitionSize) + ".txt";
        MediumDis = string(DataPath) + dataset + "/" + dataset + "_OD_MediumDis_"+ to_string(PartitionSize) +".txt";
        ShortDis = string(DataPath) + dataset + "/" + dataset + "_OD_ShortDis_"+ to_string(PartitionSize) +".txt";
        RandomDis = string(DataPath) +  dataset + "/" + dataset + "_OD_Random_"+ to_string(PartitionSize) +".txt";

        ///Shortest path query processing
        if(qtype == "S"){
            QueryType = "S";
            total_time = MC_Multipass_one(ShortDis,algo_choice);//Efficiency evaluation on short distance
        }else if(qtype == "M"){
            QueryType = "M";
            total_time = MC_Multipass_one(MediumDis,algo_choice);//Efficiency evaluation on medium distance
        }else if(qtype == "L"){
            QueryType = "L";
            total_time = MC_Multipass_one(LongDis,algo_choice);//Efficiency evaluation on long distance
        }else if(qtype == "R"){
            QueryType = "S";
            total_time += MC_Multipass_one(ShortDis, algo_choice);//Efficiency evaluation on short distance
            QueryType = "M";
            total_time += MC_Multipass_one(MediumDis, algo_choice);//Efficiency evaluation on medium distance
            QueryType = "L";
            total_time += MC_Multipass_one(LongDis, algo_choice);//Efficiency evaluation on long distance
        }
        return total_time;
    }
    //Function for multi-pass MCSP
    double EMMCGraph::MC_Multipass_one(const string& filename, int algo_choice) {
        string fname;
        vector<pair<int, int>> ODpair;  //OD pairs for querying
        int num, ID1, ID2;
        uint ave_visited = 0;    //average number of visited vertices
        uint ave_io_num = 0;  //average IO number
        double ave_p_num = 0;   //average read partition number
        double ave_time = 0;    //average time
        double ave_io_time = 0; //average io time
        vector<Distance> ave_cost(num_of_cri,0);    //average shortest path distance
        query_time = 0;     //overall query time
        double ave_min_io = 0;  //average minimal io number
        double ave_min_cluster = 0; //average minimal partition number need to read


        ///Read OD pairs
        ifstream inFile(filename, ios::in);
        if (!inFile) {
            cout << "File opening failed." << endl;
            assert(inFile);
        }
        inFile >> num;
        for (int i = 0; i < run_times; ++i) {
            inFile >> ID1 >> ID2;
            ODpair.emplace_back(make_pair(ID1, ID2));
        }
        inFile.close();

        /// Preprocessing
        if(algo_choice == EMDijk){
            EM_Preprocess();//Preprocessing for EM_Dijk
        }else if(algo_choice == DijkstraIO){
            uint basic_mem = node_num*18/(1024*1024) + 32;
            cout << "*** The basic memory consumption of Dijkstra_IO is: " << basic_mem << " MB. ***"<<  endl;
        }else if(algo_choice == BiDijkstraIO){
            uint basic_mem = node_num*23/(1024*1024) + 32;
            cout << "*** The basic memory consumption of BiDijkstra_IO is: " << basic_mem << " MB. ***"<<  endl;
        }
        cout << "Query type: " << QueryType << " \tRun times: " << run_times << endl;
        /// MCSPs Query processing
        for (int i = 0; i < run_times; ++i) {
            ID1 = ODpair[i].first;
            ID2 = ODpair[i].second;
//            cout << "Query " << i << ": " << ID1 << " " << ID2 << endl;
            //Initiation for each round
            CommonInitiation_IO();
            if(algo_choice == EMDijk){
                ///Shortest path query one criterion by one criterion
                //statistics of IO
                Stats stats_b, stats_e;
                stats_b.start();
//                Timer tt3;
//                tt3.start();
                for(sc_i=0;sc_i<num_of_cri;++sc_i) {
                    ave_cost[sc_i] += EM_Dijk(ID1 - 1, ID2 - 1);
                    //reset done bit to 0
                    for (int j = 0; j < node_num; ++j) {
                        node_to_category[sc_i][j].set(0, false);
                    }
                }
                //io record
//                tt3.stop();
//                cout << "Query time: " << tt3.GetRuntime() << " s." << endl;
                stats_e = stats_b.get_stats();
                if (stats_e.m_noOfWrites > 0 || stats_e.m_noOfReads > 0) {
                    cout << "There are " << stats_e.m_noOfWrites + stats_e.m_noOfReads << " io for stxxl." << endl;
                    //cout<<stats_e<<endl;
                }
            }else if(algo_choice == DijkstraIO){
                //statstics of IO
                Stats stats_b, stats_e;
                stats_b.start();
//                Timer tt3;
//                tt3.start();
                ///Shortest path query one criterion by one criterion
                for(sc_i=0;sc_i<num_of_cri;++sc_i) {
                    EMEdgesIndex.assign(node_num,make_tuple(false,-1,-1));
                    VectorMCEdgesEMTuple_DijkIO EMMCEdges;
                    ave_cost[sc_i] += Dijkstra_IO(ID1 - 1, ID2 - 1,EMMCEdges);
                    EMEdgesIndex.clear();
                    EMMCEdges.clear();
                }
//                tt3.stop();
//                cout << "Query time: " << tt3.GetRuntime() << " s." << endl;
                //io record
                stats_e = stats_b.get_stats();
                if (stats_e.m_noOfWrites > 0 || stats_e.m_noOfReads > 0) {
                    cout << "There are " << stats_e.m_noOfWrites + stats_e.m_noOfReads << " io for stxxl." << endl;
                    //cout<<stats_e<<endl;
                }
            }else if(algo_choice == BiDijkstraIO){
                //statstics of IO
                Stats stats_b, stats_e;
                stats_b.start();
                ///Shortest path query one criterion by one criterion
                for(sc_i=0;sc_i<num_of_cri;++sc_i) {
                    EMEdgesIndex.assign(node_num,make_tuple(false,-1,-1));
                    VectorMCEdgesEMTuple_BiDijkIO EMMCEdges;
                    ave_cost[sc_i] += BiDijkstra_IO(ID1 - 1, ID2 - 1,EMMCEdges);
                    EMEdgesIndex.clear();
                    EMMCEdges.clear();
                }
                //io record
                stats_e = stats_b.get_stats();
                if (stats_e.m_noOfWrites > 0 || stats_e.m_noOfReads > 0) {
                    cout << "There are " << stats_e.m_noOfWrites + stats_e.m_noOfReads << " io for stxxl." << endl;
                    //cout<<stats_e<<endl;
                }
            }
            //record information
            ave_io_num += io_num;
            ave_io_time += io_time;
            ave_p_num += p_num;
            ave_visited += visited_number;
            ave_min_io += Minimal_IO();
            ave_min_cluster += set_readCluster.size();
            //clear
            set_readCluster.clear();
//            set_visited.clear();
        }
        //clear
        if(algo_choice == EMDijk){
            node_to_category.clear();
        }
        ave_time = query_time / run_times;
        ave_visited /= run_times;
        ave_io_num /= run_times;
        ave_io_time /= run_times;
        ave_min_io /= run_times;
        ave_p_num /= run_times;
        ave_min_cluster /= run_times;
        if(algo_choice == EMDijk){
            cout << "Average performance of Multi-Pass EM_Dijk (on " << run_times << " " << QueryType << " OD pairs) is shown below." << endl;
        }else if(algo_choice == DijkstraIO){
            cout << "Average performance of Multi-Pass Dijkstra_IO  (on " << run_times << " " << QueryType << " OD pairs) is shown below." << endl;
        }else if(algo_choice == BiDijkstraIO){
            cout << "Average performance of Multi-Pass BiDijkstra_IO  (on " << run_times << " " << QueryType << " OD pairs) is shown below." << endl;
        }
        for(int i=0;i<num_of_cri;++i) {
            ave_cost[i] /= run_times;
            cout << "Average distance of criteria " << mc_criteria[i] <<": " << ave_cost[i] << endl;
        }
        cout << "! Average query time: " << ave_time/1000 << " s." << endl;
        cout << "! Average io time: " << ave_io_time/1000 << " s." << endl;
        cout << "! Average number of IO: " << ave_io_num << endl;
        cout << "! Average number of nodes visited: " << ave_visited << endl;
        cout << "! Average number of read partitions: " << ave_p_num << endl;
        if(ifOptimal){
            cout << "! Minimal #IO for reading: " << ave_min_io << endl;
            cout << "! Minimal #Partition for reading: " << ave_min_cluster << endl;
        }
        cout << "Average proportion of io time: " << 100*ave_io_time/ave_time << " %" << endl;

        return query_time;
    }

    //// Functions for EM_Dijk
    //Function for the preprocessing of EM_Dijk
    void EMMCGraph::EM_Preprocess(){
        string file_r_graph = string(DataPath) + dataset+"/"+dataset+"_MCEdges_io_" + to_string(PartitionSize)+".txt";
        file_r_graph = string(DataPath) + dataset+"/"+dataset+"_MCEdgesS_" + to_string(PartitionSize)+".txt";
        string file_r_graph1 = string(DataPath) + this->dataset + "/" + this->dataset + "_EM_Dijk_Index_"+ to_string(PartitionSize)+".txt";
        num_pool.assign(num_of_cri,RangeId());
        pow_low.assign(num_of_cri,vector<EdgeWeight>());
        Timer tt;
        tt.start();

        //Construct the edge category index for all vertices
        int ID1,ID2;
        int degree,weight;
        int temp_int;
        int num_cri;
        string temp_cri;
        int temp_id=0;

        ifstream inFile1(file_r_graph1, ios::in);
        if (!inFile1) {//if the file does not exist
            cout << "File opening failed." << endl;
//            assert(inFile1);
//            exit(1);
            //Open file
            ifstream inFile(file_r_graph, ios::in);
            if (!inFile) {
                cout << "File opening failed." << endl;
                exit(1);
            }
            cout << "Graph Data loading for EM_Dijk preprocessing..." << endl;
            inFile >> node_num >> edge_num >> num_cri;
            for(int i=0;i<num_cri;++i){
                inFile >> temp_cri;
            }
            max_w.assign(num_cri,0);
            // get node-to-category mapping
            node_to_category.assign(num_cri,vector<bitset1<WeightPowMax+1,EdgeWeight>>(node_num));//bitset<WeightPowMax+1>
            while(inFile){
                inFile>>ID1>>ID2;
                assert(ID1<=node_num);
                assert(ID2<=node_num);
                for(int k=0;k<num_cri;++k){
                    inFile>>weight;
                    max_w[k] = max(max_w[k],weight);
                    node_to_category[k][ID1-1].set(Weight_to_category(weight)+1, true);//record edge category Weight_to_category(weight)+1
                }
            }
            inFile.close();
            /// generate category index file
            ofstream outFile(file_r_graph1, ios::out);
            if (!outFile) {
                cout << "Write File opening failed." << endl;
                assert(outFile);
                exit(1);
            }
            cout<<"Writing index of node_to_category to disk...\t";
            outFile << num_cri << endl;
            for(int i=0;i<max_w.size();++i){
                if(i<max_w.size()-1){
                    outFile << max_w[i] << " ";
                }else{
                    outFile << max_w[i]<<endl;
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
        }else{//if it can be opened
            cout << "Loading preprocessed data of EM_Dijk..." << endl;
            inFile1 >> num_cri;
            max_w.assign(num_cri,0);
            assert(num_cri >= num_of_cri);
            //get the maximal edge weights
            for(int i=0;i<num_cri;++i){
                inFile1 >> temp_int;
                if(i<num_of_cri){
                    max_w[i]=temp_int;
                }
            }
            // get node-to-category mapping
            node_to_category.assign(num_of_cri,vector<bitset1<WeightPowMax+1,EdgeWeight>>(node_num));//bitset<WeightPowMax+1>
            for(int i=0;i<node_num;++i){
                for(int j=0;j<num_cri;++j){
                    inFile1 >> temp_int;
                    if(j<num_of_cri){
                        node_to_category[j][i].m_val = temp_int;
                    }

                }
            }
            inFile1.close();
        }

        tt.stop();
        cout << "Done. The time used for data preprocessing is: " << tt.GetRuntime() <<" s."<<endl;

        int temp_weightPowMax = 0;
        // get num_pool and pow_low
        for(int i=0;i<num_of_cri;++i){
            num_pool[i] = Weight_to_category(max_w[i])+1;
            if(temp_weightPowMax < num_pool[i])
                temp_weightPowMax = num_pool[i];
            for(int j=0;j<num_pool[i];++j){
                pow_low[i].push_back(pow(2,j));
            }
        }


        HotPool_EdgeSZ = MemGraph_EMDijk*1024*1024/(Edge_SZ*temp_weightPowMax);//update the hot pool edge size
        uint basic_mem = node_num*28/(1024*1024) + WeightPowMax*8 + 32;
        cout << "*** The actual basic memory consumption of EM_Dijk is: " << basic_mem << " MB. ***"<<  endl;
        cout << "*** The avoidable memory consumption caused by node_to_category is: " << node_num*16/(1024*1024) << " MB. ***" << endl;

    }
    //function of initialising multi-pass algorithms
    void EMMCGraph::EM_Dijk_Initiation(vector<HotPool>& hotPools){
        double nsize=0;
        hotPools.assign(num_pool[sc_i],HotPool());//reserve space for hot pools
        for(int i=0;i<hotPools.size();++i){
            hotPools[i].EM_Map.assign(partition_number,make_pair(0,0));
            if(sc_i<4){//for correlated criterion
                nsize = HotPool_EdgeSZ*4.0/(abs(i-10)+1);
                hotPools[i].set_capacity(nsize);
            }else{//for random criterion
                hotPools[i].set_capacity(HotPool_EdgeSZ);
            }
        }
    }
    //Function of Semi-external Memory Dijkstra: new version powered by GraphPool
    Distance EMMCGraph::EM_Dijk(int node_start, int node_end){//immediate pop
        if(node_start==node_end){
            cout<<"Same source and target vertex id!"<<endl;
            return 0;
        }
        //statistics of IO
        Stats stats_b, stats_e;
        stats_b.start();
        //Initiation
        vector<HotPool> hotPools;
        EM_Dijk_Initiation(hotPools);

        //Variables
        PriorityQueue pqueue(PQ_Pool);
        NodeId item_id,temp_id;
        Distance item_dis,temp_dis;
        VertexCost item_;      //recording top element of priority queue
        vector<EdgePairW> item_adj; //retrieve result of adjacency list
        bool flag_empty;        //flag of whether the priority queue will be empty if current top element is popped
        Distance relax_i_dis;   //the distance of the first element of relax i
        Distance nextTopDis;    //the distance of the next top element
        Distance min_cost;      //the shortest path distance
        vector<Distance> cost(node_num,INF);  //In-memory vertex cost vector
        //vector<NodeId> pre(node_num,-1);     //In-memory predecessor id
        VertexCost relax_top;  //the front element of relax_i
        vector<em_queue> relaxPools(num_pool[sc_i]);//em_queue()
        partitions_read.assign(partition_number, false);
        int cluster_id;
        vector<bool> clusterIM(partition_number,false); //flag of whether the partition is in memory

        Timer tt1;
        tt1.start();

        //Initiation of start node
        cost[node_start] = 0;
        pqueue.push(VertexCost(node_start,0));

        //Iteration
        while (!pqueue.empty()) {//for every node in pqueue
            //Step 1: top element in priority queue
            item_ = pqueue.top(); pqueue.pop();
            if(node_to_category[sc_i][item_.id].test(0)){//if already closed
                continue;
            }
            assert(item_.id<node_num);
            node_to_category[sc_i][item_.id].set(0,true);
            if (item_.id == node_end) {//if reach target node
                min_cost = cost[item_.id];
                ++visited_number;
                cout<<"Minimal cost of criterion " << sc_i << ": "<<min_cost<<endl;
                break;
            }
//            set_visited.insert(item_.id);
            ++visited_number;//update #visited_number
            //flag initiation for each round
            flag_empty = false;

            //eradicate the closed items in priority queue for next round
            while(!pqueue.empty()){
                if(node_to_category[sc_i][pqueue.top().id].test(0))
                    pqueue.pop();
                else
                    break;
            } //if closed

            //judge empty
            if(pqueue.empty()){//judge whether will priority queue be empty
                flag_empty = true;
                nextTopDis = INF;
            }else{
                nextTopDis = pqueue.top().cost;
            }

            //Step 2: Add information to relax pools
            for(int i=0;i<pow_low[sc_i].size();++i){
                if(node_to_category[sc_i][item_.id].test(i+1))
                    relaxPools[i].push(item_);
            }
            //Relax the nodes in relax pools
            for(int category_id=0;category_id<num_pool[sc_i];++category_id){
                if(!relaxPools[category_id].empty()){//if relax pool is not empty
                    //Step 3: check delayed relaxation condition
                    relax_i_dis = relaxPools[category_id].front().cost;//the oldest cost in relax_i
                    if(relax_i_dis + pow_low[sc_i][category_id] < nextTopDis || flag_empty){//if satisfy delayed relaxation condition
                        //Step 4: Relax the vertices in pool
                        while(!relaxPools[category_id].empty()){
                            relax_top = relaxPools[category_id].front();
                            //fetch the adjacency list of category i for vertex in category i
                            item_id = relax_top.id; item_dis = cost[item_id];
                            assert(item_id<node_num);
                            /// Relaxation
                            cluster_id = node_to_cluster[item_id];
                            if(!clusterIM[cluster_id]){//if the partition is in disk
                                if(ifOptimal){
                                    set_readCluster.insert(cluster_id);
                                }
                                ++p_num;
                                item_adj.clear();
                                EM_ReadCluster(category_id,item_id,item_adj,hotPools);
                                clusterIM[cluster_id] = true;
                                assert(!item_adj.empty());
                                //relaxation
                                for(auto it=item_adj.begin();it!=item_adj.end();++it){
                                    temp_id = it->ID2;
                                    assert(temp_id<node_num);
                                    if(node_to_category[sc_i][temp_id].test(0)){//if closed
                                        continue;
                                    }
                                    temp_dis = item_dis + it->w;
                                    if (cost[temp_id] > temp_dis) {//slack operation
                                        cost[temp_id] = temp_dis;
                                        //pre[temp_id] = item_id;
                                        pqueue.push(VertexCost(temp_id, temp_dis));
                                    }
                                }
                            }else{//if the partition is in memory or external vector
                                // check memory
                                if(hotPools[category_id].Edges_IM.find(item_id)!=hotPools[category_id].Edges_IM.end()){//if found in memory
                                    auto iter = hotPools[category_id].Edges_IM.equal_range(item_id);
                                    //relaxation
                                    for(auto it=iter.first;it!=iter.second;++it){
                                        temp_id = it->second.ID2;
                                        assert(temp_id<node_num);
                                        if(node_to_category[sc_i][temp_id].test(0)){//if closed
                                            continue;
                                        }
                                        temp_dis = item_dis + it->second.w;
                                        if (cost[temp_id] > temp_dis) {//slack operation
                                            cost[temp_id] = temp_dis;
                                            //pre[temp_id] = item_id;
                                            pqueue.push(VertexCost(temp_id, temp_dis));
                                        }
                                    }
                                    hotPools[category_id].Edges_IM.erase(item_id);//delete visited edges
                                }else if(hotPools[category_id].EM_Map[cluster_id].second > 0) {//check the external vector
                                    bool flag_finish = false;
                                    // linear search
                                    for (int i = hotPools[category_id].EM_Map[cluster_id].first; i < hotPools[category_id].EM_Map[cluster_id].second; ++i) {
                                        if (hotPools[category_id].Edges_EM[i].ID1 == item_id) {//if it is the target edge
                                            //relaxation
                                            for (auto it = item_adj.begin(); it != item_adj.end(); ++it) {
                                                temp_id = hotPools[category_id].Edges_EM[i].ID2;
                                                assert(temp_id < node_num);
                                                if (node_to_category[sc_i][temp_id].test(0)) {//if closed
                                                    continue;
                                                }
                                                temp_dis = item_dis + hotPools[category_id].Edges_EM[i].w;
                                                if (cost[temp_id] > temp_dis) {//slack operation
                                                    cost[temp_id] = temp_dis;
                                                    //pre[temp_id] = item_id;
                                                    pqueue.push(VertexCost(temp_id, temp_dis));
                                                }
                                            }
                                            flag_finish = true;
                                        } else if (flag_finish) {
                                            break;
                                        }
                                    }
                                    assert(flag_finish);
                                }else{
                                    cout << "Wrong!!!" << endl;
                                    assert(!clusterIM[cluster_id]);
                                }
                            }
                            relaxPools[category_id].pop();
                        }
                    }
                }
            }
        }
        tt1.stop();
        query_time += tt1.GetRuntime() * 1000;//
        while(!pqueue.empty()){
            pqueue.pop();
        }
        //io record
        stats_e = stats_b.get_stats();
        if (stats_e.m_noOfWrites > 0 || stats_e.m_noOfReads > 0) {
            //cout << "There are " << stats_e.m_noOfWrites + stats_e.m_noOfReads << " io for stxxl." << endl;
            //cout<<stats_e<<endl;
            EM_IORecord(stats_e);
        }
        //clear
        for(int i=0;i<num_pool[sc_i];++i){
            while(!relaxPools[i].empty()){
                relaxPools[i].pop();
            }
            hotPools[i].clear();
        }
        partitions_read.clear();
        return min_cost;
    }
    //function of reading partition from disk
    void EMMCGraph::EM_ReadCluster(int category_id,NodeId ID1,vector<EdgePairW>& resultVE, vector<HotPool>& hotPools){
        Timer tt;
        tt.start();
        unordered_map<int,vector<EdgePairW>> tempEdges;//from category to Edges
        int target_p_id = node_to_cluster[ID1];
        int partition_id;
        int temp_id;
        int temp_category;
        EdgePairW temp_edge;

        if(!partitions_read[target_p_id]){//if the partition has not been read
            partitions_read[target_p_id] = true;
        }else{//if true
            cout<<"!!!Wrong! The partition has already been read!"<<endl;
            assert(!partitions_read[target_p_id]);
            exit(-1);
        }
        //read partition
        char filePath[100];
        strcpy(filePath, DataPath);
        strcat(filePath, dataset.c_str());
        strcat(filePath, "/Partitions/");
        strcat(filePath, dataset.c_str());
        strcat(filePath, "_Partition_");
        strcat(filePath, to_string(PartitionSize).c_str());
        strcat(filePath, "_");
        strcat(filePath, to_string(target_p_id).c_str());
        strcat(filePath, ".bin");
        FILE* file = fopen(filePath, "rb");

        ReadBuffer read_io(file);

        int u, u_deg;
        int v, weight;
        //read partition
        while(!read_io.isend) {//if not read out
            read_io.read(&u);//get node id
            temp_id = u-1;
            assert(u<=node_num);
            read_io.read(&u_deg);//get node degree
            tempEdges.clear();
            //vertex store in indexNode, edge store in e
            for (int i = 0; i < u_deg; ++i) {
                read_io.read(&v);//read end node
                assert(v<=node_num);
                for (int j = 0; j < NUM_OF_CRITERIA; ++j) {
                    read_io.read(&weight);
                    if (j == sc_i){
//                        tempEdges[Weight_to_category(weight)].emplace_back(EdgePairW(u-1, v-1, weight));
                        temp_category = Weight_to_category(weight);
                        temp_edge.ID1 = temp_id; temp_edge.ID2 = v-1; temp_edge.w = weight;
                        if(temp_category==category_id && temp_id==ID1){//if it is the target
                            resultVE.emplace_back(temp_edge);
                        }else{//if it is not the target, insert into memory or external vector
//                            tempEdges[temp_category].emplace_back(temp_edge);
                            if (hotPools[temp_category].Edges_IM.size() + u_deg < hotPools[temp_category].capacityEdges) {//if not exceed the memory,push into in-memory map
                                hotPools[temp_category].Edges_IM.insert({temp_id, temp_edge});
                            } else {//if memory is full
                                //push into external vector
                                partition_id = node_to_cluster[temp_id];
                                if (hotPools[temp_category].EM_Map[partition_id].second == 0) {//if not exist
                                    hotPools[temp_category].EM_Map[partition_id].first = hotPools[temp_category].Edges_EM.size();
                                }
                                hotPools[temp_category].Edges_EM.push_back(temp_edge);
                                hotPools[temp_category].EM_Map[partition_id].second = hotPools[temp_category].Edges_EM.size();
                            }
                        }
                    }
                }
            }
        }
        io_num += read_io.io_number;//update #IO
        fclose(file);
        tt.stop();
        io_time += tt.GetRuntime() * 1000;
        if(cluster_to_IO[target_p_id] == 0){
            cluster_to_IO[target_p_id] = read_io.io_number;
        }
    }
    //function of partition-based Dijkstra's algorithm
    /*Distance EMMCGraph::Dijkstra_IO(int node_start, int node_end){
        if(node_start==node_end){
            cout<<"Same source and target vertex id!"<<endl;
            return 0;
        }
        //statstics of IO
        Stats stats_b, stats_e;
        stats_b.start();
        //Variables
        PriorityQueue pqueue(PQ_Pool);
        HotPool2 myHotPool(partition_number,PageNumber_DijkIO);

        VertexCost item_;//recording top element of priority queue
        NodeId item_id, temp_id;
        Distance item_dis, temp_dis;
        vector<bool> closed(node_num, false); //flag vector of whether closed
        vector<Distance> cost(node_num, INF);   //vector of cost
        vector<NodeId> pre(node_num, -1);       //vector of predecessor id
        Distance min_cost = INF;
        int index_begin, index_end;
        int partition_id;

        Timer tt1;
        tt1.start();
        //Initiation of start node
        cost[node_start] = 0;//cost of start node
        pqueue.push(VertexCost(node_start,0));

        //Iteration
        while (!pqueue.empty()) {//for every node in pqueue
            item_ = pqueue.top(); pqueue.pop();
            if(closed[item_.id]){//if already closed
                continue;
            }
            item_id = item_.id; item_dis = item_.cost;
            if (item_id == node_end) {//if reach target node
                min_cost = cost[item_id];
//                set_visited.insert(item_id);
                ++visited_number;
//                cout<<"Minimal cost of criterion " << sc_i << ": "<<min_cost<<endl;
                break;
            }
            //relaxation
//            set_visited.insert(item_id);//update the visited set for this criteria
            ++visited_number;//update #visited_number

            partition_id = node_to_cluster[item_id];

            if(!myHotPool.FlagIM[partition_id]){//if not in memory
                GraphReadCluster(myHotPool,partition_id);
                ++p_num;
            }

            assert(get<1>(EMEdgesIndex[item_id]) != -1);
            /// read from EMMCEdges
            for(int i = get<1>(EMEdgesIndex[item_id]); i < get<2>(EMEdgesIndex[item_id]); ++i){
                temp_id = myHotPool.HotPools[partition_id][i].ID2;
                if (closed[temp_id])//if closed
                    continue;
                temp_dis = item_dis + myHotPool.HotPools[partition_id][i].w;
                if (cost[temp_id] > temp_dis) {//slack operation
                    cost[temp_id] = temp_dis;
                    pre[temp_id] = item_id;
                    pqueue.push(VertexCost(temp_id, temp_dis));
                }
            }

            closed[item_id] = true;
        }
//        Dij_getPath(pre,node_start,node_end);//traverse the pre vector to get the shortest path
        //io record
        stats_e = stats_b.get_stats();
        if (stats_e.m_noOfWrites > 0 || stats_e.m_noOfReads > 0) {
//            cout << "There are " << stats_e.m_noOfWrites + stats_e.m_noOfReads << " io for stxxl." << endl;
            //cout<<stats_e<<endl;
            EM_IORecord(stats_e);
        }
        tt1.stop();
        query_time += tt1.GetRuntime() * 1000;//
        while(!pqueue.empty())
            pqueue.pop();
        return min_cost;
    }*/
    //function of stxxl vector based Dijkstra's algorithm
    Distance EMMCGraph::Dijkstra_IO(int node_start, int node_end, VectorMCEdgesEMTuple_DijkIO & EMMCEdges){
        if(node_start==node_end){
            cout<<"Same source and target vertex id!"<<endl;
            return 0;
        }
        //statstics of IO
        Stats stats_b, stats_e;
        stats_b.start();
        //Variables
        PriorityQueue pqueue(PQ_Pool);

        VertexCost item_;//recording top element of priority queue
        NodeId item_id, temp_id;
        Distance item_dis, temp_dis;
        vector<bool> closed(node_num, false); //flag vector of whether closed
        vector<Distance> cost(node_num, INF);   //vector of cost
        //vector<NodeId> pre(node_num, -1);       //vector of predecessor id
        Distance min_cost = INF;
        int index_begin, index_end;

        Timer tt1;
        tt1.start();
        //Initiation of start node
        cost[node_start] = 0;//cost of start node
        pqueue.push(VertexCost(node_start,0));

        //Iteration
        while (!pqueue.empty()) {//for every node in pqueue
            item_ = pqueue.top(); pqueue.pop();
            if(closed[item_.id]){//if already closed
                continue;
            }
            item_id = item_.id; item_dis = item_.cost;
            if (item_id == node_end) {//if reach target node
                min_cost = cost[item_id];
//                set_visited.insert(item_id);
                ++visited_number;
                cout<<"Minimal cost of criterion " << sc_i << ": "<<min_cost<<endl;
                break;
            }
            //relaxation
//            set_visited.insert(item_id);//update the visited set for this criteria
            ++visited_number;//update #visited_number
            if(get<1>(EMEdgesIndex[item_id]) == -1){
                ReadClusterToExternalVectorMC(node_to_cluster[item_id],EMMCEdges);
                ++p_num;
                assert(get<1>(EMEdgesIndex[item_id]) != -1);
            }
            /// read from EMMCEdges
            for(int i = get<1>(EMEdgesIndex[item_id]); i < get<2>(EMEdgesIndex[item_id]); ++i){
                assert(EMMCEdges[i].ID1 == item_id);
                temp_id = EMMCEdges[i].ID2;
                if (closed[temp_id])//if closed
                    continue;
                temp_dis = item_dis + EMMCEdges[i].getW(sc_i);
                if (cost[temp_id] > temp_dis) {//slack operation
                    cost[temp_id] = temp_dis;
                    //pre[temp_id] = item_id;
                    pqueue.push(VertexCost(temp_id, temp_dis));
                }
            }

            closed[item_id] = true;
        }
//        Dij_getPath(pre,node_start,node_end);//traverse the pre vector to get the shortest path
        //io record
        stats_e = stats_b.get_stats();
        if (stats_e.m_noOfWrites > 0 || stats_e.m_noOfReads > 0) {
//            cout << "There are " << stats_e.m_noOfWrites + stats_e.m_noOfReads << " io for stxxl." << endl;
            //cout<<stats_e<<endl;
            EM_IORecord(stats_e);
        }
        tt1.stop();
        query_time += tt1.GetRuntime() * 1000;//
        while(!pqueue.empty())
            pqueue.pop();
        return min_cost;
    }
    //function of reading partition data to external vector
    template<class T>
    void EMMCGraph::ReadClusterToExternalVectorMC(int target_p_id, T & EMMCEdges){
        char filePath[100];
        int u, u_deg;
        int v, weight;
        bool flag_sizeUpdate = false;
        int temp_id = 0;
        MCEdgeT temp_edge;

        strcpy(filePath, DataPath);
        strcat(filePath, dataset.c_str());
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

        //read partition
        while(!read_io.isend) {//if not read out
            read_io.read(&u);//get node id

            if(u != temp_id){
                get<1>(EMEdgesIndex[u-1]) = EMMCEdges.size();
                if(temp_id != 0){
                    get<2>(EMEdgesIndex[temp_id-1]) = EMMCEdges.size();
                }
                temp_id = u;
            }
            temp_edge.ID1 = u-1;
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
                EMMCEdges.push_back(temp_edge);
            }
        }
        get<2>(EMEdgesIndex[u-1]) = EMMCEdges.size();
        io_num += read_io.io_number;//update #IO
        fclose(file);
        tt.stop();
        io_time += tt.GetRuntime() * 1000;
        if(cluster_to_IO[target_p_id] == 0){
            cluster_to_IO[target_p_id] = read_io.io_number;
        }
    }
    //function of reading partition data into hot pool
    void EMMCGraph::GraphReadCluster(HotPool2 & mcHotPool,int target_p_id){
        int evict_p_id = mcHotPool.LRU_IMCluster.put(target_p_id,1);
        //// evict old partition for new one
        if(evict_p_id != -1) {//if it is necessary to evict old partition
            mcHotPool.FlagIM[evict_p_id] = false;
            mcHotPool.HotPools[evict_p_id].clear();
            mcHotPool.HotPools[evict_p_id].shrink_to_fit();
        }
        //// read partition from original disk file to internal vector
        //check status of storage
        char filePath[100];
        int u, u_deg;
        int v, weight;
        MCEdgeT temp_edge;

        strcpy(filePath, DataPath);
        strcat(filePath, dataset.c_str());
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

            read_io.read(&u_deg);//get node degree
            //vertex store in indexNode, edge store in e
            for (int i = 0; i < u_deg; ++i) {
                read_io.read(&v);//read end node
                temp_edge.ID2 = v-1;
                for (int j = 0; j < NUM_OF_CRITERIA; ++j) {
                    read_io.read(&weight);
                    if(j == sc_i){
                        mcHotPool.HotPools[target_p_id].emplace_back(EdgePairW(u-1,v-1,weight));
                    }
                }
            }
        }

        get<2>(EMEdgesIndex[u-1]) = mcHotPool.HotPools[target_p_id].size();

        io_num += read_io.io_number;//update #IO
        fclose(file);
        tt.stop();
        io_time += tt.GetRuntime() * 1000;
        cluster_to_IO[target_p_id] = read_io.io_number;

        mcHotPool.FlagIM[target_p_id] = true;
    }
    //function of stxxl vector based BiDijkstra's algorithm
    Distance EMMCGraph::BiDijkstra_IO(int node_start, int node_end, VectorMCEdgesEMTuple_BiDijkIO & EMMCEdges){
        if(node_start==node_end){
            cout<<"Same source and target vertex id!"<<endl;
            return 0;
        }
        //statstics of IO
        Stats stats_b, stats_e;
        stats_b.start();
        //Variables
        PriorityQueue pqueue(PQ_Pool);
        PriorityQueue pqueue_r(PQ_Pool);

        NodeId item_id, temp_id;
        NodeId item_id_r, temp_id_r;
        Distance item_dis, temp_dis;
        Distance item_dis_r, temp_dis_r;
        vector<bool> closed(node_num, false); //flag vector of whether closed
        vector<Distance> cost(node_num, INF);   //vector of cost
        //vector<NodeId> pre(node_num, -1);       //vector of predecessor id
        vector<bool> closed_r(node_num, false); //flag vector of whether closed
        vector<Distance> cost_r(node_num, INF);   //vector of cost
        //vector<NodeId> pre_r(node_num, -1);       //vector of predecessor id
        Distance min_cost = INF;
        int index_begin, index_end;
        int terminate_id = -1;

        Timer tt1;
        tt1.start();
        //Initiation of start node
        cost[node_start] = 0;//cost of start node
        pqueue.push(VertexCost(node_start,0));
        cost_r[node_end] = 0;//cost of start node
        pqueue_r.push(VertexCost(node_end,0));

        //Iteration
        while (!pqueue.empty() || !pqueue_r.empty()) {//for every node in pqueue

            //termination judgement
            if (pqueue.top().cost + pqueue_r.top().cost >= min_cost) {//condition of termination
                cout << node_start << " " << node_end << " " << sc_i << " "<<min_cost << endl;
                break;
            }

            /// Forward searching
            item_id = pqueue.top().id;
            item_dis = pqueue.top().cost;
            pqueue.pop();
            closed[item_id] = true;

            //relaxation
//            set_visited.insert(item_id);//update the visited set for this criteria
            ++visited_number;//update #visited_number
            if(get<1>(EMEdgesIndex[item_id]) == -1){
                ReadClusterToExternalVectorMC(node_to_cluster[item_id],EMMCEdges);
                ++p_num;
                assert(get<1>(EMEdgesIndex[item_id]) != -1);
            }
            // read from EMMCEdges
            for(int i = get<1>(EMEdgesIndex[item_id]); i < get<2>(EMEdgesIndex[item_id]); ++i){
                assert(EMMCEdges[i].ID1 == item_id);
                temp_id = EMMCEdges[i].ID2;
                temp_dis = item_dis + EMMCEdges[i].getW(sc_i);
                if (!closed[temp_id]){//if not closed
                    if (cost[temp_id] > temp_dis) {//slack operation
                        cost[temp_id] = temp_dis;
                        //pre[temp_id] = item_id;
                        pqueue.push(VertexCost(temp_id, temp_dis));
                    }
                }
                if (closed_r[temp_id] && temp_dis + cost_r[temp_id] < min_cost) {
                    min_cost = temp_dis + cost_r[temp_id];
                    terminate_id = temp_id;
                }

            }
            /// Reverse searching
            item_id_r = pqueue_r.top().id;
            item_dis_r = pqueue_r.top().cost;
            pqueue_r.pop();
            closed_r[item_id_r] = true;

            //relaxation
//            set_visited.insert(item_id);//update the visited set for this criteria
            ++visited_number;//update #visited_number
            if(get<1>(EMEdgesIndex[item_id_r]) == -1){
                ReadClusterToExternalVectorMC(node_to_cluster[item_id_r],EMMCEdges);
                ++p_num;
                assert(get<1>(EMEdgesIndex[item_id_r]) != -1);
            }
            // read from EMMCEdges
            for(int i = get<1>(EMEdgesIndex[item_id_r]); i < get<2>(EMEdgesIndex[item_id_r]); ++i){
                assert(EMMCEdges[i].ID1 == item_id_r);
                temp_id_r = EMMCEdges[i].ID2;
                temp_dis_r = item_dis_r + EMMCEdges[i].getW(sc_i);
                if (!closed_r[temp_id_r]){//if not closed
                    if (cost_r[temp_id_r] > temp_dis_r) {//slack operation
                        cost_r[temp_id_r] = temp_dis_r;
                        //pre_r[temp_id_r] = item_id_r;
                        pqueue_r.push(VertexCost(temp_id_r, temp_dis_r));
                    }
                }
                if (closed[temp_id_r] && temp_dis_r + cost[temp_id_r] < min_cost) {
                    min_cost = temp_dis_r + cost[temp_id_r];
                    terminate_id = temp_id_r;
                }
            }
            //deal with invalid top elements
            while (closed[pqueue.top().id]) {//if already closed
                pqueue.pop();
            }
            while (closed_r[pqueue_r.top().id]) {
                pqueue_r.pop();
            }
        }
//        Dij_getPath(pre,node_start,node_end);//traverse the pre vector to get the shortest path
        //io record
        stats_e = stats_b.get_stats();
        if (stats_e.m_noOfWrites > 0 || stats_e.m_noOfReads > 0) {
//            cout << "There are " << stats_e.m_noOfWrites + stats_e.m_noOfReads << " io for stxxl." << endl;
            //cout<<stats_e<<endl;
            EM_IORecord(stats_e);
        }
        tt1.stop();
        query_time += tt1.GetRuntime() * 1000;//
        while(!pqueue.empty())
            pqueue.pop();
        while(!pqueue_r.empty())
            pqueue_r.pop();
        return min_cost;
    }
    //Entry for one-pass algorithms
    double EMMCGraph::MC_OnePass(const string& qtype,int algo_choice){
        string LongDis,MediumDis,ShortDis,RandomDis;
        LongDis = string(DataPath) + dataset + "/" + dataset + "_OD_LongDis_" + to_string(PartitionSize) + ".txt";
        MediumDis = string(DataPath) + dataset + "/" + dataset + "_OD_MediumDis_"+ to_string(PartitionSize) +".txt";
        ShortDis = string(DataPath) + dataset + "/" + dataset + "_OD_ShortDis_"+ to_string(PartitionSize) +".txt";
        RandomDis = string(DataPath) +  dataset + "/" + dataset + "_OD_Random_"+ to_string(PartitionSize) +".txt";

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
            QueryType = "S";
            total_time += MC_OnePass_one(ShortDis, "SDis", algo_choice);//Efficiency evaluation on short distance
            QueryType = "M";
            total_time += MC_OnePass_one(MediumDis, "MDis", algo_choice);//Efficiency evaluation on medium distance
            QueryType = "L";
            total_time += MC_OnePass_one(LongDis, "LDis", algo_choice);//Efficiency evaluation on long distance
        }

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
                basic_mem = node_num*(num_of_cri*5 + 11)/(1024*1024) + 32*num_of_cri;
                cout << "*** The basic memory consumption of OneHopNoIO is: " << basic_mem << " MB. The graph data memory is "<< MemGraph_NoIO<<" MB. ***"<<  endl; break;
            }
            case MultiHopsNoIO:{
                basic_mem = node_num*(num_of_cri*5 + 11)/(1024*1024) + 32*num_of_cri;
                cout << "*** The basic memory consumption of MultiHopsNoIO is: " << basic_mem << " MB. The graph data memory is "<< MemGraph_NoIO<<" MB. ***"<<  endl; break;
            }
            case BiMultiHopsNoIO:{
                basic_mem = node_num*(num_of_cri*10 + 11)/(1024*1024) + 64*num_of_cri;
                cout << "*** The basic memory consumption of BiMultiHopsNoIO is: " << basic_mem << " MB. The graph data memory is "<< MemGraph_NoIO_Bi<<" MB. ***"<<  endl; break;
            }
            case OneHop:{
                basic_mem = node_num*(num_of_cri*5 + 12)/(1024*1024) + 32*num_of_cri;
                cout << "*** The basic memory consumption of OneHop is: " << basic_mem << " MB. The graph data memory is "<< MemGraph_IO<<" MB. ***"<< endl;
                cout << "Alpha: " << alpha << "\tMu: " << double(MuForEM)/10 << endl;
                break;
            }
            case MultiHops:{
                basic_mem = node_num*(num_of_cri*5 + 12)/(1024*1024) + 32*num_of_cri;
                cout << "*** The basic memory consumption of MultiHops is: " << basic_mem << " MB. The graph data memory is "<< MemGraph_IO<<" MB. ***"<<  endl;
                cout << "Alpha: " << alpha_multi <<"\tMu: " << double(MuForEM) / 10 << endl;
                break;
            }
            case BiMultiHops:{
                basic_mem = node_num*(num_of_cri*10 + 12)/(1024*1024) + 64*num_of_cri;
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
        gbxxl::PriorityQueue pqueue2(gbxxl::PQ_Pool); gbxxl::PriorityQueue pqueue3(gbxxl::PQ_Pool);
        gbxxl::PriorityQueue pqueue4(gbxxl::PQ_Pool);
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
            gbxxl::PriorityQueue pqueue0_r(gbxxl::PQ_Pool); gbxxl::PriorityQueue pqueue1_r(gbxxl::PQ_Pool);
            gbxxl::PriorityQueue pqueue2_r(gbxxl::PQ_Pool); gbxxl::PriorityQueue pqueue3_r(gbxxl::PQ_Pool);
            gbxxl::PriorityQueue pqueue4_r(gbxxl::PQ_Pool);
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
            for (int i = 0; i < run_times; ++i) {
                ID1 = ODpair[i].first;
                ID2 = ODpair[i].second;
//                cout << "Query "<<i<<" : "<<ID1-1 << " " << ID2-1 << endl;
                CommonInitiation_IO();//common initiation of each round
                OnePassInitiation_IO(algo_choice);

                if(algo_choice == BiMultiHops){
                    EM_MC_BiMultiHop(ID1-1,ID2-1);
                }else if(algo_choice == BiMultiHopsNoIO){
                    EM_MC_BiMultiHop_NoIO(ID1-1,ID2-1);
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
//                cout << "Query: "<<ID1-1 << " " << ID2-1 << endl;

                CommonInitiation_IO();//common initiation of each round
                OnePassInitiation_IO(algo_choice);
                switch (algo_choice){
                    case OneHopNoIO:{
                        EM_MC_OneHop_NoIO(ID1-1, ID2-1); break;
                    }
                    case MultiHopsNoIO:{
                        EM_MC_MultiHop_NoIO(ID1-1, ID2-1); break;
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
//                cout << "I/O number: " << io_num << endl;
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
        cout << "! Average share visited nodes: " << ave_share << endl;
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
    //One-hop algorithm without io optimization
    void EMMCGraph::EM_MC_OneHop_NoIO(int node_start, int node_end){
        if (node_start == node_end) {
            cout << "Source vertex and Target vertex are the same!!!" << endl;
            return;
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
        VectorMCEdgesEMTuple_NoIO EMMCEdges;//the external vectors of different criteria
        vector<vector<Distance>> mc_cost(num_of_cri, vector<Distance>(node_num, INF));   //vector of cost
        //vector<vector<NodeId>> mc_pre(num_of_cri, vector<int>(node_num, -1));      //vector of predecessor id
        vector<vector<bool>> mc_closed(num_of_cri, vector<bool>(node_num));    //flag of vertex closing
        vector<int> cluster_to_cri(partition_number,0);//used to map which criterion read this partition
        unordered_set<int> set_remain;//set for storing the id of unshared criteria

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
                while (mc_closed[cri_i][EM_MC_PQueue[cri_i]->top().id]) {//if found in visited set
                    EM_MC_PQueue[cri_i]->pop();
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
                    //relaxation
                    partition_id = node_to_cluster[item_id];
                    if (get<1>(EMEdgesIndex[item_id]) == -1) {
                        ReadClusterToExternalVectorMC(partition_id, EMMCEdges);
                        ++p_num;
                        assert(get<1>(EMEdgesIndex[item_id]) != -1);
                    }
                    cluster_to_cri[partition_id] = cri_i;//record which criterion leads to reading from external disk
                    if (cluster_to_cri[partition_id] != cri_i) {//if the cluster is not read by this criterion
                        ++share_number;
                    }
                    //set closed
                    mc_closed[cri_i][item_id] = true;
                    ++visited_number;//update #visited_number
                    //slack
                    /// read from EMMCEdges
                    for (int i = get<1>(EMEdgesIndex[item_id]); i < get<2>(EMEdgesIndex[item_id]); ++i) {
                        assert(i < EMMCEdges.size());
                        assert(EMMCEdges[i].ID1 == item_id);
                        temp_id = EMMCEdges[i].ID2;
                        assert(temp_id < node_num);
                        if (mc_closed[cri_i][temp_id])//if closed
                            continue;
                        temp_dis = item_dis + EMMCEdges[i].getW(cri_i);
                        if (mc_cost[cri_i][temp_id] > temp_dis) {//slack operation
                            mc_cost[cri_i][temp_id] = temp_dis;
                            //mc_pre[cri_i][temp_id] = item_id;
                            EM_MC_PQueue[cri_i]->push(VertexCost(temp_id, temp_dis));
                        }
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
        //clear
        EMMCEdges.clear();
    }
    //Multi-hops algorithm without io optimization
    void EMMCGraph::EM_MC_MultiHop_NoIO(int node_start, int node_end){
        if (node_start == node_end) {
            cout << "Source vertex and Target vertex are the same!!!" << endl;
            return;
        }
        //statstics of IO
        Stats stats_b, stats_e;
        stats_b.start();
        //variables
        VectorMCEdgesEMTuple_NoIO EMMCEdges;
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
        unordered_set<int> set_remain;//set for storing the id of unshared criteria

        Timer tt;
        tt.start();

        //Initiation of start node
        for (int i = 0; i < num_of_cri; ++i) {
            set_cri.insert(i);
            mc_cost[i][node_start] = 0;//cost of start node
            EM_MC_PQueue[i]->push(VertexCost(node_start,0));
        }
        //first read
        ReadClusterToExternalVectorMC(node_to_cluster[node_start],EMMCEdges);
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
                }
                item_id = EM_MC_PQueue[cri_i]->top().id;// top min item
                item_dis = EM_MC_PQueue[cri_i]->top().cost;

                partition_id = node_to_cluster[item_id];

                while (get<1>(EMEdgesIndex[item_id]) != -1) {//if found
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
                        //slack
                        /// read from EMMCEdges
                        for (int i = get<1>(EMEdgesIndex[item_id]); i < get<2>(EMEdgesIndex[item_id]); ++i) {
                            assert(i < EMMCEdges.size());
                            assert(EMMCEdges[i].ID1 == item_id);
                            temp_id = EMMCEdges[i].ID2;
                            assert(temp_id < node_num);
                            if (mc_closed[cri_i][temp_id])//if closed
                                continue;
                            temp_dis = item_dis + EMMCEdges[i].getW(cri_i);
                            if (mc_cost[cri_i][temp_id] > temp_dis) {//slack operation
                                mc_cost[cri_i][temp_id] = temp_dis;
                                //mc_pre[cri_i][temp_id] = item_id;
                                EM_MC_PQueue[cri_i]->push(VertexCost(temp_id, temp_dis));
                            }
                        }


                        while (mc_closed[cri_i][EM_MC_PQueue[cri_i]->top().id]) {//if already closed
                            EM_MC_PQueue[cri_i]->pop();
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
                    ReadClusterToExternalVectorMC(temp_cluster_id, EMMCEdges);
                    ++p_num;
                }
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
        EMMCEdges.clear();
    }
    //BiMulti-hops algorithm without io optimization: 2022-03-23 version 1
    void EMMCGraph::EM_MC_BiMultiHop_NoIO(int node_start, int node_end) {
        if (node_start == node_end) {
            cout << "Source vertex and Target vertex are the same!!!" << endl;
            return;
        }
        //statstics of IO
        Stats stats_b, stats_e;
        stats_b.start();
        //variables
        unordered_set<int> set_remain;//set for storing the id of unshared criteria
        VectorMCEdgesEMTuple_NoIO_Bi EMMCEdges;
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
        ReadClusterToExternalVectorMC(node_to_cluster[node_start],EMMCEdges);
        ++p_num;
        cluster_to_cri[node_to_cluster[node_start]] = 0;
        cluster_to_bi[node_to_cluster[node_start]] = FORWARD;
        if(node_to_cluster[node_end]!=node_to_cluster[node_start]){
            ReadClusterToExternalVectorMC(node_to_cluster[node_end],EMMCEdges);
            ++p_num;
            cluster_to_cri[node_to_cluster[node_end]] = 0;
            cluster_to_bi[node_to_cluster[node_end]] = REVERSE;
        }

        //Iteration
        while (!EM_JudgeEmptyBi(EM_MC_PQueue, EM_MC_PQueue_r)) {//for every node in pqueue !mc_finished[0] || !mc_finished[1] || !mc_finished[2]
            set_remain = set_cri;
            //EM_BiMultiHop_strategy_New(mc_HotPool, set_remain, EM_MC_PQueue, EM_MC_PQueue_r);//divide_and_slack
            while (!set_remain.empty()) {
                cri_i = *set_remain.begin();//pick one criterion to process
                //deal with invalid top elements
                while (mc_closed[cri_i][EM_MC_PQueue[cri_i]->top().id]) {
                    EM_MC_PQueue[cri_i]->pop();
                }
                while (mc_closed_r[cri_i][EM_MC_PQueue_r[cri_i]->top().id]) {
                    EM_MC_PQueue_r[cri_i]->pop();
                }
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
                else {//if not terminate
                    //// Forward searching
                    item_id = EM_MC_PQueue[cri_i]->top().id;// top min item
                    item_dis = EM_MC_PQueue[cri_i]->top().cost;// top min item

                    partition_id = node_to_cluster[item_id];
                    while (get<1>(EMEdgesIndex[item_id]) != -1 && !mc_closed_r[cri_i][item_id]) {//if found
                        ++cri_hops[cri_i];
                        if (cluster_to_cri[node_to_cluster[item_id]] != cri_i) {//if the cluster is not read by this criterion
                            ++share_number;
                        }
                        //set closed
                        mc_closed[cri_i][item_id] = true;
                        EM_MC_PQueue[cri_i]->pop();//pop min item
                        ++visited_number;//update #visited_number
                        if(item_id == node_end){
                            mc_finished[cri_i] = true;
                            set_cri.erase(cri_i);
                            mc_min_cost[cri_i] = mc_cost[cri_i][item_id];
//                        cout<<"Minimal cost of criterion " << cri_i << ": "<<mc_min_cost[cri_i]<<endl;
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
                        /// read from EMMCEdges
                        assert(get<2>(EMEdgesIndex[item_id]) > 0);
                        for (int i = get<1>(EMEdgesIndex[item_id]); i < get<2>(EMEdgesIndex[item_id]); ++i) {
                            assert(i < EMMCEdges.size());
                            assert(EMMCEdges[i].ID1 == item_id);
                            temp_id = EMMCEdges[i].ID2;
                            temp_dis = item_dis + EMMCEdges[i].getW(cri_i);
                            assert(temp_id < node_num);
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


                        while (mc_closed[cri_i][EM_MC_PQueue[cri_i]->top().id]) {//if already closed
                            EM_MC_PQueue[cri_i]->pop();
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
                    while (get<1>(EMEdgesIndex[item_id_r]) != -1 && !mc_closed[cri_i][item_id_r]) {//if found
                        ++cri_hops_r[cri_i];
                        if (cluster_to_cri[node_to_cluster[item_id_r]] != cri_i) {//if the cluster is not read by this criterion
                            ++share_number;
                        }
                        //set closed
                        mc_closed_r[cri_i][item_id_r] = true;
                        EM_MC_PQueue_r[cri_i]->pop();//pop min item
                        ++visited_number;//update #visited_number
                        //slack
                        /// read from EMMCEdges
                        assert(get<2>(EMEdgesIndex[item_id_r]) > 0);
                        for (int i = get<1>(EMEdgesIndex[item_id_r]); i < get<2>(EMEdgesIndex[item_id_r]); ++i) {
                            assert(i < EMMCEdges.size());
                            assert(EMMCEdges[i].ID1 == item_id_r);
                            temp_id_r = EMMCEdges[i].ID2;
                            temp_dis_r = item_dis_r + EMMCEdges[i].getW(cri_i);
                            assert(temp_id_r < node_num);
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

                        while (mc_closed_r[cri_i][EM_MC_PQueue_r[cri_i]->top().id]) {
                            EM_MC_PQueue_r[cri_i]->pop();
                        }
                        item_id_r = EM_MC_PQueue_r[cri_i]->top().id;// top min item
                        item_dis_r = EM_MC_PQueue_r[cri_i]->top().cost;// top min item

                        partition_id_r = node_to_cluster[item_id_r];
                    }
                    cri_to_cluster_r[cri_i] = node_to_cluster[item_id_r];
                    cluster_to_cri[cri_to_cluster_r[cri_i]] = cri_i;//node_to_cluster[item_id_r]
                    cri_pqueue_r.update(cri_i, cri_hops_r[cri_i]);
                }
                set_remain.erase(cri_i);
                // judge reading
                if (set_remain.empty() && !flag_mcFinished) {//if all criteria are processed and not all criteria are finished
                    /// read by synchronizing the hops
                    temp_cluster_id = cri_to_cluster[cri_pqueue.top_value()];
                    if (cluster_to_bi[temp_cluster_id] == 0 || get<1>(EMEdgesIndex[partitions_info[temp_cluster_id]]) == -1) {
                        ReadClusterToExternalVectorMC(temp_cluster_id, EMMCEdges);
                        cluster_to_bi[temp_cluster_id] = FORWARD;
                        ++p_num;
                    }
                    temp_cluster_id_r = cri_to_cluster_r[cri_pqueue_r.top_value()];
                    if (cluster_to_bi[temp_cluster_id_r] == 0 || get<1>(EMEdgesIndex[partitions_info[temp_cluster_id_r]]) == -1) {
                        if (temp_cluster_id != temp_cluster_id_r) {//if different
                            ReadClusterToExternalVectorMC(temp_cluster_id_r, EMMCEdges);
                            cluster_to_bi[temp_cluster_id_r] = REVERSE;
                            ++p_num;
                        }
                    }
                }
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
        MCHotPool<VectorMCEdgesEMTuple_IO> mc_HotPool;       //Hot pool for one-pass algorithms
        mc_HotPool.init(partition_number);
        mc_HotPool.set_capacity(Partition_N);
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
        unordered_set<int> set_remain;//set for storing the id of unshared criteria

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
                while (mc_closed[cri_i][EM_MC_PQueue[cri_i]->top().id]) {//if found in visited set
                    EM_MC_PQueue[cri_i]->pop();
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
                    if (vertex_cri[item_id] == 0) {//evict valid in-memory vertices immediately
                        if (get<0>(mc_HotPool.clusterStatus[partition_id]) == 'I') {
                            get<0>(EMEdgesIndex[item_id]) = true;
                            get<2>(mc_HotPool.clusterStatus[partition_id]) -= temp_degree;
                            if (get<2>(mc_HotPool.clusterStatus[partition_id]) <= 0) {//if empty, clear it immediately
                                assert(get<2>(mc_HotPool.clusterStatus[partition_id]) == 0);
                                for(auto it=mc_HotPool.HotPools[partition_id].begin();it!=mc_HotPool.HotPools[partition_id].end();++it){
                                    get<1>(EMEdgesIndex[it->ID1]) = -1;//hard delete
                                    get<2>(EMEdgesIndex[it->ID1]) = -1;
                                }
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
//        cout << "Query Time: " << tt.GetRuntime() << " s." <<  endl;
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
                for(auto it=mcHotPool.HotPools[evict_p_id].begin();it!=mcHotPool.HotPools[evict_p_id].end();++it){
                    get<1>(EMEdgesIndex[it->ID1]) = -1;//hard delete
                    get<2>(EMEdgesIndex[it->ID1]) = -1;
//                    get<0>(EMEdgesIndex[it->ID1]) = false;
                }
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
//                    else{
//                        get<1>(EMEdgesIndex[it->ID1]) = -1;//hard delete
//                        get<2>(EMEdgesIndex[it->ID1]) = -1;
//                    }
                }
                get<2>(EMEdgesIndex[temp_id]) = mcHotPool.MCEdges_EM.size();
                tt.stop();
                io_time += tt.GetRuntime() * 1000;

                get<0>(mcHotPool.clusterStatus[evict_p_id]) = 'S';
                mcHotPool.HotPools[evict_p_id].clear();//clear
                mcHotPool.HotPools[evict_p_id].shrink_to_fit();//free memory
//                mcHotPool.im_partitions.erase(evict_p_id);
            }
            else {
                ++partition_left;
                if (partition_left > (partition_number/2)) {
                    if (partition_left % 10 == 0) {
                        cout << "The number of partitions which have less than 4KB memory consumption: " << partition_left << endl;
                    }
                }
            }
        }
        
        //// read partition from original disk file to internal vector
        //check status of storage
        char filePath[100];
        int u, u_deg;
        int v, weight;
        MCEdgeT temp_edge;

        strcpy(filePath, DataPath);
        strcat(filePath, dataset.c_str());
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
        MCHotPool<VectorMCEdgesEMTuple_IO> mc_HotPool;       //Hot pool for one-pass algorithms
        mc_HotPool.init(partition_number);
        mc_HotPool.set_capacity(Partition_N);
        vector<vector<Distance>> mc_cost(num_of_cri, vector<Distance>(node_num, INF));   //vector of cost
        //vector<vector<NodeId>> mc_pre(num_of_cri, vector<int>(node_num, -1));      //vector of predecessor id
        vector<vector<bool>> mc_closed(num_of_cri, vector<bool>(node_num));    //flag of vertex closing
        vector<int> cluster_to_cri(partition_number,0);//used to map which criterion read this partition
        unordered_set<int> set_remain;//set for storing the id of unshared criteria
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

        Timer tt;
        tt.start();

        //Initiation of start node
        for (int i = 0; i < num_of_cri; ++i) {
            set_cri.insert(i);
            mc_cost[i][node_start] = 0;//cost of start node
            EM_MC_PQueue[i]->push(VertexCost(node_start,0));
            cri_to_cluster[i] = node_to_cluster[node_start];
            cri_pqueue.update(i,0);
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

                while (mc_closed[cri_i][EM_MC_PQueue[cri_i]->top().id]) {//if found in visited set
                    EM_MC_PQueue[cri_i]->pop();
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
                        if (vertex_cri[item_id] == 0) {//evict valid in-memory vertices immediately
                            if (get<0>(mc_HotPool.clusterStatus[partition_id]) == 'I') {
                                get<0>(EMEdgesIndex[item_id]) = true;
                                get<2>(mc_HotPool.clusterStatus[partition_id]) -= temp_degree;
                                if (get<2>(mc_HotPool.clusterStatus[partition_id]) <= 0) {//if empty, clear it immediately
                                    assert(get<2>(mc_HotPool.clusterStatus[partition_id]) == 0);
                                    for(auto it=mc_HotPool.HotPools[partition_id].begin();it!=mc_HotPool.HotPools[partition_id].end();++it){
                                        get<1>(EMEdgesIndex[it->ID1]) = -1;//hard delete
                                        get<2>(EMEdgesIndex[it->ID1]) = -1;
                                    }
                                    mc_HotPool.HotPools[partition_id].clear();
                                    mc_HotPool.HotPools[partition_id].shrink_to_fit();//free memory
//                                    mc_HotPool.im_partitions.erase(partition_id);
                                    mc_HotPool.LRU_IMCluster.erase(partition_id);
                                    get<0>(mc_HotPool.clusterStatus[partition_id]) = 'D';
                                }
                            }
                        }

                        while (mc_closed[cri_i][EM_MC_PQueue[cri_i]->top().id]) {//if already closed
                            EM_MC_PQueue[cri_i]->pop();
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
            }
        }
        tt.stop();
        query_time += tt.GetRuntime() * 1000;
//        cout << "Query Time: " << tt.GetRuntime() << " s." << endl;
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
        MCHotPool<VectorMCEdgesEMTuple_IO_Bi> mc_HotPool;       //Hot pool for one-pass algorithms
        mc_HotPool.init(partition_number);
        mc_HotPool.set_capacity(Partition_N_Bi);
        unordered_set<int> set_remain;//set for storing the id of unshared criteria
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

        cluster_to_cri[node_to_cluster[node_start]] = 0;
        cluster_to_cri[node_to_cluster[node_end]] = 0;
        //Iteration
        while (!EM_JudgeEmptyBi(EM_MC_PQueue, EM_MC_PQueue_r)) {//for every node in pqueue !mc_finished[0] || !mc_finished[1] || !mc_finished[2]
            set_remain = set_cri;
            /// read by synchronizing the hops 
            temp_cluster_id = cri_to_cluster[cri_pqueue.top_value()];
            GraphMCReadCluster_New_Bi(mc_HotPool, item_id, temp_cluster_id, cluster_to_bi, FORWARD, clusterRead);//if in disk, read corresponding partition
            cluster_to_bi[temp_cluster_id] = FORWARD;
            ++p_num;
            //                    }
            temp_cluster_id_r = cri_to_cluster_r[cri_pqueue_r.top_value()];
            if (temp_cluster_id != temp_cluster_id_r) {//if different
                GraphMCReadCluster_New_Bi(mc_HotPool, item_id_r, temp_cluster_id_r, cluster_to_bi, REVERSE, clusterRead);//if in disk, read corresponding partition
                cluster_to_bi[temp_cluster_id_r] = REVERSE;
                ++p_num;
            }
            while (!set_remain.empty()) {
                cri_i = *set_remain.begin();//pick one criterion to process
                //deal with invalid top elements
                while (mc_closed[cri_i][EM_MC_PQueue[cri_i]->top().id]) {
                    EM_MC_PQueue[cri_i]->pop();
                }
                while (mc_closed_r[cri_i][EM_MC_PQueue_r[cri_i]->top().id]) {
                    EM_MC_PQueue_r[cri_i]->pop();
                }
                //// Termination judging
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
                else {//if not terminate
                    //// Forward searching
                    item_id = EM_MC_PQueue[cri_i]->top().id;// top min item
                    item_dis = EM_MC_PQueue[cri_i]->top().cost;// top min item

                    partition_id = node_to_cluster[item_id];
                    while (get<2>(EMEdgesIndex_Bi[item_id]) != -1 && !mc_closed_r[cri_i][item_id]) {//if the vertex is in memory or stxxl vector.
                        ++cri_hops[cri_i];
                        if (item_id == node_end) {
                            mc_finished[cri_i] = true;
                            mc_min_cost[cri_i] = mc_cost[cri_i][item_id];
                            mc_terminate_id[cri_i] = item_id;
                            set_cri.erase(cri_i);
                            cri_hops[cri_i] = INF;
                            cri_pqueue.update(cri_i, INF);
                            cri_hops_r[cri_i] = INF;
                            cri_pqueue_r.update(cri_i, INF);
                            cout<<node_start<<" "<<node_end<<endl;
                            cout<<"Minimal cost of criterion "<<cri_i<<": "<<mc_min_cost[cri_i]<<endl;
                            flag_mcFinished = true;
                            for (auto it = mc_finished.begin(); it != mc_finished.end(); ++it) {
                                if (!*it) {//if none of them is false, then flag_mcFinished = true
                                    flag_mcFinished = false;
                                    break;
                                }
                            }
                            break;
                        }

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
                            if (vertex_cri_Bi[item_id].first == 0) {//set the flag of evict to true
//                                assert(get<0>(mc_HotPool.clusterStatus[partition_id]) == 'I');
                                if (get<0>(mc_HotPool.clusterStatus[partition_id]) == 'I') {
                                    get<0>(EMEdgesIndex_Bi[item_id]) = true;//soft delete
//                                    get<2>(EMEdgesIndex_Bi[item_id]) = -1; get<3>(EMEdgesIndex_Bi[item_id]) = -1;
                                    get<2>(mc_HotPool.clusterStatus[partition_id]) -= temp_degree;
                                    if (get<2>(mc_HotPool.clusterStatus[partition_id]) <= 0) {//if empty, evict the partition immediately
                                        assert(get<2>(mc_HotPool.clusterStatus[partition_id]) == 0);
                                        for (auto it = mc_HotPool.HotPools[partition_id].begin(); it != mc_HotPool.HotPools[partition_id].end(); ++it) {
                                            get<2>(EMEdgesIndex_Bi[it->ID1]) = -1;//hard delete
                                            get<3>(EMEdgesIndex_Bi[it->ID1]) = -1;
                                        }
                                        mc_HotPool.HotPools[partition_id].clear();
                                        mc_HotPool.HotPools[partition_id].shrink_to_fit();//free memory
//                                        mc_HotPool.im_partitions.erase(partition_id);
                                        mc_HotPool.LRU_IMCluster.erase(partition_id);
                                        get<0>(mc_HotPool.clusterStatus[partition_id]) = 'D';
                                    }
                                }
                            }
                        }

                        while (mc_closed[cri_i][EM_MC_PQueue[cri_i]->top().id]) {//if already closed
                            EM_MC_PQueue[cri_i]->pop();
                        }
                        item_id = EM_MC_PQueue[cri_i]->top().id;// top min item
                        item_dis = EM_MC_PQueue[cri_i]->top().cost;// top min item

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
                        ++cri_hops_r[cri_i];
                        if (item_id_r == node_start) {
                            mc_finished[cri_i] = true;
                            mc_min_cost[cri_i] = mc_cost_r[cri_i][item_id_r];
                            mc_terminate_id[cri_i] = item_id_r;
                            set_cri.erase(cri_i);
                            cri_hops[cri_i] = INF;
                            cri_pqueue.update(cri_i, INF);
                            cri_hops_r[cri_i] = INF;
                            cri_pqueue_r.update(cri_i, INF);
                            cout<<node_start<<" "<<node_end<<endl;
                            cout<<"Minimal cost of criterion "<<cri_i<<": "<<mc_min_cost[cri_i]<<endl;
                            flag_mcFinished = true;
                            for (auto it = mc_finished.begin(); it != mc_finished.end(); ++it) {
                                if (!*it) {//if none of them is false, then flag_mcFinished = true
                                    flag_mcFinished = false;
                                    break;
                                }
                            }
                            break;
                        }
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
                                if (mc_HotPool.MCEdges_EM[i].ID1 != item_id_r) {
                                    cout << item_id_r << " " << mc_HotPool.MCEdges_EM[i].ID1 << endl;
                                }
                                assert(mc_HotPool.MCEdges_EM[i].ID1 == item_id_r);
                                ++temp_degree;
                                temp_id_r = mc_HotPool.MCEdges_EM[i].ID2;
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
                            if (vertex_cri_Bi[item_id_r].second == 0) {//evict valid in-memory vertices immediately
//                                assert(get<0>(mc_HotPool.clusterStatus[partition_id_r]) == 'I');
                                if (get<0>(mc_HotPool.clusterStatus[partition_id_r]) == 'I') {
                                    get<1>(EMEdgesIndex_Bi[item_id_r]) = true;//soft delete
//                                    get<2>(EMEdgesIndex_Bi[item_id_r]) = -1; get<3>(EMEdgesIndex_Bi[item_id_r]) = -1;
                                    get<2>(mc_HotPool.clusterStatus[partition_id_r]) -= temp_degree;
                                    if (get<2>(mc_HotPool.clusterStatus[partition_id_r]) <= 0) {//if empty, clear it immediately
                                        assert(get<2>(mc_HotPool.clusterStatus[partition_id_r]) == 0);
                                        for (auto it = mc_HotPool.HotPools[partition_id_r].begin(); it != mc_HotPool.HotPools[partition_id_r].end(); ++it) {
                                            get<2>(EMEdgesIndex_Bi[it->ID1]) = -1;//hard delete
                                            get<3>(EMEdgesIndex_Bi[it->ID1]) = -1;
                                        }
                                        mc_HotPool.HotPools[partition_id_r].clear();
                                        mc_HotPool.HotPools[partition_id_r].shrink_to_fit();//free memory
//                                        mc_HotPool.im_partitions.erase(partition_id_r);
                                        mc_HotPool.LRU_IMCluster.erase(partition_id_r);
                                        get<0>(mc_HotPool.clusterStatus[partition_id_r]) = 'D';
                                    }
                                }
                            }
                        }

                        while (mc_closed_r[cri_i][EM_MC_PQueue_r[cri_i]->top().id]) {
                            EM_MC_PQueue_r[cri_i]->pop();
                        }
                        item_id_r = EM_MC_PQueue_r[cri_i]->top().id;// top min item
                        item_dis_r = EM_MC_PQueue_r[cri_i]->top().cost;// top min item

                        partition_id_r = node_to_cluster[item_id_r];
                    }
                    //cri_to_cluster_r[cri_i] = node_to_cluster[item_id_r];
                    //cluster_to_cri[cri_to_cluster_r[cri_i]] = cri_i;//node_to_cluster[item_id_r]
                    //cri_pqueue_r.update(cri_i, cri_hops_r[cri_i]);
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
                
                // judge reading
                //if (set_remain.empty() && !flag_mcFinished) {//if all criteria are processed and not all criteria are finished
                //    /// read by synchronizing the hops 
                //    temp_cluster_id = cri_to_cluster[cri_pqueue.top_value()];
                //    GraphMCReadCluster_New_Bi(mc_HotPool, item_id, temp_cluster_id, cluster_to_bi, FORWARD, clusterRead);//if in disk, read corresponding partition
                //    cluster_to_bi[temp_cluster_id] = FORWARD;
                //    ++p_num;
                //    //                    }
                //    temp_cluster_id_r = cri_to_cluster_r[cri_pqueue_r.top_value()];
                //    if (temp_cluster_id != temp_cluster_id_r) {//if different
                //        GraphMCReadCluster_New_Bi(mc_HotPool, item_id_r, temp_cluster_id_r, cluster_to_bi, REVERSE, clusterRead);//if in disk, read corresponding partition
                //        cluster_to_bi[temp_cluster_id_r] = REVERSE;
                //        ++p_num;
                //    }
                //}
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
        mc_HotPool.clear();
    }
    //function of reading partition for one-pass algorithms with io optimization
    void EMMCGraph::GraphMCReadCluster_New_Bi(MCHotPool<VectorMCEdgesEMTuple_IO_Bi> & mcHotPool,NodeId& ID1,int target_p_id,vector<int8_t>& cluster_to_bi,int8_t direction,vector<bool>& clusterRead){
        int evict_p_id = mcHotPool.LRU_IMCluster.put(target_p_id,1);
        //// evict old partition for new one
        if(evict_p_id != -1){//if it is necessary to evict old partition
            assert(get<0>(mcHotPool.clusterStatus[evict_p_id]) == 'I');
            unsigned int evict_p_size = get<2>(mcHotPool.clusterStatus[evict_p_id]);
            uint temp_sz = get<1>(mcHotPool.clusterStatus[evict_p_id])*alpha_bi;//
            assert(temp_sz>0);
            assert(evict_p_size>0);

            if(evict_p_size > temp_sz){//if the remaining size is too large or no node remain, evict all
//                get<2>(mcHotPool.clusterStatus[evict_p_id]) = 0;
                get<0>(mcHotPool.clusterStatus[evict_p_id]) = 'D';
//                cluster_to_bi[evict_p_id] = 0;
                for(auto it=mcHotPool.HotPools[evict_p_id].begin();it!=mcHotPool.HotPools[evict_p_id].end();++it){
                    get<2>(EMEdgesIndex_Bi[it->ID1]) = -1;//hard delete
                    get<3>(EMEdgesIndex_Bi[it->ID1]) = -1;
                }
                mcHotPool.HotPools[evict_p_id].clear();//clear
                mcHotPool.HotPools[evict_p_id].shrink_to_fit();//free memory
//                mcHotPool.im_partitions.erase(evict_p_id);
            }else if(evict_p_size > MCEDGE_PER_BLK){//if the remaining size is moderate, store vertices in temp file; if the remaining size is smaller than block size, remain the vertices in hop pool
                assert(mcHotPool.HotPools[evict_p_id].size()>0);
                /// write to stxxl vector
                int temp_id = -1;
                Timer tt;
                tt.start();
                /// The partition that is stored in stxxl vector has limited information for the other search, therefore, it may need to be re-read from the disk. There are two cases that the partition need to be written into stxxl vector: 1) the first time for one directional search to write, i.e. get<0>(mcHotPool.clusterStatus[evict_p_id]) = 'I'; 2) the second time for the other directional search to write, in this case, all the other useful vertices that are not in stxxl vector will be written into stxxl vector.
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
                    cout << "The second write of partition " <<evict_p_id<<endl;
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
        char filePath[100];
        int u = 0, u_deg;
        int v, weight;
        MCEdgeT temp_edge;

        strcpy(filePath, DataPath);
        strcat(filePath, dataset.c_str());
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
    bool EMMCGraph::EM_JudgeEmptyBi(vector<PriorityQueue*> &em_mc_pqueue, vector<PriorityQueue*> &em_mc_pqueue_r){
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
    }
    //function of checking if the memory setting is feasible
    void EMMCGraph::MemoryCheck(int algo_choice, uint mem, bool& flag_exit) {
        uint basic_mem;

        ClusterInfoLoad(false);
        cout << "Node number: "<<node_num<<"\tEdge number: "<<edge_num<<endl;

        switch (algo_choice){
            case EMDijk:{
                basic_mem = node_num*12/(1024*1024) + WeightPowMax*8 + 32;
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
                basic_mem = node_num*18/(1024*1024) + 32;
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
                basic_mem = node_num*23/(1024*1024) + 32;
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
                basic_mem = node_num*(num_of_cri*5 + 11)/(1024*1024) + 32*num_of_cri;
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
                basic_mem = node_num*(num_of_cri*5 + 11)/(1024*1024) + 32*num_of_cri;
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
                basic_mem = node_num*(num_of_cri*10 + 11)/(1024*1024) + 64*num_of_cri;
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
                basic_mem = node_num*(num_of_cri*5 + 12)/(1024*1024) + 32*num_of_cri;
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
                basic_mem = node_num*(num_of_cri*5 + 12)/(1024*1024) + 32*num_of_cri;
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
                basic_mem = node_num*(num_of_cri*10 + 12)/(1024*1024) + 64*num_of_cri;
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
}
#endif //MCSPS_EMGRAPH_HPP
