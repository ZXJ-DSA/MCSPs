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
        string filename = string(DataPath) + dataset + "/" + aggregateStrategy+"/"+dataset + "_Partitions_" +partMethod+"_"+ to_string(PartitionSize) + "_Info.txt";
//        string filename = string(DataPath) + dataset + "/" + dataset + "_Partitions_Info.txt";
        int ID;

        //Open file
        ifstream inFile(filename, ios::in);
        if (!inFile) {
            cout << "File opening failed. " << filename<< endl;
            exit(1);
        }
        inFile >> node_num >> edge_num;
        mc_criteria.assign(NUM_OF_CRITERIA,string());
        for(int i=0;i<NUM_OF_CRITERIA;++i){
            inFile >> mc_criteria[i];
        }
        inFile >> partition_number;
        cout << "Partition number: " << partition_number <<endl;
        cluster_to_node.assign(partition_number,vector<int>());
        cluster_to_IO.assign(partition_number,0);
        int pSize=0;

        node_to_cluster.assign(node_num,0);
        for(int pid=0;pid<partition_number;++pid){
            inFile >> pSize;
            for(int i=0;i<pSize;++i){
                inFile >> ID;
                cluster_to_node[pid].push_back(ID);
                node_to_cluster[ID] = pid;
            }
            assert(cluster_to_node[pid].size()==pSize);
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
    //Disk-based MCSP Evaluation, new
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

        if(BlockPerPage*Block_SZ > PartitionSize*1024){//page size should be smaller than partition size
            cout<<"The partition size is smaller than the block size! "<<PartitionSize<<" KB"<<endl;
        }
        /*** Read graph data ***/
        if(ifShortcut){
            ShortcutConstruction();
        }else{
            ClusterInfoLoad(true);//load partition information, using node-to-cluster mapping
            ///for one-off read testing
//        string gr_Edge = string(DataPath)  + dataset + "/" + dataset + "_MCEdges.txt";
//        ReadGraphToExternalVector(gr_Edge);
//        MCReadGraph_MCEdges(gr_Edge);
            cout << "Node number: " << node_num << "\t\tEdge number: " << edge_num << endl;
        }

        if(!ifIOOptimization){//no IO optimization
            if(ifShortcut){
                cout<<"Wrong algorithm."<<endl; exit(1);
            }
            if(!ifBidirectional){//unidirectional search
                if(algo==0){
                    /// execution of Dijkstra with external vector
                    cout << "----- Below are results of Dijkstra_IO algorithm -----" << endl;
                    time_DijkstraIO = MC_Multipass(qtype,DijkstraIO);
                    cout << "-----------------------------------------------\n" << endl;
                } else if(algo==1){
                    /// execution of MC_OneHop without IO optimization
                    cout << "----- Below are results of MC_OneHop_NoIO algorithm -----" << endl;
                    time_OneHopNoIO = MC_OnePass(qtype,OneHopNoIO);
                    cout << "-----------------------------------------------\n" << endl;
                    cout << "The total run time of MC_OneHop_NoIO algorithm is: "<< time_OneHopNoIO /1000<<" s."<<endl;
                }else if(algo==2) {//0:Multi-pass, 1:OHP, 2:MHP
                    /// execution of MC_MultiHop without IO optimization
                    cout << "----- Below are results of MC_MultiHop_NoIO algorithm -----" << endl;
                    time_MultiHopsNoIO = MC_OnePass(qtype,MultiHopsNoIO);
                    cout << "-----------------------------------------------\n" << endl;
                    cout << "The total run time of MC_MultiHop_NoIO algorithm is: "<< time_MultiHopsNoIO /1000<<" s."<<endl;
                }else{
                    cout<<"Wrong."<<endl; exit(1);
                }
            }else{//bidirectional search
                if(algo==0){
                    /// execution of BiDijkstra with external vector
                    cout << "----- Below are results of BiDijkstra_IO algorithm -----" << endl;
                    time_BiDijkstraIO = MC_Multipass(qtype,BiDijkstraIO);
                    cout << "-----------------------------------------------\n" << endl;
                }else if(algo==2) {//1:OHP, 2:MHP, 3:SCP
                    /// execution of MC_BiMultiHop without IO optimization
                    cout << "----- Below are results of MC_BiMultiHop_NoIO algorithm -----" << endl;
                    time_BiMultiHopsNoIO = MC_OnePass(qtype,BiMultiHopsNoIO);
                    cout << "-----------------------------------------------\n" << endl;
                    cout << "The total run time of MC_BiMultiHop_NoIO algorithm is: " << time_BiMultiHopsNoIO / 1000 << " s." << endl;
                }else{
                    cout<<"Wrong."<<endl; exit(1);
                }
            }
        }else{//use IO optimization
            if(!ifBidirectional){//unidirectional search
                if(!ifShortcut){
                    if(algo==0){
                        /// execution of EM_Dijk
                        cout << "----- Below are results of EM_Dijk algorithm -----" << endl;
                        time_EMDijk = MC_Multipass(qtype,EMDijk);
                        cout << "-----------------------------------------------\n" << endl;
                    }else if(algo==1){
                        /// execution of MC_OneHop
                        cout << "----- Below are results of MC_OneHop algorithm -----" << endl;
                        time_OneHop = MC_OnePass(qtype,OneHop);
                        cout << "-----------------------------------------------\n" << endl;
                        cout << "The total run time of MC_OneHop algorithm is: "<< time_OneHop /1000<<" s."<<endl;

                    }else if(algo==2) {//1:OHP, 2:MHP, 3:SCP
                        /// execution of MC_MultiHop
                        cout << "----- Below are results of MC_MultiHop algorithm -----" << endl;
                        time_MultiHops = MC_OnePass(qtype,MultiHops);
                        cout << "-----------------------------------------------\n" << endl;
                        cout << "The total run time of MC_MultiHop algorithm is: "<< time_MultiHops /1000<<" s."<<endl;
                    }else{
                        cout<<"Wrong."<<endl; exit(1);
                    }
                }else{
                    if(algo==2) {//1:OHP, 2:MHP, 3:SCP
                        /// execution of SMHP
                        cout << "----- Below are results of SMHP algorithm -----" << endl;
                        time_MultiHops = MC_OnePass(qtype,MultiHops);
                        cout << "-----------------------------------------------\n" << endl;
                        cout << "The total run time of SMHP algorithm is: "<< time_MultiHops /1000<<" s."<<endl;

                    }else{
                        cout<<"Wrong."<<endl; exit(1);
                    }
                }

            }else{//bidirectional search
                if(algo==2) {//1:OHP, 2:MHP, 3:SCP
                    if(ifShortcut){
                        /// execution of BMHPS
                        cout << "----- Below are results of BMHPS algorithm -----" << endl;
                        time_BiMultiHops = MC_OnePass(qtype,BiMultiHops);
                        cout << "-----------------------------------------------\n" << endl;
                        cout << "The total run time of BSMHP algorithm is: " << time_BiMultiHops / 1000 << " s." << endl;
                    }else{
                        /// execution of MC_BiMultiHop
                        cout << "----- Below are results of MC_BiMultiHop algorithm -----" << endl;
                        time_BiMultiHops = MC_OnePass(qtype,BiMultiHops);
                        cout << "-----------------------------------------------\n" << endl;
                        cout << "The total run time of MC_BiMultiHop algorithm is: " << time_BiMultiHops / 1000 << " s." << endl;
                    }

                }else{
                    cout<<"Wrong."<<endl; exit(1);
                }
            }
        }


        node_to_cluster.clear();
    }
    //Disk-based MCSP Evaluation
    /*void EMMCGraph::MC_Evaluate_IO(const string& qtype)
    {
        ///--variables about time record--
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
        //// Read graph data
        ClusterInfoLoad(true);//load partition information, using node-to-cluster mapping
        ///for one-off read testing
//        string gr_Edge = string(DataPath)  + dataset + "/" + dataset + "_MCEdges.txt";
//        ReadGraphToExternalVector(gr_Edge);
//        MCReadGraph_MCEdges(gr_Edge);
        cout << "Node number: " << node_num << "\t\tEdge number: " << edge_num << endl;

        ///// MultiPass algorithm
//        cout << "----- Below are results of EM_Dijk algorithm -----" << endl;
//        time_EMDijk = MC_Multipass(qtype,EMDijk);
//        cout << "-----------------------------------------------\n" << endl;
        /// execution of Dijkstra with external vector
//        cout << "----- Below are results of Dijkstra_IO algorithm -----" << endl;
//        time_DijkstraIO = MC_Multipass(qtype,DijkstraIO);
//        cout << "-----------------------------------------------\n" << endl;
        /// execution of BiDijkstra with external vector
        cout << "----- Below are results of Dijkstra_IO algorithm -----" << endl;
        time_BiDijkstraIO = MC_Multipass(qtype,BiDijkstraIO);
        cout << "-----------------------------------------------\n" << endl;
        ///// OnePass algorithms
        /// execution of MC_OneHop without IO optimization
//        cout << "----- Below are results of MC_OneHop_NoIO algorithm -----" << endl;
//        time_OneHopNoIO = MC_OnePass(qtype,OneHopNoIO);
//        cout << "-----------------------------------------------\n" << endl;
//        /// execution of MC_MultiHop without IO optimization
//        cout << "----- Below are results of MC_MultiHop_NoIO algorithm -----" << endl;
//        time_MultiHopsNoIO = MC_OnePass(qtype,MultiHopsNoIO);
//        cout << "-----------------------------------------------\n" << endl;
//        /// execution of MC_BiMultiHop without IO optimization
//        cout << "----- Below are results of MC_BiMultiHop_NoIO algorithm -----" << endl;
//        time_BiMultiHopsNoIO = MC_OnePass(qtype,BiMultiHopsNoIO);
//        cout << "-----------------------------------------------\n" << endl;
//        /// execution of MC_OneHop
//        cout << "----- Below are results of MC_OneHop algorithm -----" << endl;
//        time_OneHop = MC_OnePass(qtype,OneHop);
//        cout << "-----------------------------------------------\n" << endl;
//        /// execution of MC_MultiHop
//        cout << "----- Below are results of MC_MultiHop algorithm -----" << endl;
//        time_MultiHops = MC_OnePass(qtype,MultiHops);
//        cout << "-----------------------------------------------\n" << endl;
        /// execution of MC_BiMultiHop
        cout << "----- Below are results of MC_BiMultiHop algorithm -----" << endl;
        time_BiMultiHops = MC_OnePass(qtype,BiMultiHops);
        cout << "-----------------------------------------------\n" << endl;
        cout << "The total run time of EM-Dijk algorithm is: "<<time_EMDijk/1000<<" s."<<endl;
        cout << "The total run time of Dijktra_IO algorithm is: "<<time_DijkstraIO/1000<<" s."<<endl;
        cout << "The total run time of BiDijktra_IO algorithm is: "<<time_BiDijkstraIO/1000<<" s."<<endl;
        cout << "The total run time of MC_OneHop_NoIO algorithm is: "<< time_OneHopNoIO /1000<<" s."<<endl;
        cout << "The total run time of MC_MultiHop_NoIO algorithm is: "<< time_MultiHopsNoIO /1000<<" s."<<endl;
        cout << "The total run time of MC_BiMultiHop_NoIO algorithm is: " << time_BiMultiHopsNoIO / 1000 << " s." << endl;
        cout << "The total run time of MC_OneHop algorithm is: "<< time_OneHop /1000<<" s."<<endl;
        cout << "The total run time of MC_MultiHop algorithm is: "<< time_MultiHops /1000<<" s."<<endl;
        cout << "The total run time of MC_BiMultiHop algorithm is: " << time_BiMultiHops / 1000 << " s." << endl;
        node_to_cluster.clear();
    }*/
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
        //std::chrono::high_resolution_clock::time_point t_s, t_e;//variables for time record
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
        /// execution of MC_BiMultiHop without IO optimization
        cout << "----- Below are results of MC_BiMultiHop_NoIO algorithm -----" << endl;
        cout << "The graph data memory size is " << MemGraph_NoIO_Bi << " MB." << endl;
        time_BiMultiHopsNoIO = MC_OnePass(qtype,BiMultiHopsNoIO);
        cout << "-----------------------------------------------\n" << endl;
        /// execution of MC_MultiHop without IO optimization
        cout << "----- Below are results of MC_MultiHop_NoIO algorithm -----" << endl;
        cout << "The graph data memory size is " << MemGraph_NoIO << " MB." << endl;
        time_MultiHopsNoIO = MC_OnePass(qtype,MultiHopsNoIO);
        cout << "-----------------------------------------------\n" << endl;

        /// execution of MC_OneHop
        cout << "----- Below are results of MC_OneHop algorithm -----" << endl;
        cout << "The graph data memory size is " << MemGraph_IO << " MB." << endl;
        time_OneHop = MC_OnePass(qtype,OneHop);
        cout << "-----------------------------------------------\n" << endl;
        /// execution of MC_MultiHop
        cout << "----- Below are results of MC_MultiHop algorithm -----" << endl;
        cout << "The graph data memory size is " << MemGraph_IO << " MB." << endl;
        time_MultiHops = MC_OnePass(qtype,MultiHops);
        cout << "-----------------------------------------------\n" << endl;
        ///// execution of MC_BiMultiHop
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
        double time_OneHopNoIO = 0;
        double time_MultiHopsNoIO = 0;
        double time_BiMultiHopsNoIO = 0;
        double time_OneHop = 0;
        double time_MultiHops = 0;
        double time_BiMultiHops = 0;


        assert(BlockPerPage * Block_SZ <= PartitionSize * 1024);//page size should be smaller than partition size
        /*** Read graph data ***/
        ClusterInfoLoad(true);//load partition information, using node-to-cluster mapping

        cout << "Node number: " << node_num << "\t\tEdge number: " << edge_num << endl;

        /*** MultiPass algorithm ***/
        /// execution of Dijkstra with external vector
        cout << "----- Below are results of Dijkstra_IO algorithm -----" << endl;
        cout << "The graph data memory size is " << MemGraph_DijkIO << " MB." << endl;
        time_DijkstraIO = MC_Multipass(qtype, DijkstraIO);
        cout << "\nThe total run time of Dijktra_IO algorithm is: " << time_DijkstraIO / 1000 << " s." << endl;
        cout << "-----------------------------------------------\n" << endl;
        /// execution of BiDijkstra with external vector
        cout << "----- Below are results of BiDijkstra_IO algorithm -----" << endl;
        time_BiDijkstraIO = MC_Multipass(qtype,BiDijkstraIO);
        cout << "\nThe total run time of BiDijktra_IO algorithm is: "<<time_BiDijkstraIO/1000<<" s."<<endl;
        cout << "-----------------------------------------------\n" << endl;
        /*** OnePass algorithms ***/
        /// execution of MC_OneHop without IO optimization
        cout << "----- Below are results of MC_OneHop_NoIO algorithm -----" << endl;
        cout << "The graph data memory size is " << MemGraph_NoIO << " MB." << endl;
        time_OneHopNoIO = MC_OnePass(qtype,OneHopNoIO);
        cout << "\nThe total run time of MC_OneHop_NoIO algorithm is: "<< time_OneHopNoIO /1000<<" s."<<endl;
        cout << "-----------------------------------------------\n" << endl;
        /// execution of MC_MultiHop without IO optimization
        cout << "----- Below are results of MC_MultiHop_NoIO algorithm -----" << endl;
        cout << "The graph data memory size is " << MemGraph_NoIO << " MB." << endl;
        time_MultiHopsNoIO = MC_OnePass(qtype,MultiHopsNoIO);
        cout << "\nThe total run time of MC_MultiHop_NoIO algorithm is: "<< time_MultiHopsNoIO /1000<<" s."<<endl;
        cout << "-----------------------------------------------\n" << endl;
        /// execution of MC_BiMultiHop without IO optimization
        cout << "----- Below are results of MC_BiMultiHop_NoIO algorithm -----" << endl;
        cout << "The graph data memory size is " << MemGraph_NoIO_Bi << " MB." << endl;
        time_BiMultiHopsNoIO = MC_OnePass(qtype,BiMultiHopsNoIO);
        cout << "\nThe total run time of MC_BiMultiHop_NoIO algorithm is: " << time_BiMultiHopsNoIO / 1000 << " s." << endl;
        cout << "-----------------------------------------------\n" << endl;
        /// execution of MC_OneHop
        cout << "----- Below are results of MC_OneHop algorithm -----" << endl;
        cout << "The graph data memory size is " << MemGraph_IO << " MB." << endl;
        time_OneHop = MC_OnePass(qtype,OneHop);
        cout << "\nThe total run time of MC_OneHop algorithm is: "<< time_OneHop /1000<<" s."<<endl;
        cout << "-----------------------------------------------\n" << endl;
        /// execution of MC_MultiHop
        cout << "----- Below are results of MC_MultiHop algorithm -----" << endl;
        cout << "The graph data memory size is " << MemGraph_IO << " MB." << endl;
        time_MultiHops = MC_OnePass(qtype, MultiHops);
        cout << "\nThe total run time of MC_MultiHop algorithm is: " << time_MultiHops / 1000 << " s." << endl;
        cout << "-----------------------------------------------\n" << endl;
        ///// execution of MC_BiMultiHop
        cout << "The graph data memory size is " << MemGraph_IO_Bi << " MB." << endl;
        cout << "----- Below are results of MC_BiMultiHop algorithm -----" << endl;
        time_BiMultiHops = MC_OnePass(qtype, BiMultiHops);
        cout << "\nThe total run time of MC_BiMultiHop algorithm is: " << time_BiMultiHops / 1000 << " s." << endl;
        cout << "-----------------------------------------------\n" << endl;
        //// EM_Dijk
        cout << "----- Below are results of EM_Dijk algorithm -----" << endl;
        cout << "The graph data memory size is " << MemGraph_EMDijk << " MB." << endl;
        time_EMDijk = MC_Multipass(qtype,EMDijk);
        cout << "\nThe total run time of EM-Dijk algorithm is: " << time_EMDijk / 1000 << " s." << endl;
        cout << "-----------------------------------------------\n" << endl;
        cout << "The total run time of EM-Dijk algorithm is: " << time_EMDijk / 1000 << " s." << endl;
        cout << "The total run time of Dijktra_IO algorithm is: " << time_DijkstraIO / 1000 << " s." << endl;
        cout << "The total run time of BiDijktra_IO algorithm is: "<<time_BiDijkstraIO/1000<<" s."<<endl;
        cout << "The total run time of MC_OneHop_NoIO algorithm is: "<< time_OneHopNoIO /1000<<" s."<<endl;
        cout << "The total run time of MC_MultiHop_NoIO algorithm is: "<< time_MultiHopsNoIO /1000<<" s."<<endl;
        cout << "The total run time of MC_BiMultiHop_NoIO algorithm is: " << time_BiMultiHopsNoIO / 1000 << " s." << endl;
        cout << "The total run time of MC_OneHop algorithm is: "<< time_OneHop /1000<<" s."<<endl;
        cout << "The total run time of MC_MultiHop algorithm is: " << time_MultiHops / 1000 << " s." << endl;
        cout << "The total run time of MC_BiMultiHop algorithm is: " << time_BiMultiHops / 1000 << " s." << endl;
        node_to_cluster.clear();
    }
    //Disk-based MCSP Evaluation: Exp4-Effect of the number of criteria; Exp5-Effect of memory size; Exp6-Effect of partition size; Exp7-Effect of query types
    void EMMCGraph::MC_Evaluate_EXP456(const string& qtype)
    {
        /*--variables about time record--*/
        //double time_EMDijk = 0;
        double time_DijkstraIO = 0;
        double time_BiDijkstraIO = 0;
        double time_OneHop = 0;
        double time_MultiHops = 0;
        double time_BiMultiHops = 0;

        assert(BlockPerPage*Block_SZ <= PartitionSize*1024);//page size should be smaller than partition size
        /*** Read graph data ***/
        ClusterInfoLoad(true);//load partition information, using node-to-cluster mapping
        cout << "Node number: " << node_num << "\t\tEdge number: " << edge_num << endl;

        /*** MultiPass algorithm ***/
//        cout << "----- Below are results of EM_Dijk algorithm -----" << endl;
//        cout << "The graph data memory size is " << MemGraph_EMDijk << " MB." << endl;
//        time_EMDijk = MC_Multipass(qtype,EMDijk);
//        cout << "-----------------------------------------------\n" << endl;
//        /// execution of Dijkstra with external vector
//        cout << "----- Below are results of Dijkstra_IO algorithm -----" << endl;
//        cout << "The graph data memory size is " << MemGraph_DijkIO << " MB." << endl;
//        time_DijkstraIO = MC_Multipass(qtype,DijkstraIO);
//        cout << "-----------------------------------------------\n" << endl;
//        /// execution of BiDijkstra with external vector
//        cout << "----- Below are results of BiDijkstra_IO algorithm -----" << endl;
//        time_BiDijkstraIO = MC_Multipass(qtype,BiDijkstraIO);
//        cout << "-----------------------------------------------\n" << endl;
        /*** OnePass algorithms ***/
//        /// execution of MC_OneHop
//        cout << "----- Below are results of MC_OneHop algorithm -----" << endl;
//        cout << "The graph data memory size is " << MemGraph_IO << " MB." << endl;
//        time_OneHop = MC_OnePass(qtype,OneHop);
//        cout << "\nThe total run time of MC_OneHop algorithm is: "<< time_OneHop /1000<<" s."<<endl;
//        cout << "-----------------------------------------------\n" << endl;
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
//        cout << "The total run time of EM-Dijk algorithm is: "<<time_EMDijk/1000<<" s."<<endl;
        cout << "The total run time of Dijktra_IO algorithm is: "<<time_DijkstraIO/1000<<" s."<<endl;
        cout << "The total run time of BiDijktra_IO algorithm is: "<<time_BiDijkstraIO/1000<<" s."<<endl;
        cout << "The total run time of MC_OneHop algorithm is: "<< time_OneHop /1000<<" s."<<endl;
        cout << "The total run time of MC_MultiHop algorithm is: "<< time_MultiHops /1000<<" s."<<endl;
        cout << "The total run time of MC_BiMultiHop algorithm is: " << time_BiMultiHops / 1000 << " s." << endl;
        node_to_cluster.clear();
    }
    //Function for multipass algorithms
    double EMMCGraph::MC_Multipass(const string & qtype, int algo_choice) {
        string LongDis,MediumDis,ShortDis,RandomDis;
        double total_time = 0;
        ///File path of OD pairs
        LongDis = string(DataPath) + dataset + "/" + aggregateStrategy+"/"+dataset + "_OD_LongDis_" + to_string(PartitionSize) + ".txt";
        MediumDis = string(DataPath) + dataset + "/" + aggregateStrategy+"/"+dataset + "_OD_MediumDis_"+ to_string(PartitionSize) +".txt";
        ShortDis = string(DataPath) + dataset + "/" + aggregateStrategy+"/"+dataset + "_OD_ShortDis_"+ to_string(PartitionSize) +".txt";
        RandomDis = string(DataPath) +  dataset + "/" + aggregateStrategy+"/"+dataset + "_OD_Random_"+ to_string(PartitionSize) +".txt";

        ifshow = true;
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
            QueryType = "R";
            total_time = MC_Multipass_one(RandomDis,algo_choice);//Efficiency evaluation on random distance
        }else if (qtype == "all") {
            QueryType = "S";
            total_time += MC_Multipass_one(ShortDis, algo_choice);//Efficiency evaluation on short distance
            QueryType = "M";
            total_time += MC_Multipass_one(MediumDis, algo_choice);//Efficiency evaluation on medium distance
            QueryType = "L";
            total_time += MC_Multipass_one(LongDis, algo_choice);//Efficiency evaluation on long distance
        }
        PeakMemory();
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
        Distance temp_cost=0;

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
            uint basic_mem = node_num*18/(1024*1024) + PQMemory;
            cout << "*** The basic memory consumption of Dijkstra_IO is: " << basic_mem << " MB. ***"<<  endl;
        }else if(algo_choice == BiDijkstraIO){
            cout << "PageNumber4_BiDijk: " << PageNumber4_BiDijk << endl;
            uint basic_mem = node_num*23/(1024*1024) + PQMemory;
            cout << "*** The basic memory consumption of BiDijkstra_IO is: " << basic_mem << " MB. ***"<<  endl;
        }
        cout << "Query type: " << QueryType << " \tRun times: " << run_times << endl;
        /// MCSPs Query processing
        for (int i = 0; i < run_times; ++i) {//run_times
            ID1 = ODpair[i].first;
            ID2 = ODpair[i].second;
//            cout << "Query " << i << ": " << ID1-1 << " " << ID2-1 << endl;
            //Initiation for each round
            CommonInitiation_IO();
            if(algo_choice == EMDijk){
                ///Shortest path query one criterion by one criterion
                vector<HotPool> hotPools;
                //statistics of IO
                Stats stats_b, stats_e;
                stats_b.start();
                Timer tt3;
                tt3.start();
                for(sc_i=0;sc_i<num_of_cri;++sc_i) {
//                    temp_cost = EM_Dijk(ID1 - 1, ID2 - 1);
                    temp_cost = EM_Dijk(ID1, ID2, hotPools);
                    if(temp_cost<INF100)
                        ave_cost[sc_i] += temp_cost;
                    //reset done bit to 0
                    for (int j = 0; j < node_num; ++j) {
                        node_to_category[sc_i][j].set(0, false);
                    }
                }
                //io record
                tt3.stop();
//                cout << "Query time: " << tt3.GetRuntime() << " s." << endl;
                stats_e = stats_b.get_stats();
                if (stats_e.m_noOfWrites > 0 || stats_e.m_noOfReads > 0) {
                    cout << "There are " << stats_e.m_noOfWrites + stats_e.m_noOfReads << " io for stxxl." << endl;
                    //cout<<stats_e<<endl;
                }
                hotPools.clear();
            }else if(algo_choice == DijkstraIO){
                EMEdgesIndex.assign(node_num,make_tuple(false,-1,-1));
                VectorMCEdgesEMTuple_DijkIO EMMCEdges;
//                HotPool4<VectorMCEdgesEMTuple4_DijkIO> myHotPool(partition_number,PartitionNumber_DijkIO_new);//hot pool
//                VectorMCEdgesEMTuple4_DijkIO_new EMMCEdges;
                //statstics of IO
                Stats stats_b, stats_e;
                stats_b.start();
                ///Shortest path query one criterion by one criterion
                for(sc_i=0;sc_i<num_of_cri;++sc_i) {//num_of_cri
//                    EMEdgesIndex.assign(node_num,make_tuple(false,-1,-1));
//                    VectorMCEdgesEMTuple_DijkIO EMMCEdges;
//                    HotPool4<VectorMCEdgesEMTuple4_DijkIO> myHotPool(partition_number,PartitionNumber_DijkIO_new);
                    temp_cost = Dijkstra_IO(ID1, ID2, EMMCEdges);
//                    temp_cost = Dijkstra_IO_new(ID1, ID2,myHotPool);
                    if(temp_cost<INF100)
                        ave_cost[sc_i] += temp_cost;
//                    EMEdgesIndex.clear();
//                    EMMCEdges.clear();
//                    myHotPool.clear();
                }
                //io record
                stats_e = stats_b.get_stats();
                if (stats_e.m_noOfWrites > 0 || stats_e.m_noOfReads > 0) {
                    cout << "There are " << stats_e.m_noOfWrites + stats_e.m_noOfReads << " io for stxxl." << endl;
                    //cout<<stats_e<<endl;
                }
                EMEdgesIndex.clear();
                EMMCEdges.clear();
//                myHotPool.clear();
            }else if(algo_choice == BiDijkstraIO){
                EMEdgesIndex.assign(node_num,make_tuple(false,-1,-1));
//                VectorMCEdgesEMTuple_BiDijkIO EMMCEdges;
                HotPool4<VectorMCEdgesEMTuple4_BiDijkIO> myHotPool(partition_number,PartitionNumber_BiDijkIO_new);
                //statstics of IO
                Stats stats_b, stats_e;
                stats_b.start();
                ///Shortest path query one criterion by one criterion
                for(sc_i=0;sc_i<num_of_cri;++sc_i) {//0
//                    EMEdgesIndex.assign(node_num,make_tuple(false,-1,-1));
//                    VectorMCEdgesEMTuple_BiDijkIO EMMCEdges;
//                    HotPool4<VectorMCEdgesEMTuple4_BiDijkIO> myHotPool(partition_number,PartitionNumber_BiDijkIO_new);
//                    temp_cost = BiDijkstra_IO(ID1, ID2,EMMCEdges);
                    temp_cost = BiDijkstra_IO_new(ID1, ID2, myHotPool);
                    if(temp_cost<INF100)
                        ave_cost[sc_i] += temp_cost;
//                    EMEdgesIndex.clear();
//                    EMMCEdges.clear();
//                    myHotPool.clear();
                }
                //io record
                stats_e = stats_b.get_stats();
                if (stats_e.m_noOfWrites > 0 || stats_e.m_noOfReads > 0) {
                    cout << "There are " << stats_e.m_noOfWrites + stats_e.m_noOfReads << " io for stxxl." << endl;
                    //cout<<stats_e<<endl;
                }
                EMEdgesIndex.clear();
//                EMMCEdges.clear();
                myHotPool.clear();
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
        string file_r_graph = string(DataPath) + dataset+"/"+aggregateStrategy+"/"+dataset+"_" + to_string(PartitionSize)+".MCEdges";
        string file_r_graph1 = string(DataPath) + this->dataset + "/" + aggregateStrategy+"/"+this->dataset + "_"+ to_string(PartitionSize)+".EMDijk";
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
            string line;
            getline(inFile,line);
            vector<string> re1;
            boost::split(re1,line,boost::is_any_of(" \t"));
            if(re1.size()==2){
                node_num=stoi(re1[0]), edge_num=stoi(re1[1]);
            }else{
                cout<<"wrong syntax. "<<re1.size()<<" "<<line<<endl; exit(1);
            }
            getline(inFile,line);
            re1.clear();
            boost::split(re1,line,boost::is_any_of(" \t"));
            if(re1.size()==6){
                num_cri=stoi(re1[0]);
                cout<<"Number of criteria: "<<num_cri<<endl;
            }else{
                cout<<"wrong syntax. "<<re1.size()<<" "<<line<<endl; exit(1);
            }

            bool flag_double = true;
            if(edge_num>1000000000){
                flag_double=false;
                cout<<"Read edge once!"<<endl;
            }

            max_w.assign(num_cri,0);
            // get node-to-category mapping
            node_to_category.assign(num_cri,vector<bitset1<WeightPowMax+1,EdgeWeight>>(node_num));//bitset<WeightPowMax+1>
            while(getline(inFile,line)){
                if(line.empty())
                    continue;
                re1.clear();
                boost::split(re1,line,boost::is_any_of(" \t"));
                if(re1.size()==7){
                    ID1=stoi(re1[0]), ID2=stoi(re1[1]);
                }
                else{
                    cout<<"Wrong line. "<< line<<endl;
                    exit(1);
                }



                if(!flag_double){
                    for(int k=0;k<num_cri;++k){
                        weight=stoi(re1[k+2]);
                        max_w[k] = max(max_w[k],weight);
                        node_to_category[k][ID1].set(Weight_to_category(weight)+1, true);//record edge category Weight_to_category(weight)+1
                        node_to_category[k][ID2].set(Weight_to_category(weight)+1, true);//record edge category Weight_to_category(weight)+1
                    }
                }else{
                    for(int k=0;k<num_cri;++k){
                        weight=stoi(re1[k+2]);
                        max_w[k] = max(max_w[k],weight);
                        node_to_category[k][ID1].set(Weight_to_category(weight)+1, true);//record edge category Weight_to_category(weight)+1
                    }
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
        }
        else{//if it can be opened
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


        HotPool_EdgeSZ = 0.9*MemGraph_EMDijk*1024*1024/(Edge_SZ*temp_weightPowMax);//update the hot pool edge size
        uint basic_mem = node_num*28/(1024*1024) + WeightPowMax*QMemory/(1024*1024) + PQMemory;
//        node_num*12/(1024*1024) + WeightPowMax*4 + PQMemory;
        cout << "*** The actual basic memory consumption of EM_Dijk is: " << basic_mem << " MB. ***"<<  endl;
        cout << "*** The avoidable memory consumption caused by node_to_category is: " << node_num*16/(1024*1024) << " MB. ***" << endl;

    }
    //function of initialising multi-pass algorithms
    void EMMCGraph::EM_Dijk_Initiation(vector<HotPool>& hotPools){
        double nsize=0;
        hotPools.assign(num_pool[sc_i],HotPool());//reserve space for hot pools
        for(int i=0;i<hotPools.size();++i){
            hotPools[i].EM_Map.assign(partition_number,make_pair(0,0));
            nsize = 1+0.7*hotPools.size()*HotPool_EdgeSZ/(pow(2,i)*(2-pow(0.5,hotPools.size()-1)));
            hotPools[i].set_capacity(nsize);
//            if(dataset == "PopeElection2013" || dataset == "BostonBomb2013"){
////                nsize = 1+0.4*HotPool_EdgeSZ*hotPools.size()/pow(hotPools.size(),i);
//                nsize = 1+hotPools.size()*HotPool_EdgeSZ/(pow(2,i)*(2-pow(0.5,hotPools.size()-1)));
//                hotPools[i].set_capacity(nsize);
//            }
//            else{
//                if(sc_i<2){//for correlated criterion
//                    nsize = HotPool_EdgeSZ*4.0/(abs(i-6)+1);
//                    hotPools[i].set_capacity(nsize);
//                }else{//for random criterion
//                    hotPools[i].set_capacity(HotPool_EdgeSZ);
//                }
//            }
        }
    }
    //Function of Semi-external Memory Dijkstra: new version powered by GraphPool
    Distance EMMCGraph::EM_Dijk(int node_start, int node_end, vector<HotPool>& hotPools){//immediate pop
        if(node_start==node_end){
            cout<<"Same source and target vertex id!"<<endl;
            return 0;
        }
        //statistics of IO
        Stats stats_b, stats_e;
        stats_b.start();

        //Variables
        PriorityQueue pqueue(PQ_Pool);
        NodeId item_id,temp_id;
        Distance item_dis,temp_dis;
        VertexCost item_;      //recording top element of priority queue
        vector<EdgePairW> item_adj; //retrieve result of adjacency list
        bool flag_empty;        //flag of whether the priority queue will be empty if current top element is popped
        Distance relax_i_dis;   //the distance of the first element of relax i
        Distance nextTopDis;    //the distance of the next top element
        Distance min_cost = INF;      //the shortest path distance
        vector<Distance> cost(node_num,INF);  //In-memory vertex cost vector
        //vector<NodeId> pre(node_num,-1);     //In-memory predecessor id
        VertexCost relax_top;  //the front element of relax_i
        vector<em_queue> relaxPools(num_pool[sc_i]);//em_queue()
        partitions_read.assign(partition_number, false);
        int cluster_id;
        vector<bool> clusterIM(partition_number,false); //flag of whether the partition is in memory

        //Hot pools
//        vector<HotPool> hotPools;
        EM_Dijk_Initiation(hotPools);

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
//                cout<<"Minimal cost of criterion " << sc_i << ": "<<min_cost<<endl;
                break;
            }
//            set_visited.insert(item_.id);
            ++visited_number;//update #visited_number
            //flag initiation for each round
            flag_empty = false;

            //eradicate the closed items in priority queue for next round
            while(!pqueue.empty()){
                if(node_to_category[sc_i][pqueue.top().id].test(0)){
                    pqueue.pop();
                }
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
//                                if(ifOptimal){//calculate the minimal partition reading
//                                    set_readCluster.insert(cluster_id);
//                                }
                                ++p_num;
                                item_adj.clear();
                                EM_ReadCluster(category_id,item_id,item_adj,hotPools);
                                clusterIM[cluster_id] = true;
                                assert(!item_adj.empty());
                                //relaxation
                                for(auto it=item_adj.begin();it!=item_adj.end();++it){
                                    temp_id = it->ID2;
                                    assert(temp_id<node_num);
                                    if (it->w == INF100)
                                        continue;
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
                                    //relaxation
                                    auto iter = hotPools[category_id].Edges_IM.equal_range(item_id);
                                    for(auto it=iter.first;it!=iter.second;++it){
                                        temp_id = it->second.ID2;
                                        assert(temp_id<node_num);
                                        if (it->second.w == INF100)
                                            continue;
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
                                }else {//check the external vector if(hotPools[category_id].EM_Map[cluster_id].second > 0)
                                    bool flag_finish = false;
                                    // linear search
                                    for (int i = hotPools[category_id].EM_Map[cluster_id].first; i < hotPools[category_id].EM_Map[cluster_id].second; ++i) {
                                        if (hotPools[category_id].Edges_EM[i].ID1 == item_id) {//if it is the target edge
                                            //relaxation
                                            for (auto it = item_adj.begin(); it != item_adj.end(); ++it) {
                                                temp_id = hotPools[category_id].Edges_EM[i].ID2;
                                                assert(temp_id < node_num);
                                                if (hotPools[category_id].Edges_EM[i].w == INF100)
                                                    continue;
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
                                }
//                                else{
//                                    cout << "Wrong!!!" << endl;
//                                    assert(!clusterIM[cluster_id]);
//                                }
                            }
                            relaxPools[category_id].pop();
                        }
                    }
                }
            }
        }
        tt1.stop();
        query_time += tt1.GetRuntime() * 1000;//
        if (item_.id == node_end) {
            cout << "Query time: " << tt1.GetRuntime() << "s." << endl;
        }
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
//        unordered_map<int,vector<EdgePairW>> tempEdges;//from category to Edges
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
        char filePath[300];
        strcpy(filePath, DataPath);
        strcat(filePath, dataset.c_str());
        strcat(filePath, "/");
        strcat(filePath, aggregateStrategy.c_str());
        strcat(filePath, "/Partitions/");
        strcat(filePath, dataset.c_str());
        strcat(filePath, "_Partition_");
        strcat(filePath, partMethod.c_str());
        strcat(filePath, "_");
        strcat(filePath, to_string(PartitionSize).c_str());
        strcat(filePath, "_");
        strcat(filePath, to_string(target_p_id).c_str());
        strcat(filePath, ".bin");
        FILE* file = fopen(filePath, "rb");
        if(!file){
            cout<<"Failed to open file "<<filePath<<endl;
            exit(1);
        }

        ReadBuffer read_io(file);

        int temp;
        read_io.read(&temp);//cluster id
        read_io.read(&temp);//cluster size

        int u, u_deg;
        int v, weight;
        //read partition
        while(!read_io.isend) {//if not read out
            read_io.read(&u);//get node id
            temp_id = u;
            assert(u<node_num);
            read_io.read(&u_deg);//get node degree
//            tempEdges.clear();
            //vertex store in indexNode, edge store in e
            for (int i = 0; i < u_deg; ++i) {
                read_io.read(&v);//read end node
                assert(v<node_num);
                for (int j = 0; j < NUM_OF_CRITERIA; ++j) {
                    read_io.read(&weight);
                    if (j == sc_i){
//                        tempEdges[Weight_to_category(weight)].emplace_back(EdgePairW(u-1, v-1, weight));
                        /// sort the edge and distribute them over hotpools according to their categories.
                        if(weight == INF100)
                            continue;
                        temp_category = Weight_to_category(weight);
                        temp_edge.ID1 = temp_id; temp_edge.ID2 = v; temp_edge.w = weight;
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
    //function of LRU based Dijkstra's algorithm (LRU+map+stxxl)
    Distance EMMCGraph::Dijkstra_IO_new(int node_start, int node_end, HotPool4<VectorMCEdgesEMTuple4_DijkIO>& myHotPool){
        if(node_start==node_end){
            cout<<"Same source and target vertex id!"<<endl;
            return 0;
        }
        if(ifshow) {
            cout << "LRU+unordered_multimap+stxxl." << endl;
            ifshow= false;
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
//        vector<NodeId> pre(node_num, -1);       //vector of predecessor id
        Distance min_cost = INF;
        int index_begin, index_end;
        int partition_id;


//        string file_w = string(DataPath) + dataset + "/" + dataset + "_dijk_www2.txt";
//        ofstream outFile(file_w, ios::out);
//        if (!outFile) {
//            cout << "File opening failed." << endl;
//            assert(outFile);
//        }

        Timer tt1;
        tt1.start();
        //Initiation of start node
        cost[node_start] = 0;//cost of start node
        pqueue.push(VertexCost(node_start,0));

        //Iteration
        while (!pqueue.empty()) {//for every node in pqueue
            item_ = pqueue.top();
//            if(item_.id == 18076106){
//                cout<<"In Progress! id: " <<item_.id<<" (dis: "<<item_.cost<<", ifClosed: "<< closed[item_.id] << ") is popped from pqueue."<<endl;
//            }
//            outFile << item_.id <<" " << item_.cost << " " << closed[item_.id] << endl;
            pqueue.pop();
            if(closed[item_.id]){//if already closed
                continue;
            }

//            outFile << "In! " << item_.id <<" " << item_.cost << endl;
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
            if(myHotPool.clusterStatus[partition_id] == 'D'){//if in disk
                GraphReadCluster4(myHotPool,partition_id);
                ++p_num;
            }
            if(myHotPool.clusterStatus[partition_id] == 'I'){
                /// read from EMMCEdges
//                auto adjlist = myHotPool.HotPools[partition_id].equal_range(item_id);
//                for(auto it = adjlist.first; it != adjlist.second; ++it){
//                    temp_id = it->second.ID2;
                for(auto it = myHotPool.HotPools[partition_id][item_id].begin(); it != myHotPool.HotPools[partition_id][item_id].end(); ++it){
                    temp_id = it->ID2;
                    if (it->getW(sc_i)==INF100)
                        continue;
                    if (closed[temp_id])//if closed
                        continue;
                    temp_dis = item_dis + it->getW(sc_i);
                    if (cost[temp_id] > temp_dis) {//slack operation
                        cost[temp_id] = temp_dis;
//                        pre[temp_id] = item_id;
//                        if(temp_id == 18076106){
//                            cout<<"I! id: " <<temp_id<<" (dis: "<<temp_dis<<") is pushed into pqueue. item_id: "<<item_id<<endl;
//                        }
                        pqueue.push(VertexCost(temp_id, temp_dis));
                    }
                }
                closed[item_id] = true;
            }else if(myHotPool.clusterStatus[partition_id] == 'S'){
                for (int i = get<1>(EMEdgesIndex[item_id]); i < get<2>(EMEdgesIndex[item_id]); ++i) {
                    assert(myHotPool.MCEdges_EM[i].ID1 == item_id);
                    temp_id = myHotPool.MCEdges_EM[i].ID2;
                    if (myHotPool.MCEdges_EM[i].getW(sc_i)==INF100)
                        continue;
                    if (closed[temp_id])//if closed
                        continue;
                    temp_dis = item_dis + myHotPool.MCEdges_EM[i].getW(sc_i);
                    if (cost[temp_id] > temp_dis) {//slack operation
                        cost[temp_id] = temp_dis;
//                        pre[temp_id] = item_id;
//                        if(temp_id == 18076106){
//                            cout<<"S! id: " <<temp_id<<" (dis: "<<temp_dis<<") is pushed into pqueue. item_id: "<<item_id<<endl;
//                        }
                        pqueue.push(VertexCost(temp_id, temp_dis));
                    }
                }
                closed[item_id] = true;
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
        while(!pqueue.empty()){
//            outFile << pqueue.top().id <<" " << pqueue.top().cost << " " << closed[pqueue.top().id] << endl;
//            outFile << "After! " << pqueue.top().id <<" " << pqueue.top().cost << endl;
//            if(pqueue.top().id == 18076106){
//                cout<<"After Progress! id: " <<pqueue.top().id<<" (dis: "<<pqueue.top().cost<<", ifClosed: "<< closed[pqueue.top().id] << ") is popped from pqueue."<<endl;
//            }
            pqueue.pop();
        }
//        outFile.close();
        return min_cost;
    }
    //function of stxxl vector based Dijkstra's algorithm
    Distance EMMCGraph::Dijkstra_IO(int node_start, int node_end, VectorMCEdgesEMTuple_DijkIO & EMMCEdges){
        if(node_start==node_end){
            cout<<"Same source and target vertex id!"<<endl;
            return 0;
        }
        if(ifshow){
            cout<<"stxxl+index."<<endl;
            ifshow = false;
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
//                cout<<"Minimal cost of criterion " << sc_i << ": "<<min_cost<<endl;
                break;
            }
            //relaxation
//            set_visited.insert(item_id);//update the visited set for this criteria
            ++visited_number;//update #visited_number
            if(get<1>(EMEdgesIndex[item_id]) == -1){ // if the vertex is not in stxxl vector
                ReadClusterToExternalVectorMC(node_to_cluster[item_id],EMMCEdges);
                ++p_num;
                assert(get<1>(EMEdgesIndex[item_id]) != -1);
            }
            /// read from EMMCEdges
            for(int i = get<1>(EMEdgesIndex[item_id]); i < get<2>(EMEdgesIndex[item_id]); ++i){
                assert(EMMCEdges[i].ID1 == item_id);
                temp_id = EMMCEdges[i].ID2;
                if (EMMCEdges[i].getW(sc_i)==INF100)
                    continue;
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
        char filePath[300];
        int u, u_deg;
        int v, weight;
        bool flag_sizeUpdate = false;
        int temp_id = -1;
        MCEdgeT temp_edge;//

        strcpy(filePath, DataPath);
        strcat(filePath, dataset.c_str());
        strcat(filePath, "/");
        strcat(filePath, aggregateStrategy.c_str());
        strcat(filePath, "/Partitions/");
        strcat(filePath, dataset.c_str());
        strcat(filePath, "_Partition_");
        strcat(filePath, partMethod.c_str());
        strcat(filePath, "_");
        strcat(filePath, to_string(PartitionSize).c_str());
        strcat(filePath, "_");
        strcat(filePath, to_string(target_p_id).c_str());
        strcat(filePath, ".bin");
        Timer tt;
        tt.start();
        FILE* file = fopen(filePath, "rb");
        if(!file){
            cout<<"Failed to open file "<<filePath<<endl;
            exit(1);
        }
        assert(file != NULL);
        ReadBuffer read_io(file);

        int temp;
        read_io.read(&temp);//cluster id
        read_io.read(&temp);//cluster size

        //read partition
        while(!read_io.isend) {//if not read out
            read_io.read(&u);//get node id
            //cout << u - 1 << endl;
//            if (u >2)
//                cout << u-1 << endl;
            if(u != temp_id){
                get<1>(EMEdgesIndex[u]) = EMMCEdges.size();
                if(temp_id != -1){
                    get<2>(EMEdgesIndex[temp_id]) = EMMCEdges.size();
                }
                temp_id = u;
            }
            temp_edge.ID1 = u;
            read_io.read(&u_deg);//get node degree
            //vertex store in indexNode, edge store in e
            for (int i = 0; i < u_deg; ++i) {
                read_io.read(&v);//read end node
                temp_edge.ID2 = v;
                for (int j = 0; j < NUM_OF_CRITERIA; ++j) {
                    read_io.read(&weight);
                    if(j < num_of_cri){
                        temp_edge.putW(j,weight);
                    }
                }
                EMMCEdges.push_back(temp_edge);
            }
        }
        get<2>(EMEdgesIndex[u]) = EMMCEdges.size();
        io_num += read_io.io_number;//update #IO
        fclose(file);
        tt.stop();
        io_time += tt.GetRuntime() * 1000;
        if(cluster_to_IO[target_p_id] == 0){
            cluster_to_IO[target_p_id] = read_io.io_number;
        }
    }
    //function of reading graph data to HotPool2
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
        char filePath[300];
        int u, u_deg;
        int v, weight;
//        MCEdgeT temp_edge;

        strcpy(filePath, DataPath);
        strcat(filePath, dataset.c_str());
        strcat(filePath, "/");
        strcat(filePath, aggregateStrategy.c_str());
        strcat(filePath, "/Partitions/");
        strcat(filePath, dataset.c_str());
        strcat(filePath, "_Partition_");
        strcat(filePath, partMethod.c_str());
        strcat(filePath, "_");
        strcat(filePath, to_string(PartitionSize).c_str());
        strcat(filePath, "_");
        strcat(filePath, to_string(target_p_id).c_str());
        strcat(filePath, ".bin");
        Timer tt;
        tt.start();
        FILE* file = fopen(filePath, "rb");
        if(!file){
            cout<<"Failed to open file "<<filePath<<endl;
            exit(1);
        }
        assert(file != NULL);
        ReadBuffer read_io(file);
        bool flag_sizeUpdate = false;
        int temp_id = -1;

        int temp;
        read_io.read(&temp);//cluster id
        read_io.read(&temp);//cluster size

        //read partition
        while(!read_io.isend) {//if not read out
            read_io.read(&u);//get node id
//            temp_edge.ID1 = u-1;
            if(u != temp_id){
                get<1>(EMEdgesIndex[u]) = mcHotPool.HotPools[target_p_id].size();
                if(temp_id != -1){
                    get<2>(EMEdgesIndex[temp_id]) = mcHotPool.HotPools[target_p_id].size();
                }
                temp_id = u;
            }

            read_io.read(&u_deg);//get node degree
            //vertex store in indexNode, edge store in e
            for (int i = 0; i < u_deg; ++i) {
                read_io.read(&v);//read end node
//                temp_edge.ID2 = v-1;
                for (int j = 0; j < NUM_OF_CRITERIA; ++j) {
                    read_io.read(&weight);
                    if(j == sc_i){
                        mcHotPool.HotPools[target_p_id].emplace_back(EdgePairW(u,v,weight));
                    }
                }
            }
        }

        get<2>(EMEdgesIndex[u]) = mcHotPool.HotPools[target_p_id].size();

        io_num += read_io.io_number;//update #IO
        fclose(file);
        tt.stop();
        io_time += tt.GetRuntime() * 1000;
        cluster_to_IO[target_p_id] = read_io.io_number;

        mcHotPool.FlagIM[target_p_id] = true;
    }
    //function of reading graph data to HotPool3
    void EMMCGraph::GraphReadCluster3(HotPool3 & mcHotPool,int target_p_id){
        int evict_p_id = mcHotPool.LRU_IMCluster.put(target_p_id,1);
        //// evict old partition for new one
        if(evict_p_id != -1) {//if it is necessary to evict old partition
            mcHotPool.FlagIM[evict_p_id] = false;
            mcHotPool.HotPools[evict_p_id].clear();
        }
        //// read partition from original disk file to internal vector
        //check status of storage
        char filePath[300];
        int u, u_deg;
        int v, weight;

        strcpy(filePath, DataPath);
        strcat(filePath, dataset.c_str());
        strcat(filePath, "/");
        strcat(filePath, aggregateStrategy.c_str());
        strcat(filePath, "/Partitions/");
        strcat(filePath, dataset.c_str());
        strcat(filePath, "_Partition_");
        strcat(filePath, partMethod.c_str());
        strcat(filePath, "_");
        strcat(filePath, to_string(PartitionSize).c_str());
        strcat(filePath, "_");
        strcat(filePath, to_string(target_p_id).c_str());
        strcat(filePath, ".bin");
        Timer tt;
        tt.start();
        FILE* file = fopen(filePath, "rb");
        if(!file){
            cout<<"Failed to open file "<<filePath<<endl;
            exit(1);
        }
        assert(file != NULL);
        ReadBuffer read_io(file);
        bool flag_sizeUpdate = false;
        int temp_id = 0;

        int temp;
        read_io.read(&temp);//cluster id
        read_io.read(&temp);//cluster size

        //read partition
        while(!read_io.isend) {//if not read out
            read_io.read(&u);//get node id

            read_io.read(&u_deg);//get node degree
            //vertex store in indexNode, edge store in e
            for (int i = 0; i < u_deg; ++i) {
                read_io.read(&v);//read end node

                for (int j = 0; j < NUM_OF_CRITERIA; ++j) {
                    read_io.read(&weight);
                    if(j == sc_i){
                        mcHotPool.HotPools[target_p_id].emplace(u,Edge(v,weight));
                    }
                }
            }
        }

        io_num += read_io.io_number;//update #IO
        fclose(file);
        tt.stop();
        io_time += tt.GetRuntime() * 1000;
        cluster_to_IO[target_p_id] = read_io.io_number;

        mcHotPool.FlagIM[target_p_id] = true;
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
            for(auto it=cluster_to_node[evict_p_id].begin();it!=cluster_to_node[evict_p_id].end();++it){//store all elements to stxxl vector
                if(*it != temp_id){
                    get<1>(EMEdgesIndex[*it]) = mcHotPool.MCEdges_EM.size();
                    if(temp_id != -1){
                        get<2>(EMEdgesIndex[temp_id]) = mcHotPool.MCEdges_EM.size();
                    }
                    temp_id = *it;
                }
                assert(mcHotPool.HotPools[evict_p_id].find(*it)!=mcHotPool.HotPools[evict_p_id].end());
//                auto range=mcHotPool.HotPools[evict_p_id].equal_range(i);
//                for(auto it=range.first;it!=range.second;++it){
//                    mcHotPool.MCEdges_EM.push_back(it->second);
//                }
                for(auto it2=mcHotPool.HotPools[evict_p_id][*it].begin();it2!=mcHotPool.HotPools[evict_p_id][*it].end();++it2){
                    mcHotPool.MCEdges_EM.push_back(*it2);
                }
            }
            get<2>(EMEdgesIndex[temp_id]) = mcHotPool.MCEdges_EM.size();

            mcHotPool.clusterStatus[evict_p_id] = 'S';
            mcHotPool.HotPools[evict_p_id].clear();
        }
        //// read partition from original disk file to internal vector
        //check status of storage
        char filePath[300];
        int u, u_deg;
        int v, weight;

        strcpy(filePath, DataPath);
        strcat(filePath, dataset.c_str());
        strcat(filePath, "/");
        strcat(filePath, aggregateStrategy.c_str());
        strcat(filePath, "/Partitions/");
        strcat(filePath, dataset.c_str());
        strcat(filePath, "_Partition_");
        strcat(filePath, partMethod.c_str());
        strcat(filePath, "_");
        strcat(filePath, to_string(PartitionSize).c_str());
        strcat(filePath, "_");
        strcat(filePath, to_string(target_p_id).c_str());
        strcat(filePath, ".bin");
        Timer tt;
        tt.start();
        FILE* file = fopen(filePath, "rb");
        if(!file){
            cout<<"Failed to open file "<<filePath<<endl;
            exit(1);
        }
        assert(file != NULL);
        ReadBuffer read_io(file);
        bool flag_sizeUpdate = false;
        int temp_id = 0;
        MCEdgeT temp_edge;
        vector<MCEdgeT> temp_vedge;

        int temp;
        read_io.read(&temp);//cluster id
        read_io.read(&temp);//cluster size

        //read partition
        while(!read_io.isend) {//if not read out
            read_io.read(&u);//get node id
            temp_edge.ID1 = u;
            read_io.read(&u_deg);//get node degree
            //vertex store in indexNode, edge store in e
            for (int i = 0; i < u_deg; ++i) {
                read_io.read(&v);//read end node

                temp_edge.ID2 = v;
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
            mcHotPool.HotPools[target_p_id].insert({u,temp_vedge});
            temp_vedge.clear();
        }

        io_num += read_io.io_number;//update #IO
        fclose(file);
        tt.stop();
        io_time += tt.GetRuntime() * 1000;
        cluster_to_IO[target_p_id] = read_io.io_number;

        mcHotPool.clusterStatus[target_p_id] = 'I';
    }
    //function of stxxl vector based BiDijkstra's algorithm
    Distance EMMCGraph::BiDijkstra_IO(int node_start, int node_end, VectorMCEdgesEMTuple_BiDijkIO & EMMCEdges){
        if(node_start==node_end){
            cout<<"Same source and target vertex id!"<<endl;
            return 0;
        }
        if(ifshow){
            cout<<"stxxl+index."<<endl;
            ifshow = false;
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
        while (!pqueue.empty() && !pqueue_r.empty()) {//for every node in pqueue

            //termination judgement
            if (pqueue.top().cost + pqueue_r.top().cost >= min_cost) {//condition of termination
//                cout << node_start << " " << node_end << " " << sc_i << " "<<min_cost << endl;
//                cout<<"Minimal cost of criterion " << sc_i << ": "<<min_cost<<endl;
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
                //cout << item_id << endl;
                assert(get<1>(EMEdgesIndex[item_id]) != -1);
            }
            // read from EMMCEdges
            for(int i = get<1>(EMEdgesIndex[item_id]); i < get<2>(EMEdgesIndex[item_id]); ++i){
                assert(EMMCEdges[i].ID1 == item_id);
                temp_id = EMMCEdges[i].ID2;
                if (EMMCEdges[i].getW(sc_i)==INF100)
                    continue;
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
                if (EMMCEdges[i].getW(sc_i)==INF100)
                    continue;
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
            while (closed[pqueue.top().id] && !pqueue.empty()) {//if already closed
                pqueue.pop();
            }
            while (closed_r[pqueue_r.top().id] && !pqueue_r.empty()) {
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
    //function of LRU based BiDijkstra's algorithm (LRU+map+stxxl)
    Distance EMMCGraph::BiDijkstra_IO_new(int node_start, int node_end, HotPool4<VectorMCEdgesEMTuple4_BiDijkIO>& myHotPool){
         if(node_start==node_end){
             cout<<"Same source and target vertex id!"<<endl;
             return 0;
         }
         if(ifshow) {
             cout << "LRU+unordered_multimap+stxxl." << endl;
             ifshow= false;
         }
         //statstics of IO
         Stats stats_b, stats_e;
         stats_b.start();
         //Variables
         PriorityQueue pqueue(PQ_Pool);
         PriorityQueue2 pqueue_r(PQ_Pool2);
//         PriorityQueue2 pqueue_r(PQ_Pool2);


         NodeId item_id, temp_id;
         NodeId item_id_r, temp_id_r;
         Distance item_dis, temp_dis;
         Distance item_dis_r, temp_dis_r;
         vector<bool> closed(node_num, false); //flag vector of whether closed
         vector<Distance> cost(node_num, INF);   //vector of cost
//         vector<NodeId> pre(node_num, -1);       //vector of predecessor id
         vector<bool> closed_r(node_num, false); //flag vector of whether closed
         vector<Distance> cost_r(node_num, INF);   //vector of cost
//         vector<NodeId> pre_r(node_num, -1);       //vector of predecessor id
         Distance min_cost = INF;
         int index_begin, index_end;
         int terminate_id = -1;
         int partition_id, partition_id_r;

         bool flag_once = true;

         Timer tt1;
         tt1.start();
         //Initiation of start node
         cost[node_start] = 0;//cost of start node
         pqueue.push(VertexCost(node_start,0));
         cost_r[node_end] = 0;//cost of start node
         pqueue_r.push(VertexCost(node_end,0));
//         pqueue_r.push(make_pair(node_end,0));
//         pre[node_start] = node_start;
//         pre_r[node_end] = node_end;

//         string file_w = string(DataPath) + dataset + "/" + dataset + "_www.txt";
//         ofstream outFile(file_w, ios::out);
//         if (!outFile) {
//             cout << "File opening failed." << endl;
//             assert(outFile);
//         }
//        for(int i=1;i<4000000;++i){
//            pqueue_r.push(VertexCost(i,rand()));
//        }
//        pqueue_r.push(VertexCost(15738958,6267386));
//        int j=1;
//        while(!pqueue_r.empty()){
//            temp_id_r = pqueue_r.top().id;
//            if(temp_id_r == 15738958){
//                temp_dis_r = pqueue_r.top().cost;
//                cout<<"Test! Number "<<j<<" : "<<temp_id_r<<" (dis: "<<temp_dis_r<<") is topped from pqueue_r."<<endl;
//            }
//            pqueue_r.pop();
//            ++j;
//        }


         //Iteration
         while (!pqueue.empty() && !pqueue_r.empty()) {//for every node in pqueue

             //termination judgement
//             if (pqueue.top().cost + pqueue_r.top().second >= min_cost) {
             if (pqueue.top().cost + pqueue_r.top().cost >= min_cost) {//condition of termination pqueue_r.top().cost
//                 cout << node_start << " " << node_end << " " << sc_i << " "<<min_cost << endl;
                 break;
             }

             /// Forward searching
             item_id = pqueue.top().id;
             item_dis = pqueue.top().cost;
//             if(item_id == 15738958){
//                 cout<<"Forward! "<<item_id<<" (dis: "<<item_dis<<") is topped from pqueue."<<endl;
//             }
//             outFile << "Forward " << item_id <<" " << item_dis << endl;
             pqueue.pop();
             closed[item_id] = true;

             //relaxation
 //            set_visited.insert(item_id);//update the visited set for this criteria
             ++visited_number;//update #visited_number

             partition_id = node_to_cluster[item_id];
             if(myHotPool.clusterStatus[partition_id] == 'D'){//if in disk
                 GraphReadCluster4(myHotPool,partition_id);
                 ++p_num;
             }
             if(myHotPool.clusterStatus[partition_id] == 'I'){
                 /// read from EMMCEdges
//                 auto adjlist = myHotPool.HotPools[partition_id].equal_range(item_id);
//                     for(auto it = adjlist.first; it != adjlist.second; ++it){
                 for(auto it = myHotPool.HotPools[partition_id][item_id].begin(); it != myHotPool.HotPools[partition_id][item_id].end(); ++it){
                     temp_id = it->ID2;
                     if (it->getW(sc_i)==INF100)
                         continue;
                     temp_dis = item_dis + it->getW(sc_i);
                     if (!closed[temp_id]){//if not closed
                         if (cost[temp_id] > temp_dis) {//slack operation
                             cost[temp_id] = temp_dis;
//                             pre[temp_id] = item_id;
                             pqueue.push(VertexCost(temp_id, temp_dis));
                         }
                     }
                     if (closed_r[temp_id] && temp_dis + cost_r[temp_id] < min_cost) {
                         min_cost = temp_dis + cost_r[temp_id];
                         terminate_id = temp_id;
//                         cout<<"Forward! min_cost: "<<min_cost<<", terminate_id: "<<terminate_id<<endl;
                     }
                 }
                 closed[item_id] = true;
             }else if(myHotPool.clusterStatus[partition_id] == 'S'){
                 for (int i = get<1>(EMEdgesIndex[item_id]); i < get<2>(EMEdgesIndex[item_id]); ++i) {
                     assert(myHotPool.MCEdges_EM[i].ID1 == item_id);
                     temp_id = myHotPool.MCEdges_EM[i].ID2;
                     if (myHotPool.MCEdges_EM[i].getW(sc_i)==INF100)
                         continue;
                     temp_dis = item_dis + myHotPool.MCEdges_EM[i].getW(sc_i);
                     if (!closed[temp_id]){//if not closed
                         if (cost[temp_id] > temp_dis) {//slack operation
                             cost[temp_id] = temp_dis;
//                             pre[temp_id] = item_id;
                             pqueue.push(VertexCost(temp_id, temp_dis));
                         }
                     }
                     if (closed_r[temp_id] && temp_dis + cost_r[temp_id] < min_cost) {
                         min_cost = temp_dis + cost_r[temp_id];
                         terminate_id = temp_id;
//                         cout<<"Forward! min_cost: "<<min_cost<<", terminate_id: "<<terminate_id<<endl;
                     }
                 }
                 closed[item_id] = true;
             }

             /// Reverse searching
             item_id_r = pqueue_r.top().id;
             item_dis_r = pqueue_r.top().cost;
//             item_id_r = pqueue_r.top().first;
//             item_dis_r = pqueue_r.top().second;
//             if(item_id_r == 15738958){
//                 cout<<"Reverse! "<<item_id_r<<" (dis: "<<item_dis_r<<") is topped from pqueue_r."<<endl;
//             }
//             if(item_dis_r > 6267386 && flag_once){
//                 cout<<"Reverse! "<<item_id_r<<" (dis: "<<item_dis_r<<") is topped from pqueue_r."<<endl;
//                 flag_once = false;
//             }
//             outFile << "Reverse " << item_id_r <<" " << item_dis_r << endl;
             pqueue_r.pop();
             closed_r[item_id_r] = true;

             //relaxation
 //            set_visited.insert(item_id);//update the visited set for this criteria
             ++visited_number;//update #visited_number
             partition_id_r = node_to_cluster[item_id_r];
             if(myHotPool.clusterStatus[partition_id_r] == 'D'){//if in disk
                 GraphReadCluster4(myHotPool,partition_id_r);
                 ++p_num;
             }
             if(myHotPool.clusterStatus[partition_id_r] == 'I'){
                 /// read from EMMCEdges
//                 auto adjlist_r = myHotPool.HotPools[partition_id_r].equal_range(item_id_r);
//                 for(auto it = adjlist_r.first; it != adjlist_r.second; ++it){
                 for(auto it = myHotPool.HotPools[partition_id_r][item_id_r].begin(); it != myHotPool.HotPools[partition_id_r][item_id_r].end(); ++it){
                     temp_id_r = it->ID2;
//                     if(temp_id_r == 15738958 || temp_id_r == 15738957){
//                         cout<<"I! "<<temp_id_r<< " "<<it->getW(sc_i)<<" "<<closed_r[temp_id_r]<<endl;
//                     }
                     if (it->getW(sc_i)==INF100)
                         continue;
                     temp_dis_r = item_dis_r + it->getW(sc_i);
                     if (!closed_r[temp_id_r]){//if not closed
                         if (cost_r[temp_id_r] > temp_dis_r) {//slack operation
                             cost_r[temp_id_r] = temp_dis_r;
//                             pre_r[temp_id_r] = item_id_r;
                             pqueue_r.push(VertexCost(temp_id_r, temp_dis_r));
//                             pqueue_r.push(make_pair(temp_id_r, temp_dis_r));
//                             if(temp_id_r == 15738958){
//                                 cout<<"I! "<<temp_id_r<<" (dis: "<<temp_dis_r<<") is pushed into pqueue_r."<<endl;
//                             }
                         }
                     }
                     if (closed[temp_id_r] && temp_dis_r + cost[temp_id_r] < min_cost) {
                         min_cost = temp_dis_r + cost[temp_id_r];
                         terminate_id = temp_id_r;
//                         cout<<"Reverse! min_cost: "<<min_cost<<", terminate_id: "<<terminate_id<<endl;
                     }
                 }
                 closed_r[item_id_r] = true;
             }else if(myHotPool.clusterStatus[partition_id_r] == 'S'){
                 for (int i = get<1>(EMEdgesIndex[item_id_r]); i < get<2>(EMEdgesIndex[item_id_r]); ++i) {
                     assert(myHotPool.MCEdges_EM[i].ID1 == item_id_r);
                     temp_id_r = myHotPool.MCEdges_EM[i].ID2;
//                     if(temp_id_r == 15738958 || temp_id_r == 15738957){
//                         cout<<"S! "<<temp_id_r<< " "<<myHotPool.MCEdges_EM[i].getW(sc_i)<<" "<<closed_r[temp_id_r]<<endl;
//                     }
                     if (myHotPool.MCEdges_EM[i].getW(sc_i)==INF100)
                         continue;
                     temp_dis_r = item_dis_r + myHotPool.MCEdges_EM[i].getW(sc_i);
                     if (!closed_r[temp_id_r]){//if not closed
                         if (cost_r[temp_id_r] > temp_dis_r) {//slack operation
                             cost_r[temp_id_r] = temp_dis_r;
//                             pre_r[temp_id_r] = item_id_r;
                             pqueue_r.push(VertexCost(temp_id_r, temp_dis_r));
//                             pqueue_r.push(make_pair(temp_id_r, temp_dis_r));
//                             if(temp_id_r == 15738958){
//                                 cout<<"S! "<<temp_id_r<<" (dis: "<<temp_dis_r<<") is pushed into pqueue_r."<<endl;
//                             }
                         }
                     }
                     if (closed[temp_id_r] && temp_dis_r + cost[temp_id_r] < min_cost) {
                         min_cost = temp_dis_r + cost[temp_id_r];
                         terminate_id = temp_id_r;
//                         cout<<"Reverse! min_cost: "<<min_cost<<", terminate_id: "<<terminate_id<<endl;
                     }
                 }
                 closed_r[item_id_r] = true;
             }
             //deal with invalid top elements
             while (closed[pqueue.top().id]) {//if already closed
//                 if(pqueue.top().id == 15738958){
//                     cout<< pqueue.top().id <<" is popped accidentally!"<<endl;
//                 }
                 pqueue.pop();
             }
             while (closed_r[pqueue_r.top().id]) {
//                 if(pqueue_r.top().id == 15738958){
//                     cout<< pqueue_r.top().id <<" (dis: "<<pqueue_r.top().cost<<") is popped posthumously."<<endl;
//                 }
                 pqueue_r.pop();
             }
//             while (closed_r[pqueue_r.top().first]) {
//                 if(pqueue_r.top().first == 15738958){
//                     cout<< pqueue_r.top().first <<" (dis: "<<pqueue_r.top().second<<") is popped posthumously."<<endl;
//                 }
//                 pqueue_r.pop();
//             }
         }
//         outFile.close();
//         BiDij_getPath(pre,pre_r,node_start,node_end,terminate_id);//traverse the pre vector to get the shortest path
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
         while(!pqueue_r.empty()){
//             if(pqueue_r.top().id == 15738958){
//                 cout<< pqueue_r.top().id <<" is popped accidentally!"<<endl;
//             }
//             if(pqueue_r.top().first == 15738958){
//                 cout<< pqueue_r.top().first <<" is popped accidentally!"<<endl;
//             }
             pqueue_r.pop();
         }

         return min_cost;
     }
    //Entry for one-pass algorithms
    double EMMCGraph::MC_OnePass(const string& qtype,int algo_choice){
        string LongDis,MediumDis,ShortDis,RandomDis;
        LongDis = string(DataPath) + dataset + "/" + aggregateStrategy+"/"+dataset + "_OD_LongDis_" + to_string(PartitionSize) + ".txt";
        MediumDis = string(DataPath) + dataset + "/" + aggregateStrategy+"/"+dataset + "_OD_MediumDis_"+ to_string(PartitionSize) +".txt";
        ShortDis = string(DataPath) + dataset + "/" + aggregateStrategy+"/"+dataset + "_OD_ShortDis_"+ to_string(PartitionSize) +".txt";
        RandomDis = string(DataPath) +  dataset + "/" + aggregateStrategy+"/"+dataset + "_OD_Random_"+ to_string(PartitionSize) +".txt";

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
                cout << "Query "<<i<<" : "<<ID1 << " " << ID2<< endl;
                CommonInitiation_IO();//common initiation of each round
                OnePassInitiation_IO(algo_choice);

                if(algo_choice == BiMultiHops){
                    if(ifShortcut){
                        EM_MC_ShortcutSearch_Bi(ID1,ID2);
                    }else {
                        EM_MC_BiMultiHop(ID1, ID2);
                    }
                }else if(algo_choice == BiMultiHopsNoIO){
                    EM_MC_BiMultiHop_NoIO(ID1,ID2);
//                    EM_MC_BiMultiHop_NoIO_new(ID1,ID2);
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
                    if(mc_min_cost[j]<INF)
                        ave_cost[j] += mc_min_cost[j];
                }
                ave_visited += visited_number;
                ave_share += share_number;
                ave_io_time += io_time;
                ave_io_num += io_num;
                ave_p_num += p_num;
            }
        }
        else{
            for (int i = 0; i < run_times; ++i) {
                ID1 = ODpair[i].first;
                ID2 = ODpair[i].second;
                cout << "Query "<<i<<" : "<<ID1 << " " << ID2 << endl;

                CommonInitiation_IO();//common initiation of each round
                OnePassInitiation_IO(algo_choice);
                switch (algo_choice){
                    case OneHopNoIO:{
                        EM_MC_OneHop_NoIO(ID1, ID2); break;
//                        EM_MC_OneHop_NoIO_new(ID1, ID2); break;
                    }
                    case MultiHopsNoIO:{
                        EM_MC_MultiHop_NoIO(ID1, ID2); break;
//                        EM_MC_MultiHop_NoIO_new(ID1, ID2); break;
                    }
                    case OneHop:{
                        EM_MC_OneHop(ID1, ID2); break;
                    }
                    case MultiHops:{
                        if(ifShortcut){
                            EM_MC_ShortcutSearch(ID1,ID2);
                        }else{
                            EM_MC_MultiHop(ID1, ID2);
                        }
                        break;
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
    //One-hop algorithm without io optimization (stxxl+index)
    void EMMCGraph::EM_MC_OneHop_NoIO(int node_start, int node_end){
        if (node_start == node_end) {
            cout << "Source vertex and Target vertex are the same!!!" << endl;
            return;
        }
        if(ifshow){
            cout<<"stxxl+index."<<endl;
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
        VectorMCEdgesEMTuple_NoIO EMMCEdges;//the external vectors of different criteria
        vector<vector<Distance>> mc_cost(num_of_cri, vector<Distance>(node_num, INF));   //vector of cost
        //vector<vector<NodeId>> mc_pre(num_of_cri, vector<int>(node_num, -1));      //vector of predecessor id
        vector<vector<bool>> mc_closed(num_of_cri, vector<bool>(node_num));    //flag of vertex closing
        vector<int> cluster_to_cri(partition_number,0);//used to map which criterion read this partition
        set<int> set_remain;//set for storing the id of unshared criteria

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
                        if (EMMCEdges[i].getW(cri_i)==INF100)
                            continue;
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
    //function of reading partition
    void EMMCGraph::GraphMCReadCluster_NoIO(MCHotPool2 & mcHotPool,int target_p_id)
    {
        int evict_p_id = mcHotPool.LRU_IMCluster.put(target_p_id,1);
        //// evict old partition for new one
        if(evict_p_id != -1) {//if it is necessary to evict old partition
            mcHotPool.FlagIM[evict_p_id] = false;
            mcHotPool.HotPools[evict_p_id].clear();
            mcHotPool.HotPools[evict_p_id].shrink_to_fit();
        }
        //// read partition from original disk file to internal vector
        //check status of storage
        char filePath[300];
        int u, u_deg;
        int v, weight;
//        MCEdgeT temp_edge;

        strcpy(filePath, DataPath);
        strcat(filePath, dataset.c_str());
        strcat(filePath, "/");
        strcat(filePath, aggregateStrategy.c_str());
        strcat(filePath, "/Partitions/");
        strcat(filePath, dataset.c_str());
        strcat(filePath, "_Partition_");
        strcat(filePath, partMethod.c_str());
        strcat(filePath, "_");
        strcat(filePath, to_string(PartitionSize).c_str());
        strcat(filePath, "_");
        strcat(filePath, to_string(target_p_id).c_str());
        strcat(filePath, ".bin");
        Timer tt;
        tt.start();
        FILE* file = fopen(filePath, "rb");
        if(!file){
            cout<<"Failed to open file "<<filePath<<endl;
            exit(1);
        }
        assert(file != NULL);
        ReadBuffer read_io(file);
        bool flag_sizeUpdate = false;
        int temp_id = -1;
        MCEdgeCSR temp_edge;

        int temp;
        read_io.read(&temp);//cluster id
        read_io.read(&temp);//cluster size

        //read partition
        while(!read_io.isend) {//if not read out
            read_io.read(&u);//get node id
//            temp_edge.ID1 = u-1;
            if(u != temp_id){
                get<1>(EMEdgesIndex[u]) = mcHotPool.HotPools[target_p_id].size();
                if(temp_id != -1){
                    get<2>(EMEdgesIndex[temp_id]) = mcHotPool.HotPools[target_p_id].size();
                }
                temp_id = u;
            }

            read_io.read(&u_deg);//get node degree
            //vertex store in indexNode, edge store in e
            for (int i = 0; i < u_deg; ++i) {
                read_io.read(&v);//read end node
                temp_edge.ID2 = v;
                for (int j = 0; j < NUM_OF_CRITERIA; ++j) {
                    read_io.read(&weight);
                    if(j < num_of_cri){
                        temp_edge.putW(j,weight);
                    }
                }
                mcHotPool.HotPools[target_p_id].emplace_back(temp_edge);
            }
        }

        get<2>(EMEdgesIndex[u]) = mcHotPool.HotPools[target_p_id].size();

        io_num += read_io.io_number;//update #IO
        fclose(file);
        tt.stop();
        io_time += tt.GetRuntime() * 1000;
        cluster_to_IO[target_p_id] = read_io.io_number;

        mcHotPool.FlagIM[target_p_id] = true;
    }
    //Multi-hops algorithm without io optimization
    void EMMCGraph::EM_MC_MultiHop_NoIO(int node_start, int node_end){
        if (node_start == node_end) {
            cout << "Source vertex and Target vertex are the same!!!" << endl;
            return;
        }
        if(ifshow){
            cout<<"stxxl+index."<<endl;
            ifshow= false;
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
        set<int> set_remain;//set for storing the id of unshared criteria


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

                while (mc_closed[cri_i][EM_MC_PQueue[cri_i]->top().id]) {//if found in visited set && !EM_MC_PQueue[cri_i]->empty()
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
//                if(EM_MC_PQueue[cri_i]->empty()){
//                    mc_finished[cri_i] = true;
//                    set_cri.erase(cri_i);
//                    mc_min_cost[cri_i] = INF;
//                    cri_hops[cri_i] = INF;
//                    cri_pqueue.update(cri_i, INF);
//                    set_remain.erase(cri_i);
//                    continue;
//                }

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
                            if (EMMCEdges[i].getW(cri_i)==INF100)
                                continue;
                            if (mc_closed[cri_i][temp_id])//if closed
                                continue;
                            temp_dis = item_dis + EMMCEdges[i].getW(cri_i);
                            if (mc_cost[cri_i][temp_id] > temp_dis) {//slack operation
                                mc_cost[cri_i][temp_id] = temp_dis;
                                //mc_pre[cri_i][temp_id] = item_id;
                                EM_MC_PQueue[cri_i]->push(VertexCost(temp_id, temp_dis));
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
                    /// read by synchronizing the hops
                    temp_cluster_id = cri_to_cluster[cri_pqueue.top_value()];
                    assert(get<1>(EMEdgesIndex[cluster_to_node[temp_cluster_id][0]]) == -1);
                    ReadClusterToExternalVectorMC(temp_cluster_id, EMMCEdges);
                    ++p_num;
                }
                jump1: "there is a jump";
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
                    assert(get<1>(EMEdgesIndex[cluster_to_node[temp_cluster_id][0]]) == -1);
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
    //BiMulti-hops algorithm without io optimization: 2022-03-23 version 1
    void EMMCGraph::EM_MC_BiMultiHop_NoIO(int node_start, int node_end) {
        if (node_start == node_end) {
            cout << "Source vertex and Target vertex are the same!!!" << endl;
            return;
        }
        if(ifshow){
            cout<<"stxxl+index."<<endl;
            ifshow= false;
        }
        //statstics of IO
        Stats stats_b, stats_e;
        stats_b.start();
        //variables
        set<int> set_remain;//set for storing the id of unshared criteria
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
                    if(EM_MC_PQueue_r[cri_i]->empty()){
                        mc_finished[cri_i] = true;
                        set_cri.erase(cri_i);
                        cri_hops[cri_i] = INF;
                        cri_pqueue.update(cri_i, INF);
                        cri_hops_r[cri_i] = INF;
                        cri_pqueue_r.update(cri_i, INF);
                        goto jump1;
                    }
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
//                            cout<<"Minimal cost of criterion " << cri_i << ": "<<mc_min_cost[cri_i]<<endl;
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
                            if (EMMCEdges[i].getW(cri_i)==INF100)
                                continue;
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
                            if (EMMCEdges[i].getW(cri_i)==INF100)
                                continue;
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

                        while (mc_closed_r[cri_i][EM_MC_PQueue_r[cri_i]->top().id]  && !EM_MC_PQueue_r[cri_i]->empty()) {
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
                }
                set_remain.erase(cri_i);
                // judge reading
                if (set_remain.empty() && !flag_mcFinished) {//if all criteria are processed and not all criteria are finished
                    /// read by synchronizing the hops
                    temp_cluster_id = cri_to_cluster[cri_pqueue.top_value()];
                    if (cluster_to_bi[temp_cluster_id] == 0 || get<1>(EMEdgesIndex[cluster_to_node[temp_cluster_id][0]]) == -1) {//partitions_info[temp_cluster_id]
                        ReadClusterToExternalVectorMC(temp_cluster_id, EMMCEdges);
                        cluster_to_bi[temp_cluster_id] = FORWARD;
                        ++p_num;
                    }
                    temp_cluster_id_r = cri_to_cluster_r[cri_pqueue_r.top_value()];
                    if (cluster_to_bi[temp_cluster_id_r] == 0 || get<1>(EMEdgesIndex[cluster_to_node[temp_cluster_id_r][0]]) == -1) {
                        if (temp_cluster_id != temp_cluster_id_r) {//if different
                            ReadClusterToExternalVectorMC(temp_cluster_id_r, EMMCEdges);
                            cluster_to_bi[temp_cluster_id_r] = REVERSE;
                            ++p_num;
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
                        if (cluster_to_bi[temp_cluster_id] == 0 || get<1>(EMEdgesIndex[cluster_to_node[temp_cluster_id][0]]) == -1) {
//                        ReadClusterToExternalVectorMC(temp_cluster_id, EMMCEdges);
//                        GraphMCReadCluster_NoIO(myHotPool,temp_cluster_id);
                            GraphReadCluster4(myHotPool,temp_cluster_id);
                            cluster_to_bi[temp_cluster_id] = FORWARD;
                            ++p_num;
                        }
                    }else{
                        if (cluster_to_bi[temp_cluster_id_r] == 0 || get<1>(EMEdgesIndex[cluster_to_node[temp_cluster_id_r][0]]) == -1) {
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
                                for(auto it=cluster_to_node[partition_id].begin();it!=cluster_to_node[partition_id].end();++it){
                                    get<1>(EMEdgesIndex[*it]) = -1;//hard delete
                                    get<2>(EMEdgesIndex[*it]) = -1;
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
                for(auto it=cluster_to_node[evict_p_id].begin();it!=cluster_to_node[evict_p_id].end();++it){
                    get<1>(EMEdgesIndex[*it]) = -1;//hard delete
                    get<2>(EMEdgesIndex[*it]) = -1;
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
        char filePath[300];
        int u, u_deg;
        int v, weight;
        MCEdgeT temp_edge;

        strcpy(filePath, DataPath);
        strcat(filePath, dataset.c_str());
        strcat(filePath, "/");
        strcat(filePath, aggregateStrategy.c_str());
        strcat(filePath, "/Partitions/");
        strcat(filePath, dataset.c_str());
        strcat(filePath, "_Partition_");
        strcat(filePath, partMethod.c_str());
        strcat(filePath, "_");
        strcat(filePath, to_string(PartitionSize).c_str());
        strcat(filePath, "_");
        strcat(filePath, to_string(target_p_id).c_str());
        strcat(filePath, ".bin");
        Timer tt;
        tt.start();
        FILE* file = fopen(filePath, "rb");
        if(!file){
            cout<<"Failed to open file "<<filePath<<endl;
            exit(1);
        }
        assert(file != NULL);
        ReadBuffer read_io(file);
        bool flag_sizeUpdate = false;
        int temp_id = -1;
        int temp;
        read_io.read(&temp);//cluster id
        read_io.read(&temp);//cluster size


        if(get<1>(mcHotPool.clusterStatus[target_p_id]) == 0){//if this is the first read
            flag_sizeUpdate = true;
            //read partition
            while(!read_io.isend) {//if not read out
                read_io.read(&u);//get node id
                temp_edge.ID1 = u;
                if(u != temp_id){
                    get<1>(EMEdgesIndex[u]) = mcHotPool.HotPools[target_p_id].size();
                    if(temp_id != -1){
                        get<2>(EMEdgesIndex[temp_id]) = mcHotPool.HotPools[target_p_id].size();
                    }
                    temp_id = u;
                }

//                if(vertex_cri[u-1] == -1){//if not found
                assert(vertex_cri[u] == -1);
                    vertex_cri[u] = num_of_cri;
//                }
                read_io.read(&u_deg);//get node degree
                //vertex store in indexNode, edge store in e
                for (int i = 0; i < u_deg; ++i) {
                    read_io.read(&v);//read end node
                    temp_edge.ID2 = v;
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
        }
        else{//if this is not the first read
            //read partition
            while(!read_io.isend) {//if not read out
                read_io.read(&u);//get node id
                temp_edge.ID1 = u;
                if(u != temp_id){
                    get<1>(EMEdgesIndex[u]) = mcHotPool.HotPools[target_p_id].size();
                    if(temp_id != -1){
                        get<2>(EMEdgesIndex[temp_id]) = mcHotPool.HotPools[target_p_id].size();
                    }
                    temp_id = u;
                }

                if(vertex_cri[u] == -1){//if not found
                    vertex_cri[u] = num_of_cri;
                }
                read_io.read(&u_deg);//get node degree
                //vertex store in indexNode, edge store in e
                for (int i = 0; i < u_deg; ++i) {
                    read_io.read(&v);//read end node
                    temp_edge.ID2 = v;
                    for (int j = 0; j < NUM_OF_CRITERIA; ++j) {
                        read_io.read(&weight);
                        if(j < num_of_cri){
                            temp_edge.putW(j,weight);
                        }
                    }
                    mcHotPool.HotPools[target_p_id].emplace_back(temp_edge);

                }
                if (!get<0>(EMEdgesIndex[u]))//if false, if it has not been evicted
                    get<2>(mcHotPool.clusterStatus[target_p_id]) += u_deg;
            }
            io_num += 1;//update #IO
        }


        get<2>(EMEdgesIndex[u]) = mcHotPool.HotPools[target_p_id].size();
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
                        cout<<"Minimal cost of criterion " << cri_i << ": "<<mc_min_cost[cri_i]<<endl;
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
                                    for(auto it=cluster_to_node[partition_id].begin();it!=cluster_to_node[partition_id].end();++it){
                                        get<1>(EMEdgesIndex[*it]) = -1;//hard delete
                                        get<2>(EMEdgesIndex[*it]) = -1;
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
    //Function of Bidirectional Multi-Hop algorithm, version 2
    void EMMCGraph::EM_MC_BiMultiHop(int node_start, int node_end) {
        if (node_start == node_end) {
            cout << "Source vertex and Target vertex are the same!!!" << endl;
            return;
        }
        cout<<"BMHP"<<endl;
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
                        cout<<"T2 Minimal cost of criterion "<<cri_i<<": "<<mc_min_cost[cri_i]<<endl;
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
                                    for(auto it=cluster_to_node[partition_id].begin();it!=cluster_to_node[partition_id].end();++it){
                                        get<2>(EMEdgesIndex_Bi[*it]) = -1;//hard delete
                                        get<3>(EMEdgesIndex_Bi[*it]) = -1;
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
                        cout<<"T2 Minimal cost of criterion "<<cri_i<<": "<<mc_min_cost[cri_i]<<endl;
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
                            cout << "Wrong partition status!!!" << endl; break;
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
                                    for(auto it=cluster_to_node[partition_id_r].begin();it!=cluster_to_node[partition_id_r].end();++it){
                                        get<2>(EMEdgesIndex_Bi[*it]) = -1;//hard delete
                                        get<3>(EMEdgesIndex_Bi[*it]) = -1;
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
                    cout<<"Minimal cost of criterion "<<cri_i<<": "<<mc_min_cost[cri_i]<<endl;
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
                for(auto it=cluster_to_node[evict_p_id].begin();it!=cluster_to_node[evict_p_id].end();++it){
                    get<2>(EMEdgesIndex_Bi[*it]) = -1;//hard delete
                    get<3>(EMEdgesIndex_Bi[*it]) = -1;
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
        char filePath[300];
        int u = 0, u_deg;
        int v, weight;
        MCEdgeT temp_edge;

        strcpy(filePath, DataPath);
        strcat(filePath, dataset.c_str());
        strcat(filePath, "/");
        strcat(filePath, aggregateStrategy.c_str());
        strcat(filePath, "/Partitions/");
        strcat(filePath, dataset.c_str());
        strcat(filePath, "_Partition_");
        strcat(filePath, partMethod.c_str());
        strcat(filePath, "_");
        strcat(filePath, to_string(PartitionSize).c_str());
        strcat(filePath, "_");
        strcat(filePath, to_string(target_p_id).c_str());
        strcat(filePath, ".bin");

        Timer tt;
        tt.start();
        FILE* file = fopen(filePath, "rb");
        if(!file){
            cout<<"Failed to open file "<<filePath<<endl;
            exit(1);
        }
        assert(file != NULL);
        ReadBuffer read_io(file);
        bool flag_sizeUpdate = false;
        int temp_id = -1;
        int temp;
        read_io.read(&temp);//cluster id
        read_io.read(&temp);//cluster size

        if(get<1>(mcHotPool.clusterStatus[target_p_id]) == 0){//if this is the first read, judged by the original size of the cluster
            flag_sizeUpdate = true;
            while(!read_io.isend) {//if not read out
                read_io.read(&u);//get node id
                temp_edge.ID1 = u;
                if(u != temp_id){
                    get<2>(EMEdgesIndex_Bi[u]) = mcHotPool.HotPools[target_p_id].size();
                    if(temp_id != -1){
                        get<3>(EMEdgesIndex_Bi[temp_id]) = mcHotPool.HotPools[target_p_id].size();
                    }
                    temp_id = u;
                }

                assert(vertex_cri_Bi[u].first == -1);
                assert(vertex_cri_Bi[u].second == -1);
                vertex_cri_Bi[u].first = num_of_cri;
                vertex_cri_Bi[u].second = num_of_cri;

                read_io.read(&u_deg);//get node degree

                for (int i = 0; i < u_deg; ++i) {
                    read_io.read(&v);//read end node
                    temp_edge.ID2 = v;
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
        }
        else{//if it is not the first read, read all
            cout<<"Partition. The second read for partition "<<target_p_id<<endl;
            mcHotPool.HotPools[target_p_id].clear();
            if(direction == FORWARD) {//if the partition now is read by forward search
                if (cluster_to_bi[target_p_id] == REVERSE) {//if it has been read by the other search, update the edge number in memory
                    get<2>(mcHotPool.clusterStatus[target_p_id]) = 0;
                }
                while (!read_io.isend) {//if not read out
                    read_io.read(&u);//get node id
                    temp_edge.ID1 = u;

                    if (u != temp_id) {
                        get<2>(EMEdgesIndex_Bi[u]) = mcHotPool.HotPools[target_p_id].size();

                        if (temp_id != -1) {
                            get<3>(EMEdgesIndex_Bi[temp_id]) = mcHotPool.HotPools[target_p_id].size();
                        }
                        temp_id = u;
                    }

                    read_io.read(&u_deg);//get node degree

                    for (int i = 0; i < u_deg; ++i) {
                        read_io.read(&v);//read end node
                        temp_edge.ID2 = v;
                        for (int j = 0; j < NUM_OF_CRITERIA; ++j) {
                            read_io.read(&weight);// read edge weight
                            if (j < num_of_cri) {
                                temp_edge.putW(j, weight);
                            }
                        }
                        mcHotPool.HotPools[target_p_id].emplace_back(temp_edge);
                    }
                    if (cluster_to_bi[target_p_id] == REVERSE) {//if it has been read by the other search
                        if (!get<0>(EMEdgesIndex_Bi[u]))//if false, if it has not been evicted
                            get<2>(mcHotPool.clusterStatus[target_p_id]) += u_deg;
                    }
                }
            }else if(direction == REVERSE) {//if the partition now is read by reverse search
                if(cluster_to_bi[target_p_id] == FORWARD) {//if it has been read by the other search
                    get<2>(mcHotPool.clusterStatus[target_p_id]) = 0;
                }
                while(!read_io.isend) {//if not read out
                    read_io.read(&u);//get node id
                    temp_edge.ID1 = u;

                    if(u != temp_id){
                        get<2>(EMEdgesIndex_Bi[u]) = mcHotPool.HotPools[target_p_id].size();

                        if(temp_id != -1){
                            get<3>(EMEdgesIndex_Bi[temp_id]) = mcHotPool.HotPools[target_p_id].size();
                        }
                        temp_id = u;
                    }

                    read_io.read(&u_deg);//get node degree

                    for (int i = 0; i < u_deg; ++i) {
                        read_io.read(&v);//read end node
                        temp_edge.ID2 = v;
                        for (int j = 0; j < NUM_OF_CRITERIA; ++j) {
                            read_io.read(&weight);
                            if(j < num_of_cri){
                                temp_edge.putW(j,weight);
                            }
                        }
                        mcHotPool.HotPools[target_p_id].emplace_back(temp_edge);
                    }
                    if(cluster_to_bi[target_p_id] == FORWARD){//if it has been read by the other search
                        if(!get<1>(EMEdgesIndex_Bi[u]))//if false, if it has not been evicted
                            get<2>(mcHotPool.clusterStatus[target_p_id]) += u_deg;
                    }
                }
            }
            io_num += 1;//update #IO
        }

        get<3>(EMEdgesIndex_Bi[u]) = mcHotPool.HotPools[target_p_id].size();
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
    //new
    void EMMCGraph::DijkstraSC(int bid, int cri_i, vector<NodeId>& borderVector, vector<NodeId>& testBorder, vector<int>& testResult, unordered_set<NodeId>& newNeighbor){
        int node_start=borderVector[bid];
        int partition_id=node_to_cluster[node_start];
//    benchmark::heap<2, int, Distance> pqueue(node_num);
        priority_queue<pair<int,int>, vector<pair<int,int>>, IntPairCompareLess> pqueue;

        unordered_set<NodeId> remainedBorder; remainedBorder.clear();
        remainedBorder.insert(testBorder.begin(),testBorder.end()) ;

        //variables
//    vector<int> cost(node_num,INF);
//    vector<bool> closed(node_num,false);
//    vector<bool> addSC(node_num,false);//whether to add the shortcut
//    vector<bool> stained(node_num,false);//whether pass a boundary vertex
//    vector<NodeId> pre(node_num,-1);
        unordered_map<NodeId,int> cost; cost.clear();
        unordered_map<NodeId,bool> closed; closed.clear();
        unordered_map<NodeId,NodeId> pre; pre.clear();
        unordered_map<NodeId,bool> addSC; addSC.clear();//whether to add the shortcut
        unordered_map<NodeId,bool> stained; stained.clear();//whether pass a boundary vertex
        for(auto it=cluster_to_node[partition_id].begin();it!=cluster_to_node[partition_id].end();++it){
            int id=*it;
            cost.insert({id,INF});
            closed.insert({id,false});
            pre.insert({id,id});
            addSC.insert({id,false});
            stained.insert({id,false});
        }

        NodeId item_id;
        bool flag_finish = false;
        NodeId temp_id;
        Distance temp_dis;
        Distance item_dis;
        int temp_number;
        int temp_degree = 0;
        int cri_empty = 0;

        Timer tt;
        tt.start();

        //Initiation of start node
        cost[node_start] = 0;//cost of start node
//    pqueue.update(node_start,0);
        pqueue.push(make_pair(node_start,0));
        pre[node_start] = node_start;

        //Iteration
        while (!flag_finish && !pqueue.empty()) {//for every node in pqueue !mc_finished[0] || !mc_finished[1] || !mc_finished[2]
//        pqueue.extract_min(item_id, item_dis);
            while(closed[pqueue.top().first]){//if true
                pqueue.pop();
                if(pqueue.empty()){
                    goto jump1;
                }
            }
            item_id=pqueue.top().first; item_dis= pqueue.top().second;
            if(item_dis!=cost[item_id]){
                cout<<"Distance inconsistent! "<<item_id<<" "<<item_dis<<" "<<cost[item_id]<<endl;
                exit(1);
            }
            pqueue.pop();

            if(node_to_cluster[item_id]!=partition_id){
                cout<<"Wrong item_id!"<<endl; exit(1);
            }

            if (BoundaryTag[item_id] != 0){//if it is boundary vertex
                if(item_id!=node_start){
                    if(!stained[pre[item_id]]){//if precedent is not stained
                        stained[item_id] = true;
                        addSC[item_id] = true;
                    }else{//if precedent is stained
                        stained[item_id] = true;
                    }

                    if (remainedBorder.find(item_id) != remainedBorder.end()) {//if reach
                        remainedBorder.erase(item_id);
//                        cout<<"Shortcut of neighbor "<<item_id<<": "<<cost[item_id]<<endl;
                        if(remainedBorder.empty()){
                            flag_finish = true;
//                            cout<<"Finish with testBorder is empty! "<< node_start << endl;
                            goto jump1;
                        }
                    }
                }

            }else{//if it is not boundary vertex
                if(stained[pre[item_id]]){//if precedent is stained
                    stained[item_id] = true;
                }
            }


            //set closed
            closed[item_id] = true;
            temp_degree = 0;
//            if(item_id==11783){
//                cout<<item_id<<" "<<get<1>(EMEdgesIndex[item_id])<<" "<<get<2>(EMEdgesIndex[item_id])<<endl;
//            }

            for (int i = get<1>(EMEdgesIndex[item_id]); i < get<2>(EMEdgesIndex[item_id]); ++i) {
                ++temp_degree;
                temp_id = NeighborsMC[partition_id][i].ID2;
                Distance temp_w = NeighborsMC[partition_id][i].getW(cri_i);
//        for (int i = 0; i < MCEdgesT[item_id].size(); ++i) {
//            ++temp_degree;
//            temp_id = MCEdgesT[item_id][i].ID2;
//            Distance temp_w = MCEdgesT[item_id][i].getW(cri_i);

                if(node_to_cluster[temp_id] != partition_id){
                    continue;
                }
                if (temp_w==INF100)
                    continue;
                if (closed[temp_id])//if closed
                    continue;
                temp_dis = item_dis + temp_w;

                if (cost[temp_id] > temp_dis) {//slack operation
                    cost[temp_id] = temp_dis;
                    pre[temp_id] = item_id;
//                pqueue.update(temp_id, temp_dis);
                    pqueue.push(make_pair(temp_id,temp_dis));
                }else if(cost[temp_id] == temp_dis){
                    if(stained[item_id]){//if precedent is stained
                        pre[temp_id] = item_id;
                    }
                }
            }
        }
        jump1:
        int stainN=0;
        for(int i=0;i<testBorder.size();++i){
            int id=testBorder[i];
            if(cost.find(id)==cost.end()){
                cout<<"Not found "<<id<<" "<< node_to_cluster[id]<<" "<<partition_id<<endl;
                exit(1);
            }
            testResult[i]=cost[id];
            if(addSC[id]){
//                cout<<"Boundary vertex "<<id<<" "<<testResult[i]<<endl;
                newNeighbor.insert(id);
                stainN++;
            }
        }
//        cout<<node_start<<" 's Boundary shortcut number: "<<stainN<<" "<<borderVector.size()<<", criterion "<< cri_i <<", partition "<<partition_id<<endl;
        tt.stop();
    }

    //new
    void EMMCGraph::ShortcutVertex(int bid, int pid, vector<vector<MCEdgeT>>& parResult){
        int ID=borderVectors[pid][bid];
        bool flag_full=true;
        if(BoundaryTag[ID]==2){
            flag_full=true;
        }else if(BoundaryTag[ID]==1){
            flag_full=false;
        }else{
            cout<<"Wrong BoundaryTag! "<<int(BoundaryTag[ID])<<endl; exit(1);
        }

        vector<NodeId> testBorder; testBorder.clear();
        unordered_set<NodeId> newNeighbor;
        int temp_id;

        for (int i = get<1>(EMEdgesIndex[ID]); i < get<2>(EMEdgesIndex[ID]); ++i) {
            temp_id = NeighborsMC[pid][i].ID2;

            if(node_to_cluster[temp_id]==pid){//if neighbor is in the same partition
                if(flag_full){
                    testBorder.push_back(temp_id);
//                parResult[bid].push_back(NeighborsMC[pid][i]);
//                newNeighbor.insert(temp_id);
                }

            }else{//if neighbor is in other partitions
                parResult[bid].push_back(NeighborsMC[pid][i]);
                newNeighbor.insert(temp_id);
            }
        }


        if(!flag_full){//if it is half-connected boundary vertex
            testBorder=borderVectors[pid];
        }


        vector<vector<int>> testResults;
        Timer tt;
        if(!testBorder.empty()){//not empty
            vector<int> temp_v(testBorder.size(),INF);
            testResults.assign(NUM_OF_CRITERIA,temp_v);
            parResult[bid].reserve(testBorder.size());
            tt.start();
            for(int cri=0;cri<NUM_OF_CRITERIA;++cri){
                DijkstraSC(bid,cri,borderVectors[pid],testBorder,testResults[cri], newNeighbor);
//            DijkstraPart(bid,cri,borderVectors[pid]);
            }
            tt.stop();
//            cout<<"Time of "<<ID<<": "<<tt.GetRuntime()<<" s."<<endl;
        }
        else{//if empty
            if(!flag_full){//if half-connected vertex
                cout<<"Wrong for half-connected vertex! BorderSet is empty! "<<flag_full<<" "<<ID<<endl; exit(1);
            }
        }

        if(newNeighbor.empty()){
            cout<<"Empty neighbor set! "<<ID<<endl; exit(1);
        }

        MCEdgeT tempE; tempE.ID1=ID;
        for(int i=0;i<testBorder.size();++i){
            if(newNeighbor.find(testBorder[i])!=newNeighbor.end()){//if it is new neighbor
                tempE.ID2=testBorder[i];
                for(int j=0;j<testResults.size();++j){
                    if(testResults[j][i] <= 0){
                        cout<<"Wrong shortcut edge! "<<ID<<"("<<int(BoundaryTag[ID])<<") "<<tempE.ID2<<"("<<int(BoundaryTag[tempE.ID2])<<") "<<testResults[j][i]<<" "<< testBorder.size()<<endl;
                        exit(1);
                    }
                    if(testResults[j][i] == INF){
//                        cout<<"Shortcut weight between "<<ID<<" and "<<tempE.ID2<<" is "<<INF<<" , criterion "<<j<<endl;
                        tempE.putW(j,INF100);
                    }else{
                        tempE.putW(j,testResults[j][i]);
                    }
                }
//                parResult[bid][i]=tempE;
                parResult[bid].emplace_back(tempE);
            }
        }
    }

    void EMMCGraph::ShortcutVertexV(vector<int>& p, int pid, vector<vector<MCEdgeT>>& parResult){
        int bid;
        for(int i=0;i<p.size();++i){
            bid=p[i];
            ShortcutVertex(bid, pid, parResult);
        }
    }
    //function of reading the original partition graph from disk
    void EMMCGraph::ReadPartition(int pid, vector<MCEdgeT>& NeighborsP){
        //check status of storage
        char filePath[300];
        int u, u_deg;
        int v, weight;
        MCEdgeT temp_edge;

        strcpy(filePath, DataPath);
        strcat(filePath, dataset.c_str());
        strcat(filePath, "/");
        strcat(filePath, aggregateStrategy.c_str());
        strcat(filePath, "/Partitions/");
        strcat(filePath, dataset.c_str());
        strcat(filePath, "_Partition_");
        strcat(filePath, partMethod.c_str());
        strcat(filePath, "_");
        strcat(filePath, to_string(PartitionSize).c_str());
        strcat(filePath, "_");
        strcat(filePath, to_string(pid).c_str());
        strcat(filePath, ".bin");
        Timer tt;
        tt.start();
        FILE* file = fopen(filePath, "rb");
        if(!file){
            cout<<"Failed to open file "<<filePath<<endl;
            exit(1);
        }
        assert(file != NULL);
        ReadBuffer read_io(file);
        bool flag_sizeUpdate = false;
        int temp_id = -1;

        int temp;
        read_io.read(&temp);//cluster id
        read_io.read(&temp);//cluster size

        flag_sizeUpdate = true;
        //read partition
        while(!read_io.isend) {//if not read out
            read_io.read(&u);//get node id
            temp_edge.ID1 = u;

//            if(temp_edge.ID1==11783){
//                cout<<"read "<<temp_edge.ID1;
//            }
            if(u != temp_id){
                get<1>(EMEdgesIndex[u]) = NeighborsP.size();
                if(temp_id != -1){
                    get<2>(EMEdgesIndex[temp_id]) = NeighborsP.size();
                }
                temp_id = u;
            }

            read_io.read(&u_deg);//get node degree
            //vertex store in indexNode, edge store in e
            for (int i = 0; i < u_deg; ++i) {
                read_io.read(&v);//read end node
                temp_edge.ID2 = v;
                for (int j = 0; j < NUM_OF_CRITERIA; ++j) {
                    read_io.read(&weight);
//                    if(j < num_of_cri){
                        temp_edge.putW(j,weight);
//                    }
                }
                NeighborsP.emplace_back(temp_edge);

            }

        }
        fclose(file);
        get<2>(EMEdgesIndex[u]) = NeighborsP.size();
        io_num += read_io.io_number;//update #IO
        tt.stop();
        io_time += tt.GetRuntime() * 1000;//ms
    }
    void EMMCGraph::ReadPartitionBi(int pid, vector<MCEdgeT>& NeighborsP){
        char filePath[300];
        int u = 0, u_deg;
        int v, weight;
        MCEdgeT temp_edge;

        strcpy(filePath, DataPath);
        strcat(filePath, dataset.c_str());
        strcat(filePath, "/");
        strcat(filePath, aggregateStrategy.c_str());
        strcat(filePath, "/Partitions/");
        strcat(filePath, dataset.c_str());
        strcat(filePath, "_Partition_");
        strcat(filePath, partMethod.c_str());
        strcat(filePath, "_");
        strcat(filePath, to_string(PartitionSize).c_str());
        strcat(filePath, "_");
        strcat(filePath, to_string(pid).c_str());
        strcat(filePath, ".bin");

        Timer tt;
        tt.start();
        FILE* file = fopen(filePath, "rb");
        if(!file){
            cout<<"Failed to open file "<<filePath<<endl;
            exit(1);
        }
        assert(file != NULL);
        ReadBuffer read_io(file);
        bool flag_sizeUpdate = false;
        int temp_id = -1;
        int temp;
        read_io.read(&temp);//cluster id
        read_io.read(&temp);//cluster size

        flag_sizeUpdate = true;
        while(!read_io.isend) {//if not read out
            read_io.read(&u);//get node id
            temp_edge.ID1 = u;
            if(u != temp_id){
                get<2>(EMEdgesIndex_Bi[u]) = NeighborsP.size();
                if(temp_id != -1){
                    get<3>(EMEdgesIndex_Bi[temp_id]) = NeighborsP.size();
                }
                temp_id = u;
            }

            assert(vertex_cri_Bi[u].first == -1);
            assert(vertex_cri_Bi[u].second == -1);
//            vertex_cri_Bi[u].first = num_of_cri;
//            vertex_cri_Bi[u].second = num_of_cri;

            read_io.read(&u_deg);//get node degree

            for (int i = 0; i < u_deg; ++i) {
                read_io.read(&v);//read end node
                temp_edge.ID2 = v;
                for (int j = 0; j < NUM_OF_CRITERIA; ++j) {
                    read_io.read(&weight);
                    if(j < num_of_cri){
                        temp_edge.putW(j,weight);
                    }
                }
                NeighborsP.emplace_back(temp_edge);
            }
        }

        fclose(file);
        get<3>(EMEdgesIndex_Bi[u]) = NeighborsP.size();
        io_num += read_io.io_number;//update #IO
        tt.stop();
        io_time += tt.GetRuntime() * 1000;//ms
    }
    void EMMCGraph::ReadPartitionMCEdges(int pid){
        //check status of storage
        char filePath[300];
        int u, u_deg;
        int v, weight;
        MCEdgeT temp_edge;

        strcpy(filePath, DataPath);
        strcat(filePath, dataset.c_str());
        strcat(filePath, "/");
        strcat(filePath, aggregateStrategy.c_str());
        strcat(filePath, "/Partitions/");
        strcat(filePath, dataset.c_str());
        strcat(filePath, "_Partition_");
        strcat(filePath, partMethod.c_str());
        strcat(filePath, "_");
        strcat(filePath, to_string(PartitionSize).c_str());
        strcat(filePath, "_");
        strcat(filePath, to_string(pid).c_str());
        strcat(filePath, ".bin");
        Timer tt;
        tt.start();
        FILE* file = fopen(filePath, "rb");
        if(!file){
            cout<<"Failed to open file "<<filePath<<endl;
            exit(1);
        }
        assert(file != NULL);
        ReadBuffer read_io(file);
        bool flag_sizeUpdate = false;
        int temp_id = -1;

        int temp;
        read_io.read(&temp);//cluster id
        read_io.read(&temp);//cluster size

        flag_sizeUpdate = true;
        //read partition
        while(!read_io.isend) {//if not read out
            read_io.read(&u);//get node id
            temp_edge.ID1 = u;

//            if(temp_edge.ID1==11783){
//                cout<<"read "<<temp_edge.ID1;
//            }
//            if(u != temp_id){
//                get<1>(EMEdgesIndex[u]) = NeighborsP.size();
//                if(temp_id != -1){
//                    get<2>(EMEdgesIndex[temp_id]) = NeighborsP.size();
//                }
//                temp_id = u;
//            }

            read_io.read(&u_deg);//get node degree
            //vertex store in indexNode, edge store in e
            for (int i = 0; i < u_deg; ++i) {
                read_io.read(&v);//read end node
                temp_edge.ID2 = v;
                for (int j = 0; j < NUM_OF_CRITERIA; ++j) {
                    read_io.read(&weight);
                    if(j < num_of_cri){
                        temp_edge.putW(j,weight);
                    }
                }
                MCEdgesT[u].emplace_back(temp_edge);

            }

        }
        fclose(file);
//        get<2>(EMEdgesIndex[u]) = NeighborsP.size();
        io_num += read_io.io_number;//update #IO
        tt.stop();
        io_time += tt.GetRuntime() * 1000;//ms
    }

    //function of reading graph data into MCEdges
    void EMMCGraph::ReadMCEdges(const string& filename){
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
//            node_num=stoi(re1[0]), edge_num=stoul(re1[1]);
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
//        mc_criteria.clear();
        int tempNum_criteria=stoi(re1[0]);
        if(re1.size()==tempNum_criteria+1){
//            for(int i=0;i<tempNum_criteria;++i){
//                temp_str = re1[i+1];
//                cout << temp_str << " " ;
//                mc_criteria.push_back(temp_str);
//            }
        }else{
            cout<<"Wrong input! re.size: "<<re1.size()<<". " <<line<<endl;
            exit(1);
        }
//        cout<<endl;


        while(getline(inFile,line)){
            if(line.empty())
                continue;
            re1.clear();
            boost::split(re1,line,boost::is_any_of(" \t"));
            if(re1.size()==tempNum_criteria+2){
                ID1=stoi(re1[0]), ID2=stoi(re1[1]);
            }
            else{
                cout<<"Wrong line. "<< line<<endl;
                exit(1);
            }


            if(ID1>=0 && ID1<node_num && ID2>=0 && ID2<node_num && ID1!=ID2){
                set_A.insert(ID1), set_A.insert(ID2);
                temp_w.clear();
                for(int i=0;i<tempNum_criteria;++i) {//get edge weights
                    weight=stoi(re1[i+2]);
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


        tt.stop();
        cout << "The time used for data loading:" << tt.GetRuntime() << " s." << endl;
        inFile.close();
    }

    void EMMCGraph::WritePartitionShortcut(int pid, vector<NodeId>& borderVector, vector<vector<MCEdgeT>>& scResult, bool ifbin){

        if(ifbin){
            char file_bin_path[300];

            strcpy(file_bin_path, DataPath);
            strcat(file_bin_path, dataset.c_str());
            strcat(file_bin_path, "/");
            strcat(file_bin_path, aggregateStrategy.c_str());
            strcat(file_bin_path, "/Partitions/");
            strcat(file_bin_path, dataset.c_str());
            strcat(file_bin_path, "_Shortcut_");
            strcat(file_bin_path, partMethod.c_str());
            strcat(file_bin_path, "_");
            strcat(file_bin_path, to_string(PartitionSize).c_str());
            strcat(file_bin_path, "_");
            strcat(file_bin_path, to_string(pid).c_str());
            strcat(file_bin_path, ".bin");

            int * write_buff = (int *)malloc(SZ_VERTEX*(VERTEX_PER_BLK+2));
            FILE * bfile = fopen(file_bin_path, "wb");
            if(!bfile){
                cout<<"Failed to open file "<<file_bin_path<<endl;
                exit(1);
            }
            int x;
            int read_number=0;
            int degree=0;

            if(scResult.size() != borderVector.size()){
                cout<<"Inconsistent border size! "<<scResult.size()<<" "<<borderVector.size()<<endl; exit(1);
            }

            for(int bid=0;bid<borderVector.size();++bid){
                int ID=borderVector[bid];

//            x=fscanf(textFile, "%d", &node);
                x=ID;
                write_buff[read_number++]=x;
                if(read_number == VERTEX_PER_BLK){
                    fwrite(write_buff,SZ_VERTEX,read_number,bfile);
                    read_number = 0;
                }

//            if(scResult.find(ID)==scResult.end()){//if not found
//                cout<<"Not found in scResult! "<<ID<<endl;
//                exit(1);
//            }

//            x=fscanf(textFile, "%d", &degree);
                degree = scResult[bid].size();
                write_buff[read_number++]=degree;
                if(read_number == VERTEX_PER_BLK){
                    fwrite(write_buff,SZ_VERTEX,read_number,bfile);
                    read_number = 0;
                }

                for(int j=0; j<degree; ++j){

//                x=fscanf(textFile, "%d", &node2);
                    x=scResult[bid][j].ID2;
                    write_buff[read_number++]=x;
                    if(read_number == VERTEX_PER_BLK){
                        fwrite(write_buff,SZ_VERTEX,read_number,bfile);
                        read_number = 0;
                    }

                    for(int i=0; i<NUM_OF_CRITERIA; ++i){
//                    x=fscanf(textFile, "%d", &weight);
                        x= scResult[bid][j].getW(i);
                        if(x<=0){
                            cout<<" Wrong shortcut weight between "<<ID<<" and "<<scResult[bid][j].ID2<<" : "<<x<<"; criterion "<<i<<", partition "<<pid<<endl;
                            exit(1);
                        }

                        write_buff[read_number++]=x;
                        if(read_number == VERTEX_PER_BLK){
                            fwrite(write_buff,SZ_VERTEX,read_number,bfile);
                            read_number = 0;
                        }
                    }
                }
            }

            if(read_number>0){
                fwrite(write_buff,SZ_VERTEX,read_number,bfile);
                read_number = 0;
            }

            fclose(bfile);
            free(write_buff);
        }
        else{
            string filename=string(DataPath)+dataset+"/"+aggregateStrategy+"/Partitions/"+dataset+"_Shortcut_"+partMethod+"_"+to_string(PartitionSize)+"_"+to_string(pid);

            ofstream outFile(filename,ios::out);
            if(!outFile.is_open()){
                cout<<"Cannot open file "<<filename<<endl; exit(1);
            }

            if(scResult.size() != borderVector.size()){
                cout<<"Inconsistent border size! "<<scResult.size()<<" "<<borderVector.size()<<endl; exit(1);
            }

            int ID,degree,weight;
            for(int bid=0;bid<borderVector.size();++bid){
                ID=borderVector[bid];
                degree = scResult[bid].size();
                outFile<<ID<<" "<<degree;
                for(int j=0; j<degree; ++j){

//                x=fscanf(textFile, "%d", &node2);
                    outFile<<" "<<scResult[bid][j].ID2;

                    for(int i=0; i<NUM_OF_CRITERIA; ++i){
//                    x=fscanf(textFile, "%d", &weight);
                        weight= scResult[bid][j].getW(i);
                        outFile<<" "<<weight;
                    }
                }
                outFile<<endl;
            }

            outFile.close();

        }


    }
    // shortcut construction for partition pid
    /*void EMMCGraph::ShortcutConstruct(int pid, vector<vector<MCEdgeT>>& scResult){
//        cout<<"Building the shortcuts of partition "<<pid<<endl;
//        ReadPartition(pid,NeighborsMC[pid]);
//            unordered_set<int> set_A; set_A.clear();
//            for(int i=partitions_info[pid];i<partitions_info[pid+1];++i){
//                set_A.insert(i);
//            }
//            DFS_CC(pid, NeighborsMC[pid], set_A, node_num);
//        unordered_map<NodeId,vector<MCEdgeT>> scResult;
//        scResult.clear();
//        cout<<"Flag 1"<<endl;
        int fullNum=0;
        int temp_id;
        unordered_set<int> testBorder; testBorder.clear();
//                borderSet.insert(borderVectors[pid].begin()+i,borderVectors[pid].end());

        if(borderVectors[pid].empty()){
            cout<<"Wrong! Boundary vertex number of partition "<<pid<<" is "<<borderVectors[pid].size();
            exit(1);
        }

        int temp_i=0;

//        scResult.assign(borderVectors[pid].size(),vector<MCEdgeT>());
        Timer tt;
        for(int bid=0;bid<borderVectors[pid].size();++bid){
//            if(temp_i%100==0){
//                cout<<"temp "<<temp_i<<" "<<pid<<endl;
//            }
//            ++temp_i;
            tt.start();
            int ID=borderVectors[pid][bid];
            bool flag_full=true;
            if(BoundaryTag[temp_id]==2){
                flag_full=true;
            }else if(BoundaryTag[temp_id]==1){
                flag_full=false;
            }else{
                cout<<"Wrong BoundaryTag! "<<int(BoundaryTag[temp_id])<<endl; exit(1);
            }
//            for (int i = get<1>(EMEdgesIndex[ID]); i < get<2>(EMEdgesIndex[ID]); ++i) {
//                temp_id = NeighborsMC[pid][i].ID2;
//                if(BoundaryTag[temp_id]==0){//if temp_id is not boundary vertex
//                    flag_full = false;
//                }
//            }

            vector<unordered_map<NodeId,Distance>> shortcuts;
            shortcuts.assign(num_of_cri,unordered_map<NodeId,Distance>());
            set<NodeId> newNeighbor;

            if(flag_full){
//                if(BoundaryTag[ID]!=2){
//                    cout<<"Wrong Full BoundaryTag! "<<BoundaryTag[ID]<<endl; exit(1);
//                }
//                cout<<"Full-connected boundary vertex! "<<ID<<endl;
                ++fullNum;
                for (int i = get<1>(EMEdgesIndex[ID]); i < get<2>(EMEdgesIndex[ID]); ++i) {
                    temp_id = NeighborsMC[pid][i].ID2;
                    newNeighbor.insert(temp_id);
                    for(int cri=0;cri<num_of_cri;++cri){
                        shortcuts[cri].insert({temp_id,NeighborsMC[pid][i].getW(cri)});
                    }
                    if(node_to_cluster[temp_id]==pid){
                        testBorder.insert(temp_id);
                    }
                }

            }
            else{//if it is half-connected boundary vertex
                if(BoundaryTag[ID]!=1){
                    cout<<"Wrong Half BoundaryTag! "<<BoundaryTag[ID]<<endl; exit(1);
                }
                testBorder.insert(borderVectors[pid].begin(),borderVectors[pid].end());
            }

            tt.stop();
//            cout<<"T1: "<<tt.GetRuntime()<<endl;
            tt.start();
            if(!testBorder.empty()){//not empty
                for(int cri=0;cri<num_of_cri;++cri){
                    DijkstraSC(bid,cri,borderVectors[pid],testBorder,shortcuts[cri], newNeighbor);
                }
            }
            else{//if empty
                if(!flag_full){//if half-connected vertex
                    cout<<"Wrong for half-connected vertex! BorderSet is empty! "<<flag_full<<" "<<ID<<endl; exit(1);
                }
            }
            tt.stop();
//            cout<<"T2: "<<tt.GetRuntime()<<endl;
            tt.start();
            if(newNeighbor.empty()){
                cout<<"Empty neighbor set! "<<ID<<endl; exit(1);
            }
//                cout<<"Neighbor number: "<<newNeighbor.size()<<endl;
            /// aggregate the shortcuts
            MCEdgeT tempEdge;
            vector<MCEdgeT> tempEdges;
            for(auto it=newNeighbor.begin();it!=newNeighbor.end();++it){
                tempEdge.ID1=ID, tempEdge.ID2=*it;
                for(int cri=0;cri<num_of_cri;++cri){
                    if(shortcuts[cri].find(tempEdge.ID2)!=shortcuts[cri].end()){//if found
                        tempEdge.putW(cri,shortcuts[cri][tempEdge.ID2]);
                    }else{
                        tempEdge.putW(cri,INF100);
                    }
                }
//                    ShortcutsMC[pid].emplace_back(tempEdge);
                tempEdges.emplace_back(tempEdge);
            }
            scResult[bid]=tempEdges;
            tt.stop();
//            cout<<"T3: "<<tt.GetRuntime()<<endl;
        }
//        cout<<"Full-connected boundary vertex number of partition " << pid <<": "<<fullNum<<" "<<borderVectors[pid].size()<<endl;
        /// write to disk
//        WritePartitionShortcut(pid,borderVectors[pid],scResult);
//        NeighborsMC[pid].clear();
    }*/
    // shortcut construction for partition pid, parallel
    void EMMCGraph::ShortcutConstruct(int pid, vector<vector<MCEdgeT>>& parResult){
//        cout<<"Building the shortcuts of partition "<<pid<<endl;
//        ReadPartition(pid,NeighborsMC[pid]);
//            unordered_set<int> set_A; set_A.clear();
//            for(int i=partitions_info[pid];i<partitions_info[pid+1];++i){
//                set_A.insert(i);
//            }
//            DFS_CC(pid, NeighborsMC[pid], set_A, node_num);
//        unordered_map<NodeId,vector<MCEdgeT>> scResult;
//        scResult.clear();
//        cout<<"Flag 1"<<endl;
        int fullNum=0;
        int temp_id;
//        unordered_set<int> testBorder; testBorder.clear();
//                borderSet.insert(borderVectors[pid].begin()+i,borderVectors[pid].end());

        if(borderVectors[pid].empty()){
            cout<<"Wrong! Boundary vertex number of partition "<<pid<<" is "<<borderVectors[pid].size();
            exit(1);
        }

//        scResult.assign(borderVectors[pid].size(),vector<MCEdgeT>());
        Timer tt;


        vector<int> vertices;
        for(int bid=0;bid<borderVectors[pid].size();++bid){
            ShortcutVertex(bid,pid,parResult);
//        vertices.push_back(bid);
        }

//        cout<<"Full-connected boundary vertex number of partition " << pid <<": "<<fullNum<<" "<<borderVectors[pid].size()<<endl;
        /// write to disk
//        WritePartitionShortcut(pid,borderVectors[pid],parResult);
//        NeighborsMC[pid].clear();
    }

    void EMMCGraph::ShortcutConstructV(vector<NodeId>& p, vector<vector<vector<MCEdgeT>>>& scResult, vector<double>& timeRecord){
        Timer tt;
        int ti=0;
        for(int i=0;i<p.size();++i){
            int pid=p[i];

            tt.start();
            ShortcutConstruct(pid, scResult[pid]);
            tt.stop();
            timeRecord[pid]=tt.GetRuntime();
            if(ti%100==0){
                cout<<"Partition "<<i<<" "<<pid<<" "<<borderVectors[pid].size()<< ": "<<tt.GetRuntime()<<" s."<< endl;
            }
            ti++;
        }
    }
    //shortcut construction for BSHP, batch processing
    void EMMCGraph::ShortcutConstruction(){
        cout<<"Begin shortcut construction..."<<endl;
        GetBoundary();

        string filename=string(DataPath)+dataset+"/"+aggregateStrategy+"/Partitions/"+dataset+"_Shortcut_"+partMethod+"_"+to_string(PartitionSize)+"_0.bin";
        ifstream istrm(filename, ios::binary);
        if (!istrm.is_open()){
//        if(true){
            istrm.close();
            std::cout << "failed to open " << filename << '\n';
            cout<<"Computing shortcuts..."<<endl;
            NeighborsMC.assign(partition_number,vector<MCEdgeT>());
            EMEdgesIndex.assign(node_num,make_tuple(false,-1,-1));
            Timer tt;
            tt.start();

            vector<double> timeRecord(partition_number,0);

            //allocate space
            vector<vector<vector<MCEdgeT>>> scResult;
            scResult.assign(partition_number,vector<vector<MCEdgeT>>());

            vector<vector<NodeId>> processID;
            processID.assign(threadnum, vector<NodeId>());
            vector<NodeId> vertices;
//        MCEdgesT.assign(node_num,vector<MCEdgeT>());
//            ReadMCEdges(string(DataPath)+dataset+"/"+dataset+".MCEdges");
            for(int pid=0;pid<partition_number;++pid){
                vertices.emplace_back(pid);
                ReadPartition(pid,NeighborsMC[pid]);
//            ReadPartitionMCEdges(pid);
                scResult[pid].assign(borderVectors[pid].size(), vector<MCEdgeT>());
            }
            /// multi thread
            ThreadDistribute(vertices, processID);
            cout<<"Batch number: "<<processID[0].size()<<endl;
            boost::thread_group thread;
            for(auto j=0;j<processID.size();++j){
                thread.add_thread(new boost::thread(&EMMCGraph::ShortcutConstructV, this, boost::ref(processID[j]) , boost::ref(scResult), boost::ref(timeRecord)));
//                cout<<"ProcessID "<<j<<endl;
//            ShortcutConstructV(processID[j], scResult, timeRecord);
            }
            thread.join_all();

        //// single thread
//        for(auto j=0;j<processID.size();++j){
//            ShortcutConstructV(processID[j], scResult, timeRecord);
//        }

//        vector<int> pp;
//        for(int i=0;i<5;++i){
//            int pid=rand()%partition_number;
//            pp.push_back(pid);
//        }
//        ShortcutConstructV(pp, scResult, timeRecord);
            tt.stop();

            double t_max=0;
            double t_min=INF;
            double t_ave=0;
            int temp_num=0;
            for(int pid=0;pid<partition_number;++pid){
                t_max = max(t_max,timeRecord[pid]);
                if(timeRecord[pid]!=0){
                    t_min = min(t_min,timeRecord[pid]);
                    temp_num++;
                }
                t_ave += timeRecord[pid];
            }
            t_ave/=temp_num;
            cout<<"Single partition's shortcut construction. Average time: "<<t_ave<<" s; Maximum time: "<<t_max<<" s; Minimum time: "<<t_min<<" s."<<endl;

//            cout<<"Flag 1"<<endl;
            NeighborsMC.clear();
            EMEdgesIndex.clear();
            for(int pid=0;pid<partition_number;++pid){
                WritePartitionShortcut(pid,borderVectors[pid],scResult[pid], true);
                if(pid%100==0){
                    WritePartitionShortcut(pid,borderVectors[pid],scResult[pid], false);
                }
            }
            scResult.clear();
//            cout<<"Flag 2"<<endl;

            cout<<"Time for shortcut construction: "<<tt.GetRuntime()<<" s."<<endl;
        }
        else{
            istrm.close();
        }

//        exit(0);
    }
    //compute the boundnary vertex
    void EMMCGraph::GetBoundary() {
//        vector<vector<int>> boundaryVertex;
        ClusterInfoLoad(true);//load partition information, using node-to-cluster mapping
        cout << "Node number: " << node_num << "\t\tEdge number: " << edge_num << endl;

        string w_file=string(DataPath)+dataset+"/" + aggregateStrategy+"/"+dataset+"_"+partMethod+"_"+to_string(PartitionSize)+".boundary";
        ifstream inFile(w_file, ios::in);
        if (!inFile) {
//        if(true){
//            inFile.close();
            cout << "File opening failed. " <<w_file<< endl;
            cout<<"Compute boundary vertices..."<<endl;

            borderVectors.assign(partition_number,vector<int>());
            uint borderNum=0;
            unordered_set<int> borderSet; borderSet.clear();
            uint maxBorder=0;
            BoundaryTag.assign(node_num, 0);

            for(int pid=0;pid<partition_number;++pid){
                char filePath[300];
                int u = 0, u_deg;
                int v, weight;
                int ID1, ID2;

                strcpy(filePath, DataPath);
                strcat(filePath, dataset.c_str());
                strcat(filePath, "/");
                strcat(filePath, aggregateStrategy.c_str());
                strcat(filePath, "/Partitions/");
                strcat(filePath, dataset.c_str());
                strcat(filePath, "_Partition_");
                strcat(filePath, partMethod.c_str());
                strcat(filePath, "_");
                strcat(filePath, to_string(PartitionSize).c_str());
                strcat(filePath, "_");
                strcat(filePath, to_string(pid).c_str());
                strcat(filePath, ".bin");

                FILE* file = fopen(filePath, "rb");
                if(!file){
                    cout<<"Failed to open file "<<filePath<<endl;
                    exit(1);
                }
                assert(file != NULL);
                ReadBuffer read_io(file);

                unordered_map<NodeId,vector<NodeId>> bNeighbors; //record in-partition neighbors
                bNeighbors.clear();
                vector<NodeId> Neigh;

                int temp,clusterSize;
                read_io.read(&temp);//cluster id
                read_io.read(&clusterSize);//cluster size

                while (!read_io.isend) {//if not read out
                    read_io.read(&u);//get node id
                    ID1 = u;

                    Neigh.clear();

                    read_io.read(&u_deg);//get node degree
                    bool flag_b= false;

                    for (int i = 0; i < u_deg; ++i) {
                        read_io.read(&v);//read end node
                        ID2 = v;
                        if(node_to_cluster[ID2] == pid){
                            Neigh.emplace_back(ID2);
                        }

                        if(node_to_cluster[ID2] != pid && !flag_b){
                            borderVectors[pid].emplace_back(ID1);
                            borderSet.insert(ID1);
                            BoundaryTag[ID1] = 1;
                            flag_b = true;
                        }

                        for (int j = 0; j < NUM_OF_CRITERIA; ++j) {
                            read_io.read(&weight);// read edge weight
                        }

                    }

                    bNeighbors.insert({ID1,Neigh});

                }
                if(borderVectors[pid].empty()){
                    cout<<"Wrong! "<<pid<<" "<<borderVectors[pid].size()<<" "<<clusterSize<<endl;
                    exit(1);
                }

                bool flag_full=true;
                int fullNum=0;
                for(int bid=0;bid<borderVectors[pid].size();++bid){
                    ID1=borderVectors[pid][bid];
                    flag_full = true;
                    if(bNeighbors.find(ID1)==bNeighbors.end()){
                        cout<<"Wrong! Cannot find "<<ID1<<endl; exit(1);
                    }
                    if(bNeighbors[ID1].empty()){
                        BoundaryTag[ID1] = 2;
                        fullNum++;
                    }else{
                        for(int i=0;i<bNeighbors[ID1].size();++i){
                            ID2=bNeighbors[ID1][i];
                            if(BoundaryTag[ID2] == 0){
                                flag_full = false;
                            }
                        }
                        if(flag_full){
                            BoundaryTag[ID1] = 2;
                            fullNum++;
                        }
                    }

                }

                borderNum += borderVectors[pid].size();
                if(maxBorder<borderVectors[pid].size()){
                    maxBorder=borderVectors[pid].size();
                }
                fclose(file);
            }
            cout<<"Boundary vertex number: "<<borderNum<<" "<<borderSet.size()<<". Average boundary vertex number: "<<borderNum/partition_number<<". Max boundary number for all partitions: "<<maxBorder <<endl;


            /// write boundary vertex to disk
            ofstream outFile(w_file, ios::out);
            if (!outFile) {
                cout << "File opening failed. " <<w_file<< endl;
                exit(1);
            }
            outFile<<partition_number<<endl;
            for(int pid=0;pid<partition_number;++pid){
                outFile<<borderVectors[pid].size();
                for(int i=0;i<borderVectors[pid].size();++i){
                    int ID=borderVectors[pid][i];
                    outFile<<" "<<ID<<" "<<int(BoundaryTag[ID]);
                }
                outFile<<endl;
            }
            outFile.close();
        }
        else{
            inFile.close();
            ReadBoundary(w_file);
        }

    }
    //read the boundary vertices
    void EMMCGraph::ReadBoundary(string filename){
        ifstream inFile(filename, ios::in);
        if (!inFile) {
            cout << "File opening failed. " <<filename<< endl;
            exit(1);
        }
        cout<<"Reading boundary information."<<endl;
        string line;
        getline(inFile,line);
        vector<string> vs;
        boost::split(vs,line,boost::is_any_of(" \t"),boost::token_compress_on);
        int partNum=stoi(vs[0]);
        assert(partition_number==partNum);
        borderVectors.assign(partNum,vector<int>());
        BoundaryTag.assign(node_num, 0);
        int borderNum=0;
        int ID;
        int maxBorder=0;
        int pborder=0;
        for(int pid=0;pid<partNum;++pid){
            pborder=0;
            vs.clear();
            getline(inFile,line);
            boost::split(vs,line,boost::is_any_of(" \t"),boost::token_compress_on);
            int bNum=stoi(vs[0]);
            assert(vs.size()==2*bNum+1);
            for(int i=0;i<2*bNum;i+=2){
                ID=stoi(vs[i+1]);
                borderVectors[pid].emplace_back(ID);
                BoundaryTag[ID] = stoi(vs[i+2]);
                ++borderNum;
                ++pborder;
            }
            maxBorder = max(maxBorder, pborder);
        }
        inFile.close();
        cout<<"Boundary vertex number: "<<borderNum<<" . Average boundary vertex number: "<<borderNum/partition_number<<". Max boundary number for all partitions: "<<maxBorder <<endl;
    }

    //function of vertex allocation
    void EMMCGraph::ThreadDistribute(vector<NodeId>& vertices, vector<vector<NodeId>>& processID){
//        processID.assign(threadnum, vector<NodeId>());
        int pid=0;
        for(int i=0;i<vertices.size();++i){
            pid=i%processID.size();
            processID[pid].emplace_back(vertices[i]);
        }
    }

    pair<int, unsigned long long> EMMCGraph::DFS_CC(int pid, vector<MCEdgeT> & Edges, unordered_set<int> & set_A, int nodenum) {//set_A: the vertex set
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
//                for (auto it = Edges[item_id].begin(); it != Edges[item_id].end(); ++it) {
                for (int i = get<1>(EMEdgesIndex[item_id]); i < get<2>(EMEdgesIndex[item_id]); ++i) {
                    temp_id = Edges[i].ID2;
                    if(node_to_cluster[temp_id] != pid){
                        continue;
                    }

//                    temp_id = it->first;
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

        std::sort(CCs.begin(), CCs.end());

        return make_pair(LCC.first.size(),LCC.second);//Node number and edge number of LCC
    }


    //SMHP
    void EMMCGraph::EM_MC_ShortcutSearch(int node_start, int node_end){
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
        MCHotPool<VectorMCEdgesEMTuple_IO> mc_HotPool(Partition_N_SC);       //Hot pool for one-pass algorithms
        mc_HotPool.init(partition_number);
        mc_HotPool.set_capacity(Partition_N_SC);
//        cout<<"Partition capacity: "<<Partition_N_SC<<endl;
        Timer tt;
        tt.start();

        //Initiation of start node
        for (int i = 0; i < num_of_cri; ++i) {
            set_cri.insert(i);
            mc_cost[i][node_start] = 0;//cost of start node
//            EM_MC_PQueue[i]->push(VertexCost(node_start,0));
//            cri_to_cluster[i] = node_to_cluster[node_start];
            cri_pqueue.update(i, 0);
        }

        //Search the source partition
        int pid_start=node_to_cluster[node_start];
        int pid_end=node_to_cluster[node_end];
        cout<<"pid_start: "<<pid_start<<", pid_end: "<<pid_end<<endl;
//        GraphMCReadCluster_New(mc_HotPool, node_start, pid_start, MultiHops);
        vector<pair<NodeId,vector<Distance>>> localDisV;
        DijkstraPartitionMC(node_start,node_end, localDisV);
        ++p_num;
        for(auto it=localDisV.begin();it!=localDisV.end();++it){
            int bid=it->first;

            for(cri_i=0; cri_i<num_of_cri; ++cri_i){
                if(it->second[cri_i]<INF){
                    mc_cost[cri_i][bid] = it->second[cri_i];
                    EM_MC_PQueue[cri_i]->push(VertexCost(bid, mc_cost[cri_i][bid]));
                }
            }
        }
        localDisV.clear();

        for(cri_i=0;cri_i<num_of_cri;++cri_i){
            cri_to_cluster[cri_i] = node_to_cluster[EM_MC_PQueue[cri_i]->top().id];
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
//            GraphMCReadCluster_New(mc_HotPool, item_id, cri_to_cluster[cri_pqueue.top_value()], MultiHops);//if in disk, read corresponding partition
            partition_id = cri_to_cluster[cri_pqueue.top_value()];
            if(partition_id != pid_end){
                ReadPartitionShortcut(partition_id, mc_HotPool);
            }else{//if target partition
                GraphMCReadCluster_New(mc_HotPool, item_id, partition_id, MultiHops);//if in disk, read corresponding partition
            }

            ++p_num;
            while (!set_remain.empty()) {
                cri_i = *set_remain.begin();//deal with each criterion

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

//                cout<<cri_i<<" "<<item_id<<" "<<partition_id<<endl;

                while (get<0>(mc_HotPool.clusterStatus[partition_id]) == 'I' || get<0>(mc_HotPool.clusterStatus[partition_id]) == 'S') {//if found
                    if (item_id == node_end) {//if reach target node, this criterion ends.
                        //                mc_pqueue[cri_i].clear();
                        mc_finished[cri_i] = true;
                        set_cri.erase(cri_i);
                        mc_min_cost[cri_i] = mc_cost[cri_i][item_id];
                        cout<<"Minimal cost of criterion " << cri_i << ": "<<mc_min_cost[cri_i]<<endl;
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
                                    for(auto it=cluster_to_node[partition_id].begin();it!=cluster_to_node[partition_id].end();++it){
                                        get<1>(EMEdgesIndex[*it]) = -1;//hard delete
                                        get<2>(EMEdgesIndex[*it]) = -1;
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
    //function of conducting the Dijkstra's search on the partition of node_start
    void EMMCGraph::DijkstraPartitionMC(int node_start, int node_end, vector<pair<NodeId,vector<Distance>>>& localDisV){
        int partition_id=node_to_cluster[node_start];
        cout<<partition_id<<endl;
        vector<MCEdgeT> NeighborsP;
        ReadPartitionBi(partition_id,NeighborsP);
//        unordered_set<int> set_A; set_A.clear();
//        for(int i=partitions_info[partition_id];i<partitions_info[partition_id+1];++i){
//            set_A.insert(i);
//        }
//        DFS_CC(partition_id, NeighborsP, set_A, node_num);
        unordered_map<NodeId,vector<Distance>> localDis;
        localDis.clear();
        for(auto it=cluster_to_node[partition_id].begin();it!=cluster_to_node[partition_id].end();++it){
            localDis.insert({*it,vector<Distance>(num_of_cri,INF)});
        }
        int item_id, temp_id;

        /// multi-thread
        boost::thread_group thread;
        for(auto cri_i=0;cri_i<num_of_cri;++cri_i){
            thread.add_thread(new boost::thread(&EMMCGraph::DijkstraPartitionBi, this, node_start, cri_i, boost::ref(NeighborsP), boost::ref(localDis)));
        }
        thread.join_all();
        // single thread
//        for(auto cri_i=0;cri_i<num_of_cri;++cri_i){
//            DijkstraPartition(node_start, cri_i, NeighborsP, localDis);
//        }

        NeighborsP.clear();
//        EMEdgesIndex_Bi.assign(node_num, make_tuple(false,false,-1,-1));
        for(auto it=cluster_to_node[partition_id].begin();it!=cluster_to_node[partition_id].end();++it){
//            cout<<"redress "<<*it<<endl;
            get<2>(EMEdgesIndex_Bi[*it])=-1; get<3>(EMEdgesIndex_Bi[*it])=-1;
        }

        vector<NodeId> candidates;
        for(auto it=borderVectors[partition_id].begin();it!=borderVectors[partition_id].end();++it){
            item_id=*it;
            candidates.emplace_back(item_id);//the boundary vertices of partition_id are popped out after DijkstraPatition
            if(localDis.find(*it)!=localDis.end()){//if found
//                if(localDis[*it][0]<INF){
//                    cout<<*it;
//                    for(auto it2=localDis[*it].begin();it2!=localDis[*it].end();++it2){
//                        cout<<" "<<*it2;
//                    }
//                    cout<<endl;
//                }
                localDisV.emplace_back(*it,localDis[*it]);

            }else{
                cout<<"Not found this candidate. "<<*it<<endl;
                exit(1);
            }
        }

        if(localDis.find(node_end)!=localDis.end()){//if found
            localDisV.emplace_back(node_end,localDis[node_end]);
        }
    }
    //function of conducting the Dijkstra's search on the partition of node_start, with criterion cri_i
    void EMMCGraph::DijkstraPartition(int node_start, int cri_i, vector<MCEdgeT>& NeighborsP, unordered_map<NodeId,vector<Distance>>& localDis){
        int partition_id=node_to_cluster[node_start];
//        gbxxl::PriorityQueue pqueue(gbxxl::PQ_Pool);
        priority_queue<pair<int,int>, vector<pair<int,int>>, IntPairCompareLess> pqueue;

        //variables
        unordered_map<NodeId,Distance> cost; cost.clear();
        unordered_map<NodeId,bool> closed; closed.clear();
        for(auto it=cluster_to_node[partition_id].begin();it!=cluster_to_node[partition_id].end();++it){
            int id=*it;
            cost.insert({id,INF});
            closed.insert({id,false});
        }

        NodeId item_id;
        bool flag_finish = false;
        NodeId temp_id;
        Distance temp_dis;
        Distance item_dis;
        int temp_number;
        int temp_degree = 0;

        Timer tt;
        tt.start();

        //Initiation of start node
        cost[node_start] = 0;//cost of start node
        pqueue.push(make_pair(node_start,0));

        //Iteration
        while (!flag_finish && !pqueue.empty()) {//for every node in pqueue !mc_finished[0] || !mc_finished[1] || !mc_finished[2]
            while (closed[pqueue.top().first]) {//if found in visited set  && !EM_MC_PQueue[cri_i]->empty()
                pqueue.pop();
                if(pqueue.empty()){
                    flag_finish = true;
                    goto jump1;
                }
            }
            item_id = pqueue.top().first;// top min item
            item_dis = pqueue.top().second;// top min item

            if(node_to_cluster[item_id]!=partition_id){
                cout<<"Wrong item_id of DijkstraPartition! "<< item_id<<" "<<node_to_cluster[item_id]<<endl; exit(1);
            }

            //set closed
            closed[item_id] = true;
            if(localDis.find(item_id)!=localDis.end()){//if found
                localDis[item_id][cri_i] = cost[item_id];//record the local distance
            }else{
                cout<<"Wrong item_id: "<<item_id<<" "<<node_to_cluster[item_id]<<endl; exit(1);
            }


            pqueue.pop();//pop min item
            ++visited_number;//update #visited_number
            temp_degree = 0;
            if(node_to_cluster[item_id] != partition_id){
//                cout<<item_id<<" "<<node_to_cluster[item_id]<<" "<<localDis[item_id][cri_i]<<endl;
                continue;
            }
            for (int i = get<1>(EMEdgesIndex[item_id]); i < get<2>(EMEdgesIndex[item_id]); ++i) {
                ++temp_degree;
                temp_id = NeighborsP[i].ID2;
                Distance temp_w = NeighborsP[i].getW(cri_i);
                if(node_to_cluster[temp_id] != partition_id){
                    continue;
//                    if(cost.find(temp_id) == cost.end()){//if not found
////                        cout<<"Insert new temp_id "<<temp_id<<" "<<node_to_cluster[temp_id]<<endl;
//                        cost.insert({temp_id,INF});
//                        closed.insert({temp_id,false});
//                    }
                }
                if (temp_w==INF100)
                    continue;
                if (closed[temp_id])//if closed
                    continue;
                temp_dis = item_dis + temp_w;
                if (cost[temp_id] > temp_dis) {//slack operation
                    cost[temp_id] = temp_dis;
                    pqueue.push(make_pair(temp_id, temp_dis));
                }
            }


        }
        jump1:
//        cout<<"Boundary shortcut number: "<<stainN<<" "<<borderVector.size()<<", criterion: "<< cri_i <<endl;
        tt.stop();
    }
    //function of conducting the Dijkstra's search on the partition of node_start, with criterion cri_i
    void EMMCGraph::DijkstraPartitionBi(int node_start, int cri_i, vector<MCEdgeT>& NeighborsP, unordered_map<NodeId,vector<Distance>>& localDis){
        int partition_id=node_to_cluster[node_start];
//        gbxxl::PriorityQueue pqueue(gbxxl::PQ_Pool);
        priority_queue<pair<int,int>, vector<pair<int,int>>, IntPairCompareLess> pqueue;

        //variables
        unordered_map<NodeId,Distance> cost; cost.clear();
        unordered_map<NodeId,bool> closed; closed.clear();
        for(auto it=cluster_to_node[partition_id].begin();it!=cluster_to_node[partition_id].end();++it){
            int id=*it;
            cost.insert({id,INF});
            closed.insert({id,false});
        }

        NodeId item_id;
        bool flag_finish = false;
        NodeId temp_id;
        int temp_w;
        Distance temp_dis;
        Distance item_dis;
        int temp_number;
        int temp_degree = 0;

        Timer tt;
        tt.start();

        //Initiation of start node
        cost[node_start] = 0;//cost of start node
        pqueue.push(make_pair(node_start,0));

        //Iteration
        while (!flag_finish && !pqueue.empty()) {//for every node in pqueue !mc_finished[0] || !mc_finished[1] || !mc_finished[2]
            while (closed[pqueue.top().first]) {//if found in visited set  && !EM_MC_PQueue[cri_i]->empty()
                pqueue.pop();
                if(pqueue.empty()){
                    flag_finish = true;
                    goto jump1;
                }
            }
            item_id = pqueue.top().first;// top min item
            item_dis = pqueue.top().second;// top min item

            if(node_to_cluster[item_id]!=partition_id){
                cout<<"Wrong item_id of DijkstraPartition! "<< item_id<<" "<<node_to_cluster[item_id]<<endl; exit(1);
            }

            //set closed
            closed[item_id] = true;
            if(localDis.find(item_id)!=localDis.end()){//if found
                localDis[item_id][cri_i] = cost[item_id];//record the local distance
            }else{
                cout<<"Wrong item_id: "<<item_id<<" "<<node_to_cluster[item_id]<<endl; exit(1);
            }


            pqueue.pop();//pop min item
            ++visited_number;//update #visited_number
            temp_degree = 0;
            if(node_to_cluster[item_id] != partition_id){
//                cout<<item_id<<" "<<node_to_cluster[item_id]<<" "<<localDis[item_id][cri_i]<<endl;
                continue;
            }
            for (int i = get<2>(EMEdgesIndex_Bi[item_id]); i < get<3>(EMEdgesIndex_Bi[item_id]); ++i) {
                ++temp_degree;
                temp_id = NeighborsP[i].ID2;
                temp_w = NeighborsP[i].getW(cri_i);
                if(node_to_cluster[temp_id] != partition_id){
                    continue;
//                    if(cost.find(temp_id) == cost.end()){//if not found
////                        cout<<"Insert new temp_id "<<temp_id<<" "<<node_to_cluster[temp_id]<<endl;
//                        cost.insert({temp_id,INF});
//                        closed.insert({temp_id,false});
//                    }
                }
                if (temp_w==INF100)
                    continue;
                if (closed[temp_id])//if closed
                    continue;
                temp_dis = item_dis + temp_w;
                if (cost[temp_id] > temp_dis) {//slack operation
                    cost[temp_id] = temp_dis;
                    pqueue.push(make_pair(temp_id, temp_dis));
                }
            }


        }
        jump1:
//        cout<<"Boundary shortcut number: "<<stainN<<" "<<borderVector.size()<<", criterion: "<< cri_i <<endl;
        tt.stop();
    }
    //function of reading the shortcuts of partition
    void EMMCGraph::ReadPartitionShortcut(int target_p_id, MCHotPool<VectorMCEdgesEMTuple_IO> & mcHotPool){
        int evict_p_id = mcHotPool.LRU_IMCluster.put(target_p_id,1);
        //// evict old partition for new one
        if(evict_p_id != -1){//if it is necessary to evict old partition
            assert(get<0>(mcHotPool.clusterStatus[evict_p_id]) == 'I');
            unsigned int evict_p_size = get<2>(mcHotPool.clusterStatus[evict_p_id]);
            uint temp_sz = 0;
            temp_sz = get<1>(mcHotPool.clusterStatus[evict_p_id])*alpha_multi;
//            if(algo_choice == OneHop){
//                temp_sz = get<1>(mcHotPool.clusterStatus[evict_p_id])*alpha;
//            }else if(algo_choice == MultiHops){
//                temp_sz = get<1>(mcHotPool.clusterStatus[evict_p_id])*alpha_multi;
//            }
            assert(temp_sz>0);
            assert(evict_p_size>0);

            if(evict_p_size > temp_sz){//if the remaining size is too large or no node remain, evict all
                get<2>(mcHotPool.clusterStatus[evict_p_id]) = 0;
                get<0>(mcHotPool.clusterStatus[evict_p_id]) = 'D';
                for(auto it=cluster_to_node[evict_p_id].begin();it!=cluster_to_node[evict_p_id].end();++it){
                    get<1>(EMEdgesIndex[*it]) = -1;//hard delete
                    get<2>(EMEdgesIndex[*it]) = -1;
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
        char filePath[300];
        int u, u_deg;
        int v, weight;
        MCEdgeT temp_edge;

        strcpy(filePath, DataPath);
        strcat(filePath, dataset.c_str());
        strcat(filePath, "/");
        strcat(filePath, aggregateStrategy.c_str());
        strcat(filePath, "/Partitions/");
        strcat(filePath, dataset.c_str());
        strcat(filePath, "_Shortcut_");
        strcat(filePath, partMethod.c_str());
        strcat(filePath, "_");
        strcat(filePath, to_string(PartitionSize).c_str());
        strcat(filePath, "_");
        strcat(filePath, to_string(target_p_id).c_str());
        strcat(filePath, ".bin");
        Timer tt;
        tt.start();
        FILE* file = fopen(filePath, "rb");
        if(!file){
            cout<<"Failed to open file "<<filePath<<endl;
            exit(1);
        }
        assert(file != NULL);

//        cout<<"Reading shortcut file "<<filePath<<endl;

        ReadBuffer read_io(file);
        bool flag_sizeUpdate = false;
        int temp_id = -1;
//        cout<<"Read flag 1"<<endl;
        if(get<1>(mcHotPool.clusterStatus[target_p_id]) == 0){//if this is the first read
            flag_sizeUpdate = true;
//            cout<<"First read."<<endl;
            //read partition
            while(!read_io.isend) {//if not read out
                read_io.read(&u);//get node id
//                cout<<u<<endl;
                temp_edge.ID1 = u;
                if(u != temp_id){
                    get<1>(EMEdgesIndex[u]) = mcHotPool.HotPools[target_p_id].size();
                    if(temp_id != -1){
                        get<2>(EMEdgesIndex[temp_id]) = mcHotPool.HotPools[target_p_id].size();
                    }
                    temp_id = u;
                }

//                if(vertex_cri[u-1] == -1){//if not found
                assert(vertex_cri[u] == -1);
                vertex_cri[u] = num_of_cri;
//                }
                read_io.read(&u_deg);//get node degree
                //vertex store in indexNode, edge store in e
                for (int i = 0; i < u_deg; ++i) {
                    read_io.read(&v);//read end node
                    temp_edge.ID2 = v;
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
        }
        else{//if this is not the first read
//            cout<<"Second read."<<endl;
            //read partition
            while(!read_io.isend) {//if not read out
                read_io.read(&u);//get node id
                temp_edge.ID1 = u;
                if(u != temp_id){
                    get<1>(EMEdgesIndex[u]) = mcHotPool.HotPools[target_p_id].size();
                    if(temp_id != 0){
                        get<2>(EMEdgesIndex[temp_id]) = mcHotPool.HotPools[target_p_id].size();
                    }
                    temp_id = u;
                }

                if(vertex_cri[u] == -1){//if not found
                    vertex_cri[u] = num_of_cri;
                }
                read_io.read(&u_deg);//get node degree
                //vertex store in indexNode, edge store in e
                for (int i = 0; i < u_deg; ++i) {
                    read_io.read(&v);//read end node
                    temp_edge.ID2 = v;
                    for (int j = 0; j < NUM_OF_CRITERIA; ++j) {
                        read_io.read(&weight);
                        if(j < num_of_cri){
                            temp_edge.putW(j,weight);
                        }
                    }
                    mcHotPool.HotPools[target_p_id].emplace_back(temp_edge);

                }
                if (!get<0>(EMEdgesIndex[u]))//if false, if it has not been evicted
                    get<2>(mcHotPool.clusterStatus[target_p_id]) += u_deg;
            }
            io_num += 1;//update #IO
        }


        get<2>(EMEdgesIndex[u]) = mcHotPool.HotPools[target_p_id].size();
//        cout<<mcHotPool.HotPools[target_p_id].size()<<endl;
//        mcHotPool.im_partitions.insert(target_p_id);
        //io_num += read_io.io_number;//update #IO
        fclose(file);
//        cout<<"Read done."<<endl;
        tt.stop();
        io_time += tt.GetRuntime() * 1000;
        cluster_to_IO[target_p_id] = read_io.io_number;

        get<0>(mcHotPool.clusterStatus[target_p_id]) = 'I';
    }

    //BMHPS
    /*void EMMCGraph::EM_MC_ShortcutSearch_Bi(int node_start, int node_end){
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
            mc_cost_r[i][node_end] = 0;
//            EM_MC_PQueue[i]->push(VertexCost(node_start, 0));
//            EM_MC_PQueue_r[i]->push(VertexCost(node_end, 0));
//            cri_to_cluster[i] = node_to_cluster[node_start];
//            cri_to_cluster_r[i] = node_to_cluster[node_end];
            cri_pqueue.update(i,0);
            cri_pqueue_r.update(i,0);
        }

        //Search the source partition
        int pid_start=node_to_cluster[node_start];
        int pid_end=node_to_cluster[node_end];
        cout<<"pid_start: "<<pid_start<<", pid_end: "<<pid_end<<endl;
//        GraphMCReadCluster_New(mc_HotPool, node_start, pid_start, MultiHops);
        vector<pair<NodeId,vector<Distance>>> localDisV;
        DijkstraPartitionMC(node_start, node_end, localDisV);
        ++p_num;
        for(auto it=localDisV.begin();it!=localDisV.end();++it){
            int bid=it->first;

            for(cri_i=0; cri_i<num_of_cri; ++cri_i){
                if(it->second[cri_i]<INF){
                    mc_cost[cri_i][bid] = it->second[cri_i];
                    EM_MC_PQueue[cri_i]->push(VertexCost(bid, mc_cost[cri_i][bid]));
                }
            }
        }
        localDisV.clear();
        //Search the target partition
        DijkstraPartitionMC(node_end, node_start, localDisV);
        ++p_num;
        for(auto it=localDisV.begin();it!=localDisV.end();++it){
            int bid=it->first;

            for(cri_i=0; cri_i<num_of_cri; ++cri_i){
                if(it->second[cri_i]<INF){
                    mc_cost_r[cri_i][bid] = it->second[cri_i];
                    EM_MC_PQueue_r[cri_i]->push(VertexCost(bid, mc_cost_r[cri_i][bid]));
                }
            }
        }
        localDisV.clear();

        for(cri_i=0;cri_i<num_of_cri;++cri_i){
            cri_to_cluster[cri_i] = node_to_cluster[EM_MC_PQueue[cri_i]->top().id];
            cri_to_cluster_r[cri_i] = node_to_cluster[EM_MC_PQueue_r[cri_i]->top().id];
        }

        cluster_to_cri[node_to_cluster[node_start]] = 0;
        cluster_to_cri[node_to_cluster[node_end]] = 0;
        //Iteration
        while (!EM_JudgeEmptyBi(EM_MC_PQueue, EM_MC_PQueue_r)) {//for every node in pqueue !mc_finished[0] || !mc_finished[1] || !mc_finished[2]
            set_remain = set_cri;
            /// read by synchronizing the hops
            temp_cluster_id = cri_to_cluster[cri_pqueue.top_value()];
            temp_cluster_id_r = cri_to_cluster_r[cri_pqueue_r.top_value()];
            if(cri_hops[cri_pqueue.top_value()] <= cri_hops_r[cri_pqueue_r.top_value()]){
                if(temp_cluster_id!=pid_end){
                    ReadPartitionShortcut_Bi(temp_cluster_id, mc_HotPool, cluster_to_bi, FORWARD, clusterRead);//if in disk, read corresponding partition item_id,
                }else{
                    GraphMCReadCluster_New_Bi(mc_HotPool, temp_cluster_id, cluster_to_bi, FORWARD, clusterRead);//if in disk, read corresponding partition item_id,
                }
                cluster_to_bi[temp_cluster_id] = FORWARD;
                ++p_num;
            }else{
                if(temp_cluster_id_r!=pid_start){
                    ReadPartitionShortcut_Bi(temp_cluster_id_r, mc_HotPool, cluster_to_bi, REVERSE, clusterRead);//if in disk, read corresponding partition item_id_r,
                }else{
                    GraphMCReadCluster_New_Bi(mc_HotPool, temp_cluster_id_r, cluster_to_bi, REVERSE, clusterRead);//if in disk, read corresponding partition item_id_r,
                }
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
                        cout<<"T2 Minimal cost of criterion "<<cri_i<<": "<<mc_min_cost[cri_i]<<endl;
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
                                    for(auto it=cluster_to_node[partition_id].begin();it!=cluster_to_node[partition_id].end();++it){
                                        int id=*it;
                                        get<2>(EMEdgesIndex_Bi[id]) = -1;//hard delete
                                        get<3>(EMEdgesIndex_Bi[id]) = -1;
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
                        cout<<"T2 Minimal cost of criterion "<<cri_i<<": "<<mc_min_cost[cri_i]<<endl;
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
                                    for(auto it=cluster_to_node[partition_id_r].begin();it!=cluster_to_node[partition_id_r].end();++it){
                                        get<2>(EMEdgesIndex_Bi[*it]) = -1;//hard delete
                                        get<3>(EMEdgesIndex_Bi[*it]) = -1;
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
                    cout<<"Minimal cost of criterion "<<cri_i<<": "<<mc_min_cost[cri_i]<<endl;
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
    }*/
    //BMHPS, new
    void EMMCGraph::EM_MC_ShortcutSearch_Bi(int node_start, int node_end){
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
        int cri_i=0;
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
            mc_cost_r[i][node_end] = 0;
//            EM_MC_PQueue[i]->push(VertexCost(node_start, 0));
//            EM_MC_PQueue_r[i]->push(VertexCost(node_end, 0));
//            cri_to_cluster[i] = node_to_cluster[node_start];
//            cri_to_cluster_r[i] = node_to_cluster[node_end];
            cri_pqueue.update(i,0);
            cri_pqueue_r.update(i,0);
        }

        //Search the source partition
        int pid_start=node_to_cluster[node_start];
        int pid_end=node_to_cluster[node_end];
        cout<<"pid_start: "<<pid_start<<", pid_end: "<<pid_end<<endl;
//        GraphMCReadCluster_New(mc_HotPool, node_start, pid_start, MultiHops);

        vector<pair<NodeId,vector<Distance>>> localDisV;
        if(BoundaryTag[node_start]!=0 && BoundaryTag[node_end]!=0){
            cout<<"s and t are both boundary vertices."<<endl;
            for(cri_i=0; cri_i<num_of_cri; ++cri_i) {
                EM_MC_PQueue[cri_i]->push(VertexCost(node_start, 0));
                EM_MC_PQueue_r[cri_i]->push(VertexCost(node_end, 0));
            }
        }else {
            localDisV.clear();
            if(BoundaryTag[node_start]==0){//non-boundary vertex
                cout<<"s is non-boundary vertex; ";
                DijkstraPartitionMC(node_start, node_end, localDisV);
                ++p_num;
                for(auto it=localDisV.begin();it!=localDisV.end();++it){
                    int bid=it->first;

                    for(cri_i=0; cri_i<num_of_cri; ++cri_i){
                        if(it->second[cri_i]<INF){
                            mc_cost[cri_i][bid] = it->second[cri_i];
                            EM_MC_PQueue[cri_i]->push(VertexCost(bid, mc_cost[cri_i][bid]));
                        }
                    }
                }
            }else{
                cout<<"s is boundary vertex; ";
                for(cri_i=0; cri_i<num_of_cri; ++cri_i) {
                    EM_MC_PQueue[cri_i]->push(VertexCost(node_start, 0));
                }
            }
            localDisV.clear();
            if(BoundaryTag[node_end]==0){//non-boundary vertex
                cout<<"t is non-boundary vertex."<<endl;
                //Search the target partition
                DijkstraPartitionMC(node_end, node_start, localDisV);
                ++p_num;
                for(auto it=localDisV.begin();it!=localDisV.end();++it){
                    int bid=it->first;

                    for(cri_i=0; cri_i<num_of_cri; ++cri_i){
                        if(it->second[cri_i]<INF){
                            mc_cost_r[cri_i][bid] = it->second[cri_i];
                            EM_MC_PQueue_r[cri_i]->push(VertexCost(bid, mc_cost_r[cri_i][bid]));
                        }
                    }
                }

            }else{
                cout<<"t is boundary vertex."<<endl;
                for(cri_i=0; cri_i<num_of_cri; ++cri_i) {
                    EM_MC_PQueue_r[cri_i]->push(VertexCost(node_end, 0));
                }
            }
            localDisV.clear();
        }


        for(cri_i=0;cri_i<num_of_cri;++cri_i){
            if(!EM_MC_PQueue[cri_i]->empty()){
                cri_to_cluster[cri_i] = node_to_cluster[EM_MC_PQueue[cri_i]->top().id];
            }
            else{
                cout<<"F. Priority queue of criterion "<<cri_i<<" is empty!"<<endl;
                if(set_cri.find(cri_i)!=set_cri.end()){
                    set_cri.erase(cri_i);
                    mc_finished[cri_i] = true;
                    mc_min_cost[cri_i] = INF;
                    cri_hops[cri_i] = INF;
                    cri_pqueue.update(cri_i, INF);
                    cri_hops_r[cri_i] = INF;
                    cri_pqueue_r.update(cri_i, INF);
                    ++cri_empty;
                }
            }
            if(!EM_MC_PQueue_r[cri_i]->empty()){
                cri_to_cluster_r[cri_i] = node_to_cluster[EM_MC_PQueue_r[cri_i]->top().id];
            }
            else{
                cout<<"R. Priority queue of criterion "<<cri_i<<" is empty!"<<endl;
                if(set_cri.find(cri_i)!=set_cri.end()){
                    set_cri.erase(cri_i);
                    mc_finished[cri_i] = true;
                    mc_min_cost[cri_i] = INF;
                    cri_hops[cri_i] = INF;
                    cri_pqueue.update(cri_i, INF);
                    cri_hops_r[cri_i] = INF;
                    cri_pqueue_r.update(cri_i, INF);
                    ++cri_empty;
                }
            }
        }

        cluster_to_cri[node_to_cluster[node_start]] = 0;
        cluster_to_cri[node_to_cluster[node_end]] = 0;
        //Iteration
        while (!EM_JudgeEmptyBi(EM_MC_PQueue, EM_MC_PQueue_r)) {//for every node in pqueue !mc_finished[0] || !mc_finished[1] || !mc_finished[2]
            set_remain = set_cri;
            /// read by synchronizing the hops
            temp_cluster_id = cri_to_cluster[cri_pqueue.top_value()];
            temp_cluster_id_r = cri_to_cluster_r[cri_pqueue_r.top_value()];
            if(cri_hops[cri_pqueue.top_value()] <= cri_hops_r[cri_pqueue_r.top_value()]){
                if(temp_cluster_id!=pid_end){
//                    cout<<"F. Read partition shortcut "<<temp_cluster_id<<endl;
                    ReadPartitionShortcut_Bi(temp_cluster_id, mc_HotPool, cluster_to_bi, FORWARD, clusterRead);//if in disk, read corresponding partition item_id,
                }else{
//                    cout<<"F. Read partition graph "<<temp_cluster_id<<endl;
                    GraphMCReadCluster_New_Bi(mc_HotPool, temp_cluster_id, cluster_to_bi, FORWARD, clusterRead);//if in disk, read corresponding partition item_id,
                }
                cluster_to_bi[temp_cluster_id] = FORWARD;
                ++p_num;
            }else{
                if(temp_cluster_id_r!=pid_start){
//                    cout<<"R. Read partition shortcut "<<temp_cluster_id_r<<endl;
                    ReadPartitionShortcut_Bi(temp_cluster_id_r, mc_HotPool, cluster_to_bi, REVERSE, clusterRead);//if in disk, read corresponding partition item_id_r,
                }else{
//                    cout<<"R. Read partition graph "<<temp_cluster_id_r<<endl;
                    GraphMCReadCluster_New_Bi(mc_HotPool, temp_cluster_id_r, cluster_to_bi, REVERSE, clusterRead);//if in disk, read corresponding partition item_id_r,
                }
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
                        cout<<"T2 Minimal cost of criterion "<<cri_i<<": "<<mc_min_cost[cri_i]<<endl;
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
                                if (mc_HotPool.HotPools[partition_id][i].getW(cri_i)==INF100){
//                                    cout<<"F I. Continue between "<<item_id<<" and "<<temp_id<<"("<<i<<") "<<INF100<<", criterion "<<cri_i<<", partition "<<partition_id<<endl;
                                    continue;
                                }

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
                                if (mc_HotPool.MCEdges_EM[i].getW(cri_i)==INF100){
                                    cout<<"F S. Continue between "<<item_id<<" and "<<temp_id<<"("<<i<<") "<<INF100<<", criterion "<<cri_i<<", partition "<<partition_id<<endl;
                                    continue;
                                }

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
                            cout << "F. Wrong partition status!!! " <<item_id<<"("<<get<2>(EMEdgesIndex_Bi[item_id])<<") "<<partition_id<<"("<<get<0>(mc_HotPool.clusterStatus[partition_id])<<") "<<cri_i<< endl; exit(1);
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
                                    for(auto it=cluster_to_node[partition_id].begin();it!=cluster_to_node[partition_id].end();++it){
                                        int id=*it;
                                        get<2>(EMEdgesIndex_Bi[id]) = -1;//hard delete
                                        get<3>(EMEdgesIndex_Bi[id]) = -1;
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
                        cout<<"T2 Minimal cost of criterion "<<cri_i<<": "<<mc_min_cost[cri_i]<<endl;
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
                                if (mc_HotPool.HotPools[partition_id_r][i].getW(cri_i)==INF100){
//                                    cout<<"R I. Continue between "<<item_id_r<<" and "<<temp_id_r<<"("<<i<<") "<<INF100<<", criterion "<<cri_i<<", partition "<<partition_id_r<<endl;
                                    continue;
                                }

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
                                if (mc_HotPool.MCEdges_EM[i].getW(cri_i)==INF100){
                                    cout<<"R S. Continue between "<<item_id_r<<" and "<<temp_id_r<<"("<<i<<") "<<INF100<<", criterion "<<cri_i<<" partition "<<partition_id_r<<endl;
                                    continue;
                                }

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
                            cout << "R. Wrong partition status!!! " << item_id_r<<"("<<get<2>(EMEdgesIndex_Bi[item_id_r])<<") "<<partition_id_r<<"("<< get<0>(mc_HotPool.clusterStatus[partition_id_r])<<") "<<cri_i<<endl; exit(1);
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
                                    for(auto it=cluster_to_node[partition_id_r].begin();it!=cluster_to_node[partition_id_r].end();++it){
                                        get<2>(EMEdgesIndex_Bi[*it]) = -1;//hard delete
                                        get<3>(EMEdgesIndex_Bi[*it]) = -1;
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
                if (EM_MC_PQueue[cri_i]->top().cost + EM_MC_PQueue_r[cri_i]->top().cost >= mc_min_cost[cri_i]) {//condition of termination !EM_MC_PQueue[cri_i]->empty() && !EM_MC_PQueue_r[cri_i]->empty() &&
                    mc_finished[cri_i] = true;
                    set_cri.erase(cri_i);
                    cri_hops[cri_i] = INF;
                    cri_pqueue.update(cri_i, INF);
                    cri_hops_r[cri_i] = INF;
                    cri_pqueue_r.update(cri_i, INF);
//                    cout<<node_start<<" "<<node_end<<endl;
                    cout<<"Minimal cost of criterion "<<cri_i<<": "<<mc_min_cost[cri_i]<<endl;
//                    cout<<EM_MC_PQueue[cri_i]->top().id<<" ("<< EM_MC_PQueue[cri_i]->top().cost<<") ; "<<EM_MC_PQueue_r[cri_i]->top().id<<" ("<<EM_MC_PQueue_r[cri_i]->top().cost<<")"<<endl;
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
    //BMHPS, new, 2023-08-22
    /*void EMMCGraph::EM_MC_ShortcutSearch_Bi(int node_start, int node_end){
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
        int cri_i=0;
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
            mc_cost_r[i][node_end] = 0;
            EM_MC_PQueue[i]->push(VertexCost(node_start, 0));
            EM_MC_PQueue_r[i]->push(VertexCost(node_end, 0));
            cri_to_cluster[i] = node_to_cluster[node_start];
            cri_to_cluster_r[i] = node_to_cluster[node_end];
            cri_pqueue.update(i,0);
            cri_pqueue_r.update(i,0);
        }

        //Search the source partition
        int pid_start=node_to_cluster[node_start];
        int pid_end=node_to_cluster[node_end];
        cout<<"pid_start: "<<pid_start<<", pid_end: "<<pid_end<<endl;
//        GraphMCReadCluster_New(mc_HotPool, node_start, pid_start, MultiHops);

        vector<pair<NodeId,vector<Distance>>> localDisV;

        bool firstRead=false;
        if(BoundaryTag[node_start]!=0 && BoundaryTag[node_end]!=0){
            cout<<"s and t are both boundary vertices."<<endl;
        }else {
            if(pid_start == pid_end){
                GraphMCReadCluster_New_Bi(mc_HotPool, pid_start, cluster_to_bi, FORWARD, clusterRead);//if in disk, read corresponding partition item_id,
                cluster_to_bi[temp_cluster_id] = FORWARD;
                ++p_num;
            }else{
                GraphMCReadCluster_New_Bi(mc_HotPool, pid_start, cluster_to_bi, FORWARD, clusterRead);//if in disk, read corresponding partition item_id,
                cluster_to_bi[temp_cluster_id] = FORWARD;
                ++p_num;
                GraphMCReadCluster_New_Bi(mc_HotPool, pid_end, cluster_to_bi, REVERSE, clusterRead);//if in disk, read corresponding partition item_id_r,
                cluster_to_bi[temp_cluster_id_r] = REVERSE;
                ++p_num;
            }
            firstRead = true;
        }


        for(cri_i=0;cri_i<num_of_cri;++cri_i){
            cri_to_cluster[cri_i] = node_to_cluster[EM_MC_PQueue[cri_i]->top().id];
            cri_to_cluster_r[cri_i] = node_to_cluster[EM_MC_PQueue_r[cri_i]->top().id];
        }

        cluster_to_cri[node_to_cluster[node_start]] = 0;
        cluster_to_cri[node_to_cluster[node_end]] = 0;
        //Iteration
        while (!EM_JudgeEmptyBi(EM_MC_PQueue, EM_MC_PQueue_r)) {//for every node in pqueue !mc_finished[0] || !mc_finished[1] || !mc_finished[2]
            set_remain = set_cri;
            /// read by synchronizing the hops
            if(!firstRead){
                temp_cluster_id = cri_to_cluster[cri_pqueue.top_value()];
                temp_cluster_id_r = cri_to_cluster_r[cri_pqueue_r.top_value()];
                if(cri_hops[cri_pqueue.top_value()] <= cri_hops_r[cri_pqueue_r.top_value()]){
                    if(temp_cluster_id!=pid_end){
//                    cout<<"F. Read partition shortcut "<<temp_cluster_id<<endl;
                        ReadPartitionShortcut_Bi(temp_cluster_id, mc_HotPool, cluster_to_bi, FORWARD, clusterRead);//if in disk, read corresponding partition item_id,
                    }else{
//                    cout<<"F. Read partition graph "<<temp_cluster_id<<endl;
                        GraphMCReadCluster_New_Bi(mc_HotPool, temp_cluster_id, cluster_to_bi, FORWARD, clusterRead);//if in disk, read corresponding partition item_id,
                    }
                    cluster_to_bi[temp_cluster_id] = FORWARD;
                    ++p_num;
                }else{
                    if(temp_cluster_id_r!=pid_start){
//                    cout<<"R. Read partition shortcut "<<temp_cluster_id_r<<endl;
                        ReadPartitionShortcut_Bi(temp_cluster_id_r, mc_HotPool, cluster_to_bi, REVERSE, clusterRead);//if in disk, read corresponding partition item_id_r,
                    }else{
//                    cout<<"R. Read partition graph "<<temp_cluster_id_r<<endl;
                        GraphMCReadCluster_New_Bi(mc_HotPool, temp_cluster_id_r, cluster_to_bi, REVERSE, clusterRead);//if in disk, read corresponding partition item_id_r,
                    }
                    cluster_to_bi[temp_cluster_id_r] = REVERSE;
                    ++p_num;
                }
            }else{//if has first read already
                firstRead = false;
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
                        cout<<"T2 Minimal cost of criterion "<<cri_i<<": "<<mc_min_cost[cri_i]<<endl;
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
                                if (mc_HotPool.HotPools[partition_id][i].getW(cri_i)==INF100){
//                                    cout<<"F I. Continue between "<<item_id<<" and "<<temp_id<<"("<<i<<") "<<INF100<<", criterion "<<cri_i<<", partition "<<partition_id<<endl;
                                    continue;
                                }

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
                                if (mc_HotPool.MCEdges_EM[i].getW(cri_i)==INF100){
                                    cout<<"F S. Continue between "<<item_id<<" and "<<temp_id<<"("<<i<<") "<<INF100<<", criterion "<<cri_i<<", partition "<<partition_id<<endl;
                                    continue;
                                }

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
                                    for(auto it=cluster_to_node[partition_id].begin();it!=cluster_to_node[partition_id].end();++it){
                                        int id=*it;
                                        get<2>(EMEdgesIndex_Bi[id]) = -1;//hard delete
                                        get<3>(EMEdgesIndex_Bi[id]) = -1;
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
                        cout<<"T2 Minimal cost of criterion "<<cri_i<<": "<<mc_min_cost[cri_i]<<endl;
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
                                if (mc_HotPool.HotPools[partition_id_r][i].getW(cri_i)==INF100){
//                                    cout<<"R I. Continue between "<<item_id_r<<" and "<<temp_id_r<<"("<<i<<") "<<INF100<<", criterion "<<cri_i<<", partition "<<partition_id_r<<endl;
                                    continue;
                                }

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
                                if (mc_HotPool.MCEdges_EM[i].getW(cri_i)==INF100){
                                    cout<<"R S. Continue between "<<item_id_r<<" and "<<temp_id_r<<"("<<i<<") "<<INF100<<", criterion "<<cri_i<<" partition "<<partition_id_r<<endl;
                                    continue;
                                }

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
                                    for(auto it=cluster_to_node[partition_id_r].begin();it!=cluster_to_node[partition_id_r].end();++it){
                                        get<2>(EMEdgesIndex_Bi[*it]) = -1;//hard delete
                                        get<3>(EMEdgesIndex_Bi[*it]) = -1;
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
                if (EM_MC_PQueue[cri_i]->top().cost + EM_MC_PQueue_r[cri_i]->top().cost >= mc_min_cost[cri_i]) {//condition of termination !EM_MC_PQueue[cri_i]->empty() && !EM_MC_PQueue_r[cri_i]->empty() &&
                    mc_finished[cri_i] = true;
                    set_cri.erase(cri_i);
                    cri_hops[cri_i] = INF;
                    cri_pqueue.update(cri_i, INF);
                    cri_hops_r[cri_i] = INF;
                    cri_pqueue_r.update(cri_i, INF);
//                    cout<<node_start<<" "<<node_end<<endl;
                    cout<<"Minimal cost of criterion "<<cri_i<<": "<<mc_min_cost[cri_i]<<endl;
                    cout<<EM_MC_PQueue[cri_i]->top().id<<" ("<< EM_MC_PQueue[cri_i]->top().cost<<") ; "<<EM_MC_PQueue_r[cri_i]->top().id<<" ("<<EM_MC_PQueue_r[cri_i]->top().cost<<")"<<endl;
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
    }*/
    void EMMCGraph::ReadPartitionShortcut_Bi(int target_p_id, MCHotPool<VectorMCEdgesEMTuple_IO_Bi> & mcHotPool, vector<int8_t>& cluster_to_bi, int8_t direction, vector<bool>& clusterRead){
        int evict_p_id = mcHotPool.LRU_IMCluster.put(target_p_id,1);
        //// evict old partition for new one
        if(evict_p_id != -1){//if it is necessary to evict old partition
//            cout<<"There is an eviction! "<< evict_p_id<<endl;
            assert(get<0>(mcHotPool.clusterStatus[evict_p_id]) == 'I');
            unsigned int evict_p_size = get<2>(mcHotPool.clusterStatus[evict_p_id]);
            uint temp_sz = get<1>(mcHotPool.clusterStatus[evict_p_id])*alpha_bi;//
            assert(temp_sz>0);
            assert(evict_p_size>0);

            if(evict_p_size > temp_sz){//if the remaining size is too large or no node remain, evict all
//                get<2>(mcHotPool.clusterStatus[evict_p_id]) = 0;
                get<0>(mcHotPool.clusterStatus[evict_p_id]) = 'D';
//                cluster_to_bi[evict_p_id] = 0;
                for(auto it=cluster_to_node[evict_p_id].begin();it!=cluster_to_node[evict_p_id].end();++it){
                    get<2>(EMEdgesIndex_Bi[*it]) = -1;//hard delete
                    get<3>(EMEdgesIndex_Bi[*it]) = -1;
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
        char filePath[300];
        int u = 0, u_deg;
        int v, weight;
        MCEdgeT temp_edge;

        strcpy(filePath, DataPath);
        strcat(filePath, dataset.c_str());
        strcat(filePath, "/");
        strcat(filePath, aggregateStrategy.c_str());
        strcat(filePath, "/Partitions/");
        strcat(filePath, dataset.c_str());
        strcat(filePath, "_Shortcut_");
        strcat(filePath, partMethod.c_str());
        strcat(filePath, "_");
        strcat(filePath, to_string(PartitionSize).c_str());
        strcat(filePath, "_");
        strcat(filePath, to_string(target_p_id).c_str());
        strcat(filePath, ".bin");

        Timer tt;
        tt.start();
        FILE* file = fopen(filePath, "rb");
        if(!file){
            cout<<"Failed to open file "<<filePath<<endl;
            exit(1);
        }
        assert(file != NULL);
        ReadBuffer read_io(file);
        bool flag_sizeUpdate = false;
        int temp_id = -1;

        if(get<1>(mcHotPool.clusterStatus[target_p_id]) == 0){//if this is the first read, judged by the original size of the cluster
            flag_sizeUpdate = true;
            while(!read_io.isend) {//if not read out
                read_io.read(&u);//get node id
                temp_edge.ID1 = u;
                if(u != temp_id){
                    get<2>(EMEdgesIndex_Bi[u]) = mcHotPool.HotPools[target_p_id].size();
                    if(temp_id != -1){
                        get<3>(EMEdgesIndex_Bi[temp_id]) = mcHotPool.HotPools[target_p_id].size();
                    }
                    temp_id = u;
                }

                assert(vertex_cri_Bi[u].first == -1);
                assert(vertex_cri_Bi[u].second == -1);
                vertex_cri_Bi[u].first = set_cri.size();//num_of_cri;
                vertex_cri_Bi[u].second = set_cri.size();//num_of_cri;

                read_io.read(&u_deg);//get node degree

//                unordered_set<int> temp_set;

                for (int i = 0; i < u_deg; ++i) {
                    read_io.read(&v);//read end node
//                    if(temp_set.find(v)==temp_set.end()){//not found
//                        temp_set.insert(v);
//                    }else{
//                        cout<<"Repeated neighbor! "<<v<<" "<<temp_edge.ID1<<endl;
//                        exit(1);
//                    }
                    temp_edge.ID2 = v;
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
        }
        else{//if it is not the first read, read all
            cout<<"Shortcut. The second read for partition "<<target_p_id<<endl;
            mcHotPool.HotPools[target_p_id].clear();
            if(direction == FORWARD) {//if the partition now is read by forward search
                if (cluster_to_bi[target_p_id] == REVERSE) {//if it has been read by the other search, update the edge number in memory
                    get<2>(mcHotPool.clusterStatus[target_p_id]) = 0;
                }
                while (!read_io.isend) {//if not read out
                    read_io.read(&u);//get node id
                    temp_edge.ID1 = u;

                    if (u != temp_id) {
                        get<2>(EMEdgesIndex_Bi[u]) = mcHotPool.HotPools[target_p_id].size();

                        if (temp_id != -1) {
                            get<3>(EMEdgesIndex_Bi[temp_id]) = mcHotPool.HotPools[target_p_id].size();
                        }
                        temp_id = u;
                    }

                    read_io.read(&u_deg);//get node degree

                    for (int i = 0; i < u_deg; ++i) {
                        read_io.read(&v);//read end node
                        temp_edge.ID2 = v;
                        for (int j = 0; j < NUM_OF_CRITERIA; ++j) {
                            read_io.read(&weight);// read edge weight
                            if (j < num_of_cri) {
                                temp_edge.putW(j, weight);
                            }
                        }
                        mcHotPool.HotPools[target_p_id].emplace_back(temp_edge);
                    }
                    if (cluster_to_bi[target_p_id] == REVERSE) {//if it has been read by the other search
                        if (!get<0>(EMEdgesIndex_Bi[u]))//if false, if it has not been evicted
                            get<2>(mcHotPool.clusterStatus[target_p_id]) += u_deg;
                    }
                }
            }else if(direction == REVERSE) {//if the partition now is read by reverse search
                if(cluster_to_bi[target_p_id] == FORWARD) {//if it has been read by the other search
                    get<2>(mcHotPool.clusterStatus[target_p_id]) = 0;
                }
                while(!read_io.isend) {//if not read out
                    read_io.read(&u);//get node id
                    temp_edge.ID1 = u;

                    if(u != temp_id){
                        get<2>(EMEdgesIndex_Bi[u]) = mcHotPool.HotPools[target_p_id].size();

                        if(temp_id != -1){
                            get<3>(EMEdgesIndex_Bi[temp_id]) = mcHotPool.HotPools[target_p_id].size();
                        }
                        temp_id = u;
                    }

                    read_io.read(&u_deg);//get node degree

                    for (int i = 0; i < u_deg; ++i) {
                        read_io.read(&v);//read end node
                        temp_edge.ID2 = v;
                        for (int j = 0; j < NUM_OF_CRITERIA; ++j) {
                            read_io.read(&weight);
                            if(j < num_of_cri){
                                temp_edge.putW(j,weight);
                            }
                        }
                        mcHotPool.HotPools[target_p_id].emplace_back(temp_edge);
                    }
                    if(cluster_to_bi[target_p_id] == FORWARD){//if it has been read by the other search
                        if(!get<1>(EMEdgesIndex_Bi[u]))//if false, if it has not been evicted
                            get<2>(mcHotPool.clusterStatus[target_p_id]) += u_deg;
                    }
                }
            }
            io_num += 1;//update #IO
        }

        get<3>(EMEdgesIndex_Bi[u]) = mcHotPool.HotPools[target_p_id].size();
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
        cout << "Partition strategy: " << aggregateStrategy << "." <<endl;
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
