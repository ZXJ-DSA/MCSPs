/*
 * Filename:    graph.hpp
 * Created:     18 Sep 2021
 * Authors:     Xinjie ZHOU
 */

#ifndef GRAPH_HPP
#define GRAPH_HPP

#include "head.h"
#include "graph.h"

namespace gbstd {

    using namespace std;

    string Dij = "Dijkstra";
    string BiDij = "BiDijkstra";
    string Ast = "Astar";

    /*** Basic Functions for Data Load ***/
    //read weighted graph
    void MCGraph::ReadGraph_W() {

        //normal distribution generator
        Timer tt;
        tt.start();
        string lineStr;
        string line_symbol;
        string temp_str;
        int ID1, ID2, weight;

        string r_file = string(DataPath)  + dataset + "/" + dataset;
//    string r_file = string(DataPath)  + dataset + "/" + dataset + "_Processed";

        ifstream inFile(r_file, ios::in);
        if (!inFile) { // if not exist
            cout << "Fail to open file" << r_file << endl;
            exit(1);
        }

        /// read graph and recording the degree of vertices
        cout << dataset<<" graph edges Data loading..." << endl;
        inFile >> node_num >> edge_num;
        vector<int> temp_w;
        MCEdges.resize(node_num);
        while (inFile) {//read each line
            inFile >> ID1 >> ID2 >> weight;
            assert(weight > 0);
            temp_w.clear();
            temp_w.push_back(weight);temp_w.push_back(0);temp_w.push_back(0);
            MCEdges[ID1].push_back(MCEdge(ID2,temp_w));
//        Nodes[ID2].adjnodes.push_back( ID1 );
//        Nodes[ID2].adjweight.push_back( weight );
        }
        inFile.close();

        tt.stop();
        cout << "The time used for data reading: " << tt.GetRuntime() << " s." << endl;
        cout << "Node number: " << node_num << ", Edge number: " << edge_num << endl;
        cout << "--------------------" << endl;

    }

    //Function for reading DIMACS-style graph edges data
    void MCGraph::MCReadGraph_Edges(const string &filename) {
        Timer tt;
        tt.start();
        string line_symbol; //line symbol in the beginning of each line
        string temp_str;    //temp variable of string
        int num_line = 0;   //recording the number of lines
        int num_show = 6;   //line number of graph data to be printed
        int num_cri = 5;
        int ID1, ID2, weight;

        //Open file
        ifstream inFile(filename, ios::in);
        if (!inFile) {
            cout << "File opening failed." << endl;
            exit(1);
        }
        cout << "Edges Data loading..." << endl;
        while (inFile) {
            inFile >> line_symbol;
            if (line_symbol == "p") {//read graph basic information
                inFile >> temp_str;
                if (temp_str == "sp") {
                    inFile >> this->node_num >> this->edge_num;
                    this->Edges.assign(node_num, vector<Edge>());
                } else if (temp_str == "criteria") {
                    inFile >> num_cri;
                    cout << "The criteria includes: ";
                    for (int i = 0; i < num_cri; ++i) {
                        inFile >> temp_str;
                        cout << temp_str << " ";
                        mc_criteria.emplace_back(temp_str);
                    }
                    cout << endl;
                }
            } else if (line_symbol == "a") {//read graph data
                inFile >> ID1 >> ID2;
                for (int i = 0; i < num_cri; ++i) {//get edge weights
                    inFile >> weight;
                    if (i == sc_i) {
                        this->Edges[ID1 - 1].emplace_back(Edge(ID2 - 1, weight));
                    }
                }
                if (num_line < num_show) {
                    cout << ID1 << "\t" << ID2 << "\t" << weight << endl;
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
    }
    void MCGraph::MCReadGraph_MCEdgesS(const string& filename){//for MCEdgesMap
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
        MCEdges.assign(node_num, vector<MCEdge>());//allocate space for MCEdges
        cout << "There are " << num_cri << " criteria in this graph." << endl;
        cout << "The criteria includes: ";
        for (i = 0; i < num_cri; ++i) {
            inFile >> temp_str;
            cout << temp_str << " ";
            if (i < num_of_cri)
                mc_criteria.emplace_back(temp_str);
        }
        cout << endl;

        while(inFile){
            inFile >> ID1 >> ID2;
            temp_w.clear();
            for(i=0;i<num_cri;++i) {//get edge weights
                inFile >> weight;
                temp_w.push_back(weight);
            }
            MCEdges[ID1 - 1].emplace_back(MCEdge(ID2 - 1, temp_w));
            if (num_line < num_show) {//show some edges
                cout << ID1 << "\t" << ID2 ;
                for(auto it=temp_w.begin();it!=temp_w.end();++it){
                    cout << "\t" << *it;
                }
                cout << endl;
            }
            else if (num_line == num_show) {
                cout << "..." << endl;
            }
            ++num_line;
        }

        cout << "Data loaded. ";
        inFile.close();
        tt.stop();
        cout << "The number used for data loading: " << tt.GetRuntime() << " s." << endl;
        cout << "Number of nodes:" << node_num << endl;
        cout << "Number of edges:" << edge_num << endl;
        cout << "--------------------" << endl;
    }
    //one-off edges reading for MC-graph
    void MCGraph::MCReadGraph_MCEdges(string filename) {
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
    }

    //Initiate graph for in-memory algorithms
    void MCGraph::Initiation() {
        set_visited.clear();
    }
    /*** In-memory SP algorithms ***/
    //Function of in-memory Dijkstra's algorithm
    Distance MCGraph::Dijkstra(NodeId node_start, NodeId node_end) { // second version, powered by benchmark::heap
        if (node_start == node_end) return 0;
        benchmark::heap<2, NodeId, Distance> pqueue(node_num);
        NodeId item_id, temp_id;
        Distance item_dis, temp_dis;
        vector<bool> closed(node_num, false); //flag vector of whether closed
        vector<Distance> cost(node_num, INF);   //vector of cost
        vector<NodeId> pre(node_num, -1);       //vector of predecessor id
        Distance min_cost = INF;
        //Initiation of start node
        cost[node_start] = 0;//cost of start node
        pqueue.update(node_start, 0);
        visited_number = 0;
        //Iteration
        while (!pqueue.empty()) {//for every node in pqueue
            pqueue.extract_min(item_id, item_dis);// top and delete min item
            if (item_id == node_end) {//if reach target node
                min_cost = cost[item_id];
//                set_visited.insert(item_id);
                ++visited_number;
                cout <<node_start<<" "<<node_end<<" : ";
                cout<<"Minimal cost of criterion " << sc_i << ": "<<min_cost<<endl;
                break;
            }
            //relaxation
//            set_visited.insert(item_id);//update the visited set for this criteria
            ++visited_number;//update #visited_number
            for (auto it = MCEdges[item_id].begin(); it != MCEdges[item_id].end(); ++it) {
                temp_id = it->u;
                if (closed[temp_id])//if closed
                    continue;
                temp_dis = item_dis + it->w[sc_i];
                if (cost[temp_id] > temp_dis) {//slack operation
                    cost[temp_id] = temp_dis;
                    pre[temp_id] = item_id;
                    pqueue.update(temp_id, temp_dis);
                }
            }
            closed[item_id] = true;
        }
//        Dij_getPath(pre,node_start,node_end);//traverse the pre vector to get the shortest path
        return min_cost;
    }

    //Function of in-memory Bidirectional Dijkstra's algorithm
    Distance MCGraph::BiDijkstra(NodeId node_start, NodeId node_end) // second version, powered by benchmark::heap
    {
        if (node_start == node_end) return 0;
        benchmark::heap<2, NodeId, Distance> pqueue(node_num), pqueue_b(node_num);
        NodeId item_id, item_id_b, temp_id, temp_id_b;
        Distance item_dis, item_dis_b, temp_dis, temp_dis_b;
        NodeId terminate_id;//termination id of bi-dijkstra
        vector<Distance> cost(node_num, INF);         //cost for current node to target node
        vector<bool> closed(node_num, false); //flag of whether having been closed
        vector<Distance> cost_b(node_num, INF);         //cost for current node to target node
        vector<bool> closed_b(node_num, false); //flag of whether having been closed
        vector<NodeId> pre(node_num, -1); //forward search vector of predecessor vertex id
        vector<NodeId> pre_b(node_num, -1); //backward search vector of predecessor vertex id
        Distance min_cost = INF;

        //Initiation of start node
        cost[node_start] = 0;//cost of start node
        closed[node_start] = true;
        cost_b[node_end] = 0;//cost of start node
        closed_b[node_end] = true;
        pre[node_start] = node_start;
        pre_b[node_end] = node_end;
        pqueue.update(node_start, 0);
        pqueue_b.update(node_end, 0);
        visited_number = 0;
        //Iteration
        while (!pqueue.empty() && !pqueue_b.empty()) {//if either priority queue is not empty
            if (pqueue.top_key() + pqueue_b.top_key() >= min_cost) {//condition of termination
                cout<<"Minimal cost of criterion " << sc_i << ": "<<min_cost<<endl;
                break;
            }
            /*-------Forward Search-------*/
            pqueue.extract_min(item_id, item_dis);// top min item
            closed[item_id] = true;
//            set_visited.insert(item_id);//update the visited set for this criteria
            ++visited_number;//update #visited_number
            //relax the adjacent nodes of node item
            for (auto it = MCEdges[item_id].begin(); it != MCEdges[item_id].end(); ++it) {
                temp_id = it->u;
                temp_dis = item_dis + it->w[sc_i];
                if (!closed[temp_id]) {//if not visited
                    if (cost[temp_id] > temp_dis) {//relaxation
                        cost[temp_id] = temp_dis;
                        pre[temp_id] = item_id;
                        pqueue.update(temp_id, cost[temp_id]);
                    }
                }
                if (closed_b[temp_id] && temp_dis + cost_b[temp_id] < min_cost) {
                    min_cost = temp_dis + cost_b[temp_id];
                    terminate_id = temp_id;
                    cout<<"Forward! min_cost: "<<min_cost<<", terminate_id: "<<terminate_id<<endl;
                }
            }
            /*-------Reverse Search-------*/
            pqueue_b.extract_min(item_id_b, item_dis_b);// top min item
            closed_b[item_id_b] = true;
//            set_visited.insert(item_id_b);//update the visited set for this criteria
            ++visited_number;//update #visited_number
            //relax the adjacent nodes of node item
            for (auto it = MCEdges[item_id_b].begin(); it != MCEdges[item_id_b].end(); ++it) {
                temp_id_b = it->u;
                temp_dis_b = item_dis_b + it->w[sc_i];
                if (!closed_b[temp_id_b]) {//if not closed
                    if (cost_b[temp_id_b] > temp_dis_b) {//slack operation
                        cost_b[temp_id_b] = temp_dis_b;
                        pre_b[temp_id_b] = item_id_b;
                        pqueue_b.update(temp_id_b, cost_b[temp_id_b]);
                    }
                }
                if (closed[temp_id_b] && temp_dis_b + cost[temp_id_b] < min_cost) {
                    min_cost = temp_dis_b + cost[temp_id_b];
                    terminate_id = temp_id_b;
                    cout<<"Reverse! min_cost: "<<min_cost<<", terminate_id: "<<terminate_id<<endl;
                }
            }
        }
//        BiDij_getPath(pre,pre_b,node_start,node_end,terminate_id);//traverse the pre vector to get the shortest path
        return min_cost;
    }
    //Function of in-memory Dijkstra's algorithm with preprocessing
    /*Distance MCGraph::DijkstraCH(NodeId node_start, NodeId node_end) { // second version, powered by benchmark::heap
        if (node_start == node_end) return 0;
        benchmark::heap<2, NodeId, Distance> pqueue(node_num);
        NodeId item_id, temp_id;
        Distance item_dis, temp_dis;
        vector<bool> closed(node_num, false); //flag vector of whether closed
        Distance min_cost = INF;
        //Initiation of start node

        for (int i = 0; i < num_of_cri; ++i) {
            mc_cost[i][node_start] = 0;//cost of start node
        }
        pqueue.update(node_start, 0);
        visited_number = 0;
        //Iteration
        while (!pqueue.empty()) {//for every node in pqueue
            pqueue.extract_min(item_id, item_dis);// top and delete min item
            if (item_id == node_end) {//if reach target node
                min_cost = mc_cost[sc_i][item_id];
//                set_visited.insert(item_id);
                ++visited_number;
                break;
            }
            closed[item_id] = true;
            //relaxation
//            set_visited.insert(item_id);//update the visited set for this criteria
            ++visited_number;//update #visited_number
            for (auto it = MCEdges[item_id].begin(); it != MCEdges[item_id].end(); ++it) {
                temp_id = it->u;
                if (closed[temp_id])//if closed
                    continue;
                temp_dis = item_dis + it->w[sc_i];
                if (mc_cost[sc_i][temp_id] > temp_dis) {//slack operation
                    mc_cost[sc_i][temp_id] = temp_dis;
                    mc_pre[sc_i][temp_id] = item_id;
                    pqueue.update(temp_id, temp_dis);
                }
                else if(!pqueue.contains(temp_id)){//if mc_cost is smaller and pqueue does not contain element, update
                    pqueue.update(temp_id, mc_cost[sc_i][temp_id]);
                }
                //CH-style preprocessing
                if (sc_i + 1 < num_of_cri) {
                    for (int j = sc_i + 1; j < num_of_cri; ++j) {
                        temp_dis = mc_cost[j][item_id] + it->w[j];
                        if (mc_cost[j][temp_id] > temp_dis) {//slack operation
                            mc_cost[j][temp_id] = temp_dis;
                            mc_pre[j][temp_id] = item_id;
                        }
                    }
                }

            }
            closed[item_id] = true;
        }
//        Dij_getPath(mc_pre[sc_i],node_start,node_end);//traverse the pre vector to get the shortest path
        return min_cost;
    }*/

    //Function for multipass algorithms
    double MCGraph::MC_Multipass(const string &qtype, bool ifIO, int algo_choice)
    {
        string LongDis, MediumDis, ShortDis, RandomDis;
        double total_time = 0;
        ///File path of OD pairs
        LongDis = string(DataPath) +  dataset + "/" + dataset + "_OD_LongDis_"+ to_string(PartitionSize) + ".txt";
        MediumDis = string(DataPath) +  dataset + "/" + dataset + "_OD_MediumDis_"+ to_string(PartitionSize) + ".txt";
        ShortDis = string(DataPath) +  dataset + "/" + dataset + "_OD_ShortDis_"+ to_string(PartitionSize) + ".txt";
        RandomDis = string(DataPath) +  dataset + "/" + dataset + "_OD_Random_"+ to_string(PartitionSize) + ".txt";
        ///Shortest path query processing
        if (qtype == "S") {
            QueryType = "S";
            total_time = MC_Multipass_one(ShortDis, ifIO, algo_choice);//Efficiency evaluation on short distance
        } else if (qtype == "M") {
            QueryType = "M";
            total_time = MC_Multipass_one(MediumDis, ifIO, algo_choice);//Efficiency evaluation on medium distance
        } else if (qtype == "L") {
            QueryType = "L";
            total_time = MC_Multipass_one(LongDis, ifIO, algo_choice);//Efficiency evaluation on long distance
        } else if (qtype == "all") {
            QueryType = "S";
            run_times = 10;
            total_time += MC_Multipass_one(ShortDis, ifIO, algo_choice);
            QueryType = "M";
            run_times = 3;
            total_time += MC_Multipass_one(MediumDis, ifIO, algo_choice);
            QueryType = "L";
            run_times = 3;
            total_time += MC_Multipass_one(LongDis, ifIO, algo_choice);
        } else if (qtype == "R") {
            QueryType = "R";
            total_time = MC_Multipass_one(RandomDis, ifIO, algo_choice);//Efficiency evaluation on long distance
        }
        return total_time;
    }

    //Function for multi-pass MCSP
    double MCGraph::MC_Multipass_one(const string &filename, bool ifIO, int algo_choice) {
        string fname;
        vector<pair<int, int>> ODpair;
        int num, ID1, ID2;
        double ave_time = 0;//average time
        vector<Distance> ave_cost(num_of_cri, 0);
        query_time = 0;
        visited_number = 0;
        uint ave_visited = 0;   //average number of visited vertices
        uint ave_union = 0;     //the union number of visited vertices

        ///Read OD pairs
        ifstream inFile(filename, ios::in);//
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

        ///Query processing one by one
        for (int i = 7; i < 8; ++i) {//run_times
            if (algo_choice == DIJKSTRACH) {
                mc_cost.assign(num_of_cri, vector<Distance>(node_num, INF));
                mc_pre.assign(num_of_cri, vector<NodeId>(node_num, -1));
            }
            ID1 = ODpair[i].first;
            ID2 = ODpair[i].second;
            cout << "Query: " << ID1 - 1 << " " << ID2 - 1 << endl;
            Initiation();
            Timer tt2;
            tt2.start();
            ///Shortest path query one criterion by one criterion
            for (sc_i = 2; sc_i < num_of_cri; ++sc_i) {
                Timer tt1;
                tt1.start();
                if (algo_choice == DIJKSTRA) {
                    ave_cost[sc_i] += Dijkstra(ID1 - 1, ID2 - 1);
                } else if (algo_choice == BI_DIJKSTRA) {
                    ave_cost[sc_i] += BiDijkstra(ID1 - 1, ID2 - 1);
                } else if (algo_choice == DIJKSTRACH) {
//                    ave_cost[sc_i] += DijkstraCH(ID1 - 1, ID2 - 1);
                }
                tt1.stop();
                query_time += tt1.GetRuntime() * 1000;//ms
                ave_visited += visited_number;
            }
            tt2.stop();
            cout << "Query time: " << tt2.GetRuntime() << "s." <<endl;
            //record information
            ave_union += set_visited.size();
        }

        ave_time = query_time / run_times;
        ave_visited /= run_times;
        ave_union /= run_times;
        if (algo_choice == DIJKSTRA) {
            cout << "Average Dijkstra performance of MultiPass (on " << run_times << " " << QueryType
                 << " OD pairs) is shown below." << endl;
        } else if (algo_choice == BI_DIJKSTRA) {
            cout << "Average BiDijkstra performance of MultiPass (on " << run_times << " " << QueryType
                 << " OD pairs) is shown below." << endl;
        } else if (algo_choice == DIJKSTRACH) {
            cout << "Average DijkstraCH performance of MultiPass (on " << run_times << " " << QueryType
                 << " OD pairs) is shown below." << endl;
        }
        for (int i = 0; i < num_of_cri; ++i) {
            ave_cost[i] /= run_times;
            cout << "Average distance of criteria " << mc_criteria[i] << ": " << ave_cost[i] << endl;
        }
        cout << "! Average query time: " << ave_time << " ms." << endl;
        cout << "! Average number of nodes visited: " << ave_visited << endl;
        cout << "! Union visited vertices of all criteria: " << ave_union << endl;
        return query_time;
    }

    //MCSP Evaluation
    void MCGraph::MC_Evaluate(string qtype) {
        /*--variables about time record--*/
        double time_multipass = 0;
        double time_multipass_bi = 0;
//        double time_multipassplus = 0;

        /*** Read graph data ***/
        string gr_Edge = string(DataPath)  + dataset + "/" + dataset + "_MCEdges.txt";
//        MCReadGraph_MCEdges(gr_Edge);
        gr_Edge = string(DataPath)  + dataset + "/" + dataset + "_MCEdgesS_"+to_string(PartitionSize)+".txt";
        MCReadGraph_MCEdgesS(gr_Edge);
        /*** MultiPass algorithm ***/
//        ///Dijkstra's algorithm
//        cout << "----- Below are results of Dijkstra's MultiPass algorithm -----" << endl;
//        time_multipass = MC_Multipass(qtype, false, DIJKSTRA);
//        cout << "-----------------------------------------------\n" << endl;
        ///BiDijkstra's algorithm
        cout << "----- Below are results of BiDijkstra's MultiPass algorithm -----" << endl;
        time_multipass_bi = MC_Multipass(qtype, false, BI_DIJKSTRA);
        cout << "-----------------------------------------------\n" << endl;
        /*** MultiPass+ algorithms ***/
//        cout << "----- Below are results of MultiPass+ algorithm -----" << endl;
//        time_multipassplus = MC_Multipass(qtype, false, DIJKSTRACH);
//        cout << "-----------------------------------------------\n" << endl;
        cout << "The total run time of Dijkstra's MultiPass algorithm is: " << time_multipass / 1000 << " s." << endl;
        cout << "The total run time of BiDijkstra's MultiPass algorithm is: " << time_multipass_bi / 1000 << " s."
             << endl;
//        cout << "The total run time of MultiPass+ algorithm is: " << time_multipassplus / 1000 << " s." << endl;
    }



    /*** Supportive functions for SP algorithms ***/
    //Function for path showing
    void MCGraph::Pathshow(list<int> list_sp) {
        for (auto it = list_sp.begin(); it != list_sp.end(); ++it) {
            if (++it == list_sp.end()) {
                it--;
                cout << *it + 1 << endl;
            } else {
                it--;
                cout << *it + 1 << "-->";
            }
        }
    }

    //Function for traverse parental vertices to get the shortest path
    list<int> MCGraph::Dij_getPath(vector<NodeId> &pre, NodeId ID1, NodeId ID2) {
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
    list<int> MCGraph::BiDij_getPath(vector<NodeId> &pre, vector<NodeId> &pre_b, NodeId ID1, NodeId ID2, NodeId terminate_id) {
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
#endif //GRAPH_HPP
