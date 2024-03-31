/*
 * Filename:    mcspem.cpp
 * Description: execution entry for mcsp external-memory graph algorithms
 * Created:     04 March 2022
 * Authors:     Xinjie ZHOU
 */
#include "emgraph.hpp"
//#include <unistd.h>

using namespace gbxxl;

int main(int argc, char** argv)
{
    if( argc < 4 || argc > 14 ){
        printf("usage:\n<arg1> name of dataset, e.g. USA\n");
        printf("<arg2> task type, e.g. 0: single algorithm test, 1: vary parameter, 2: effectiveness, 3: comparison with baselines; 4: vary #criteria; 5: vary memory size; 6: vary partition size, default: 0\n");
        printf("<arg3> (optional) search synchronization algorithm, 0:multi-pass, 1:OHP, 2:MHP, default: 2\n");
        printf("<arg4> (optional) whether bidirectional search, 0:No, 1:Yes, default: 1\n");
        printf("<arg5> (optional) whether shortcut optimization, 0:No, 1:Yes, default: 1\n");
        printf("<arg6> (optional) whether IO optimization, 0:No, 1:Yes, default: 1\n");
        printf("<arg7> (optional) number of criteria, e.g. 3, default: 3\n");
        printf("<arg8> (optional) average partition size(KB), e.g. 1024, default: 1024\n");
        printf("<arg9> (optional) partition aggregation strategy, pri: primary strategy, ave: average strategy, wave: weighted average strategy, default: wave\n");
        printf("<arg10> (optional) partition method, e.g. DC: Determine clustering, METIS, KaHyPar, default: DC\n");
        printf("<arg11> (optional) query type, e.g. 1:Short, 2:Medium, 3:Long, 4: all, default: 4\n");
        printf("<arg12> (optional) thread number for shortcut construction, default: 15\n");
        printf("<arg13> (optional) parameter testing choice, e.g. 1: alpha, 2: mu, default 1\n");
        printf("Note that the data source path and memory ratio mu need to be modified in head.h manually.\n");
        exit(0);
    }


    //output content to file
    streambuf* pOld;
    if(ifFile) {
        ofstream fout("log_results.txt");
        pOld = cout.rdbuf(fout.rdbuf());
    }

    Timer tt;
    tt.start();
    ioxxl::DiskConfigure();

    gbxxl::EMMCGraph mcg;

    int experiment_id = EXP3;//choose experiment type, e.g. 0, EXP1, EXP2, EXP3, EXP4, EXP5
    string dt1 = "USA";// PopeElection2013 BostonBomb2013 USA Orkut Wiki_en Twitter
    mcg.dataset = dt1;
    mcg.memSize = 1792;//512;//5120;//16384;//150;//49152;//
    int run_times = 2;

    if(argc > 1){
        cout << "argc: " <<argc <<  endl;
        cout << "argv[1] (dataset): " << argv[1] << endl;//dataset
        mcg.dataset = argv[1];

        cout << "argv[2] (task): " << argv[2] << endl;//task id//partition size
        experiment_id = atoi(argv[2]);

        if(argc > 3) {
            cout << "argv[3] (synchronization strategy): " << argv[3] << endl;//synchronization algorithm
            mcg.algo = atoi(argv[3]);
        }

        if(argc > 4) {
            cout << "argv[4] (bidirectional search): " << argv[4] << endl;//whether bidirectional
            mcg.ifBidirectional = atoi(argv[4]);
        }

        if(argc > 5) {
            cout << "argv[5] (IO optimization): " << argv[5] << endl;//whether IO optimization
            mcg.ifIOOptimization = atoi(argv[5]);
        }

        if(argc > 6){
            cout << "argv[6] (shortcut optimization): " << argv[6] << endl;//whether shortcut optimization
            mcg.ifShortcut = atoi(argv[6]);
        }

        if(argc > 7){
            cout << "argv[7] (criteria number): " << argv[7] << endl;//criteria number
            num_of_cri = atoi(argv[7]);
        }

        if(argc > 8){
            cout << "argv[8] (partition size): " << argv[8] << endl;//criteria number
            PartitionSize = atoi(argv[8]);
        }

        if(argc > 9){
            cout << "argv[9] (aggregation strategy): " << argv[9] << endl;//aggregation strategy
            string at = argv[9];
            if (at == "pri"){
                mcg.aggregateStrategy = "pri";
            }else if(at == "ave"){
                mcg.aggregateStrategy = "ave";
            }else if(at == "wave"){
                mcg.aggregateStrategy = "wave";
            }else{
                cout<<"Aggregation strategy is wrong! "<< at << endl;
                exit(1);
            }
        }
        if(argc > 10){
            cout << "argv[10] (partition method): " << argv[10] << endl;//query type
            string at = argv[10];
            if (at == "DC"){
                mcg.partMethod = "DC";
            }else if(at == "METIS"){
                mcg.partMethod = "Metis";
            }else if(at == "KaHyPar"){
                mcg.partMethod = "KaHyPar";
            }else{
                cout<<"Partition method is wrong! "<< at<<endl;
                exit(1);
            }
        }

        if(argc > 11){
            cout << "argv[11] (query type): " << argv[11] << endl;//query type
            int at = atoi(argv[11]);
            if (at == 1){
                mcg.query_type = "S";
            }else if(at == 2){
                mcg.query_type = "M";
            }else if(at == 3){
                mcg.query_type = "L";
            }else if(at == 4){
                mcg.query_type = "all";
            }else{
                cout<<"Aggregation strategy is wrong!"<<endl;
                exit(1);
            }
        }


        if(argc > 12){
            cout << "argv[12] (thread number): " << argv[12] << endl;//thread number
            mcg.threadnum=atoi(argv[12]);
        }

        if(argc > 13){
            cout << "argv[13]: " << argv[13] << endl;//alpha and mu
            int at=atoi(argv[13]);
            if(at == 1){
                ifMu = false;
                cout<<"ifMu: "<<ifMu<<endl;
            }else if(at == 2){
                ifMu = true;
                cout<<"ifMu: "<<ifMu<<endl;
            }else{
                cout<<"Wrong parameter! "<< at<<endl;
                exit(1);
            }
        }
    }
    if(mcg.dataset == "BostonBomb2013" || mcg.dataset == "PopeElection2013"){
        mcg.num_criteria = 3;
        NUM_OF_CRITERIA = 3;
    }else{
        mcg.num_criteria = NUM_OF_CRITERIA;
        NUM_OF_CRITERIA = 5;
    }

    ///ad-hoc test
    //mcg.InfoPrint();
//    mcg.MC_Evaluate_IO(query_type);

    cout<<"Hello world!"<<endl;
    cout<<"Aggregation strategy: "<<mcg.aggregateStrategy<<endl;
    cout<<"Thread number: "<<mcg.threadnum<<endl;


//    priority_queue<pair<int,int>, vector<pair<int,int>>, IntPairCompareLess> pqueue;
//    pqueue.push(make_pair(1,0));
//    pqueue.push(make_pair(300,4));
//    pqueue.push(make_pair(32,4));
//    while(!pqueue.empty()){
//        cout<<pqueue.top().first<<" "<<pqueue.top().second<<endl;
//        pqueue.pop();
//    }
//    exit(0);

    ///Experimental evaluation
    switch (experiment_id) {
        case 0:{
//            mcg.GetBoundary();
//            mcg.ShortcutConstruction();
//            mcg.ShortcutSearch()
            mcg.MC_Evaluate_IO(mcg.query_type);
//            mcg.CorrectnessCheck(100);
            break;
        }
        case EXP1:{//Exp1: Parameters setting
            string dt[6] = {"PopeElection2013", "BostonBomb2013", "USA","Orkut","Wiki_en","Twitter"};
            cout << "--------------------------------" << endl;
            cout << "This is Experiment 1!" << endl;
            bool flag_exit = false;

            cout << "!!! The dataset is " << mcg.dataset << endl;
            for (int k = 7; k <= 9; ++k) {
//                mcg.MemoryCheck(k, mcg.memSize, flag_exit);
            }
            if (flag_exit) exit(1);
            cout<<"ifMu: "<<ifMu<<endl;
            if(!ifMu){
                mcg.InfoPrint();
                double alpha = 1;//alpha = 0.6, 0.7, 0.8, 0.9, 1.0
                for(int j=0;j<6;++j){//5
                    cout << "!! All parameters(i.e. alpha_OHP, alpha_MHP, alpha_BMHP) are set to " << alpha << endl;
                    mcg.alpha = alpha; mcg.alpha_multi = alpha; mcg.alpha_bi = alpha;
                    mcg.MC_Evaluate_EXP1(mcg.query_type);
                    alpha -= 0.1;
                }
            }else{
                cout << "!! Mu_OHP and Mu_MHP are set to " << MuForEM << ", while Mu_BMHP is " << MuForEM_Bi << endl;
                mcg.InfoPrint();
                mcg.MC_Evaluate_EXP1(mcg.query_type);
            }
            break;
        }
        case EXP2: {//EXP2: Effectiveness of optimization
            string dt[6] = {"PopeElection2013", "BostonBomb2013", "USA","Orkut","Wiki_en","Twitter"};
            cout << "--------------------------------" << endl;
            cout << "This is Experiment 2!" << endl;
            bool flag_exit = false;
            for (int k = 1; k <= 9; ++k) {
                mcg.MemoryCheck(k, mcg.memSize, flag_exit);
            }
            if (flag_exit) exit(1);
            cout << "!!! The dataset is " << mcg.dataset << endl;
            mcg.InfoPrint();
            mcg.MC_Evaluate_EXP2(mcg.query_type);

            break;
        }
        case EXP3:{//Exp3: Effectiveness of proposed methods
            string dt[6] = {"PopeElection2013", "BostonBomb2013", "USA","Orkut","Wiki_en","Twitter"};//"CTR", ,"UK2007"
            cout << "--------------------------------" << endl;
            cout << "This is Experiment 3!" << endl;
            bool flag_exit = false;
//            for (int k = 2; k < 3; ++k) {
            for (int k = 1; k < 10; ++k) {//10
//                mcg.MemoryCheck(k, mcg.memSize, flag_exit);
            }
            if (flag_exit) exit(1);
            cout << "!!! The dataset is " << mcg.dataset << endl;
            mcg.InfoPrint();
            mcg.MC_Evaluate_EXP3(mcg.query_type);

            break;
        }
        case EXP4:{//Exp4: Varying #criteria
            string dt[2] = {"USA","Orkut"};//
            cout << "--------------------------------" << endl;
            cout << "This is Experiment 4!" << endl;
//            mcg.dataset = dt[1];
            cout << "!!! The dataset is " << mcg.dataset << endl;
            for(int j=2;j<=5;++j){
                bool flag_exit = false;
                num_of_cri = j;
                for (int k = 1; k <= 9; ++k) {//9
//                    mcg.MemoryCheck(k, mcg.memSize, flag_exit);
                }
                if (flag_exit) exit(1);
                mcg.InfoPrint();
                mcg.MC_Evaluate_EXP456(mcg.query_type);
            }

            break;
        }
        case EXP5:{//Exp5: Varying memory size
            string dt[3] = {"BostonBomb2013","USA","Orkut"};;
            cout << "--------------------------------" << endl;
            cout << "This is Experiment 5!" << endl;
            bool flag_exit = false;
            int mem[3][4] = {
                    {256,341,426,512},
                    {1280,1536,1792,2048},
                    {3072,4096,5120,6144},
//                    {6144,8192,10240,12288},
//                    {10240,12288,14336,16384}
            };
            for(int i=2;i<3;++i){//dt->size()
                flag_exit = false;
                mcg.dataset = dt[i];
                cout << "--------------------------------" << endl;
                cout << "!!! The dataset is " << mcg.dataset << endl;
                for(int j=0;j<1;++j){//0
                    mcg.memSize = mem[i][j];
                    cout << "!! Memory limitation: " << mcg.memSize << endl;
                    //check memory
                    for(int k=1;k<=9;++k){
                        mcg.MemoryCheck(k,mem[i][j],flag_exit);
                    }
                    if(flag_exit) exit(1);
                    mcg.InfoPrint();
                    mcg.MC_Evaluate_EXP456(mcg.query_type);
                }
            }
            break;
        }
        case EXP6: {//EXP6: Effect of partition size
            string dt[2] = {"USA","Orkut"};
            int mem_partition[2][4] = {
                    {256,512,768,1024},
                    {1024,1536,2048,2560},
            };
            cout << "--------------------------------" << endl;
            cout << "This is Experiment 6!" << endl;
            int dt_id = 0;
            cout << "!!! The dataset is " << mcg.dataset << endl;
            for (int j = 1; j < 4; ++j) {
                PartitionSize = mem_partition[dt_id][j];
                cout << "!! Partition Size: " << PartitionSize << endl;
                Partition_N = (MemGraph_IO - PageNumberIO / 4) * 1024 / PartitionSize;
                Partition_N_Bi = (MemGraph_IO_Bi - PageNumberIO_Bi / 4) * 1024 / PartitionSize;
                mcg.InfoPrint();
                mcg.MC_Evaluate_EXP456(mcg.query_type);
            }

            break;
        }
        case EXP7: {//EXP7: Effect of partition aggregation strategy
            string dt[2] = {"USA","Orkut"};
            cout << "--------------------------------" << endl;
            cout << "This is Experiment 7!" << endl;
            bool flag_exit = false;
//            for (int k = 2; k < 3; ++k) {
            for (int k = 1; k < 10; ++k) {//10
//                mcg.MemoryCheck(k, mcg.memSize, flag_exit);
            }
            if (flag_exit) exit(1);
            cout << "!!! The dataset is " << mcg.dataset << endl;
            mcg.InfoPrint();
            mcg.MC_Evaluate_EXP456(mcg.query_type);

            break;
        }
        case EXP8:{//Exp8: Results on Flash Memory
//            string dt[4] = {"BostonBomb2013", "USA","Orkut","Wiki_en"};
//            cout << "--------------------------------" << endl;
//            cout << "This is Experiment 6!" << endl;
//            string qtype[3] = {"S","M","L"};
//            for(int i=0;i<dt->size();++i){
//                mcg.dataset = dt[i];
//                cout << "!!! The dataset is " << mcg.dataset << endl;
//                for(int j=0;j<qtype->size();++j){
//                    query_type = qtype[j];
//                    mcg.InfoPrint();
////                    mcg.MC_Evaluate_EXP456(query_type);
//                }
//            }
            break;
        }
        default:
            break;
    }

    tt.stop();
    cout<<"\nTotal execution time: "<<tt.GetRuntime()<<" s."<<endl;
    PeakMemory();

    if(ifFile){
        cout.rdbuf(pOld);
    }
//    while(1);
    return 0;
}
