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
    if( argc < 3 || argc > 7 ){
        printf("usage:\n<arg1> name of dataset, e.g. USA\n");
        printf("<arg2> average partition size(KB), e.g. 1024, default: 1024\n");
        printf("<arg3> (optional) number of testing query, e.g. 20, default: 100\n");
        printf("<arg4> (optional) number of criteria, e.g. 3, default: 3\n");
        printf("<arg5> (optional) partition strategy, e.g. 2: PRI, 2: AVE, 3: WAVE, default 2\n");
        printf("<arg6> (optional) memory size, e.g. 5120, default 1790\n");
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
    mcg.memSize = 1792;//5120;//512;//16384;//150;//49152;//

    if(argc > 1){
        cout << "argc: " <<argc <<  endl;
        cout << "argv[1]: " << argv[1] << endl;//dataset
        cout << "argv[2]: " << argv[2] << endl;//partition size
        mcg.dataset = argv[1];

        PartitionSize = atoi(argv[2]);

        if(argc > 3){
            cout << "argv[3]: " << argv[3] << endl;//run time
            run_times = atoi(argv[3]);
        }

        if(argc > 4){
            cout << "argv[4]: " << argv[4] << endl;//criteria number
            num_of_cri = atoi(argv[4]);
        }

        if(argc > 5){
            cout << "argv[5]: " << argv[5] << endl;//aggregation strategy
            int at = atoi(argv[5]);
            if (at == 1){
                partition_type = "pri";
            }else if(at == 2){
                partition_type = "ave";
            }else if(at == 3){
                partition_type = "wave";
            }else{
                cout<<"Aggregation strategy is wrong!"<<endl;
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

    ///Experimental evaluation

    cout<<"Hello world!"<<endl;

    ///ad-hoc test
    bool flag_exit = false;
    for (int k = 4; k < 10; ++k) {//10
         mcg.MemoryCheck(k, mcg.memSize, flag_exit);
    }
    if (flag_exit) exit(1);
    cout << "!!! The dataset is " << mcg.dataset << endl;
    mcg.InfoPrint();
    mcg.MC_Evaluate_IO(query_type);

    tt.stop();
    cout<<"\nTotal execution time: "<<tt.GetRuntime()<<" s."<<endl;
    PeakMemory();

    if(ifFile){
        cout.rdbuf(pOld);
    }

    return 0;
}
