/*
 * Filename:    main.cpp
 * Description: execution entry for MCSPs algorithms
 * Created:     04 March 2022
 * Authors:     Xinjie ZHOU
 */
#include "emgraph.hpp"
//#include <unistd.h>

using namespace gbxxl;

int main(int argc, char** argv)
{
   if( argc < 5 || argc > 8 ){
       printf("usage:\n<arg1> name of dataset, e.g. USA\n");
       printf("<arg2> average partition size (KB), e.g. 256, default: 256\n");
       printf("<arg3> number of criteria, e.g. 3, default: 3\n");
       printf("<arg4> (optional) algorithm, e.g. 1: OHP, 2: MHP, 3: BMHP, 4: OHP(w), 5: MHP(w), 6: BMHP(w), default: 3\n");
       printf("<arg5> (optional) run time, e.g. 10, default: 100\n");
       printf("<arg6> (optional) query type, e.g. S, M, L, default: S\n");
       printf("<arg7> (optional) memory size (MB), e.g. 1760, default: 1760\n");
       printf("Other parameters have to be modified manually in head.h file. Sorry for the inconvenience.\n");
       exit(0);
   }

    //ofstream fout("log_results.txt");
    //streambuf* pOld = cout.rdbuf(fout.rdbuf());

    cout << "Testing of Stop & Share algorithm!" << endl;

    Timer tt;
    tt.start();
    ioxxl::DiskConfigure();

    gbxxl::EMMCGraph mcg;

    string dt1 = "USA"; // USA Orkut UK2002 Wiki_en
    mcg.dataset = dt1;
    mcg.memSize = 1792; // 1792 5120
    int algo_choice = BiMultiHops;

    if(argc > 1){
        cout << "argc: " <<argc <<  endl;
        cout << "argv[1]: " << argv[1] << endl;//dataset
        cout << "argv[2]: " << argv[2] << endl;//partition size
        cout << "argv[3]: " << argv[3] << endl;//number of criteria

        mcg.dataset = argv[1];
        PartitionSize = atoi(argv[2]);
        num_of_cri = atoi(argv[3]);

        if(argc > 4){//algorithm
            algo_choice = atoi(argv[4]);
        }
        if(argc > 5){//run times
            run_times = atoi(argv[5]);
        }
        if(argc > 6){//query type
            query_type = argv[6];
        }
        if(argc > 7){//memory size
            mcg.memSize = atoi(argv[7]);
        }
    }

    ///ad-hoc test
    bool flag_exit = false;
    mcg.MemoryCheck(algo_choice, mcg.memSize, flag_exit);
    if (flag_exit) exit(1);
    mcg.InfoPrint();
    mcg.MC_Evaluate_IO(query_type,algo_choice);

    tt.stop();
    cout<<"\nTotal execution time: "<<tt.GetRuntime()<<" s."<<endl;

    //cout.rdbuf(pOld);
    return 0;
}
