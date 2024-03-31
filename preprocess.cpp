/*
 * Filename:    preprocess.cpp
 * Description: execution entry for graph data preprocessing
 * Created:     12 Jan 2022
 * Author:      Xinjie ZHOU
 */

#include "preprocess.hpp"

using namespace std;

void createDirectory(char* dirPath, struct stat& info){
    if( stat( dirPath, &info ) != 0 ){
        printf( "The directory %s does not exist.\n", dirPath );
        int check=mkdir(dirPath,S_IRWXU);//read,write,execute 0777
        if(!check){
            printf( "Automatically created.\n");
        }else{
            printf("Unable to create directory\n");
            exit(1);
        }
    }else if( info.st_mode & S_IFDIR )  // S_ISDIR() doesn't exist on my windows
        printf( "%s is a directory\n", dirPath );
    else
        printf( "%s is no directory\n", dirPath );
}

int main(int argc, char** argv)
{
    if( argc < 3 || argc > 8 ){
        printf("usage:\n<arg1> name of dataset, e.g. USA\n");
        printf("<arg2> average partition size(KB), e.g. 1024, default: 1024\n");
        printf("<arg3> task type; 1: preprocess without original criterion, 2: graph partitioning, 3: query generation; 4: preprocessing for partitioning methods; 5: postprocess for partitioning methods. default: 1\n");
        printf("<arg4> (optional) partition aggregation strategy; pri: primary strategy, ave: average strategy, wave: weighted average strategy, default: wave\n");
        printf("<arg5> (optional) partition method, e.g. KaHyPar, METIS\n");
        printf("<arg6> (optional) partition number, e.g. 256\n");
        printf("<arg7> (optional) vertex number, e.g. 23947347\n");
        exit(0);
    }


    Timer tt;
    tt.start();

    uint n_num = 4377184;//14081816;

    //controlling variables
    bool ifAverage = true;//whether to use the average edge weight for partitioning
    int task_type = 3;//1: preprocess without original criterion, 2: preprocess with original criterion, 3: only partition; 4: only generate cover-based queries; 5: convert MCEdges to simplified MCEdges file; 6: read file; 7: read and preprocess multiplex data
    int partitionNum=256;

    ioxxl::DiskConfigure();
    gbpre::Graph_pre g_pre;
    g_pre.dataset = "USA";// PopeElection2013 BostonBomb2013 USA Orkut Wiki_en Twitter
    g_pre.criteria = "Distance";
    g_pre.strategy = "W";

    if(argc > 1){
        cout << "argc: " <<argc <<  endl;
        cout << "argv[1] (dataset): " << argv[1] << endl;//dataset
        cout << "argv[2] (partition size): " << argv[2] << endl;//partition size
        g_pre.dataset = argv[1];
        PartitionSize = atoi(argv[2]);

        if(argc > 3){
            cout << "argv[3] (task id): " << argv[3] << endl;//task id
            task_type = atoi(argv[3]);
        }
        if(argc > 4){
            cout << "argv[4] (aggregation strategy): " << argv[4] << endl;//aggregation strategy
            string at = argv[4];
            if (at == "pri"){
                g_pre.strategy = "pri";
            }else if(at == "ave"){
                g_pre.strategy = "ave";
            }else if(at == "wave"){
                g_pre.strategy = "wave";
            }else{
                cout<<"Aggregation strategy is wrong!"<<endl;
                exit(1);
            }
        }
        if(argc > 5){
            cout << "argv[5] (partition method): " << argv[5] << endl;//partition method
            string at = argv[5];
            if (at == "DC"){
                g_pre.partMethod = "DC";
            }else if(at == "METIS"){
                g_pre.partMethod = "Metis";
            }else if(at == "KaHyPar"){
                g_pre.partMethod = "KaHyPar";
            }else{
                cout<<"Partition method is wrong!"<<endl;
                exit(1);
            }
        }
        if(argc > 6){
            cout << "argv[6] (partition number): " << argv[6] << endl;//partition number
            partitionNum = atoi(argv[6]);
        }
        if(argc > 7){
            cout << "argv[7] (vertex number): " << argv[7] << endl;//vertex number
            n_num = atoi(argv[7]);
        }

    }

    //check if the directory exist
    struct stat info;
    char dirPath[300];
    strcpy(dirPath, DataPath);
    strcat(dirPath, g_pre.dataset.c_str());
    strcat(dirPath, "/");
    strcat(dirPath, g_pre.strategy.c_str());
    createDirectory(dirPath,info);
    strcat(dirPath, "/Partitions/");
    createDirectory(dirPath,info);


    if(g_pre.dataset =="Louisiana" || g_pre.dataset =="California" || g_pre.dataset =="USA" || g_pre.dataset =="CTR" || g_pre.dataset == "WUSA"){// || g_pre.dataset =="NewYork"
        g_pre.num_criteria = 2;
    }else{
        g_pre.num_criteria = 1;
    }

    cout << "Dataset: " << g_pre.dataset << endl;
    cout << "Task type: " << task_type << endl;
    cout << "Partition size: " << PartitionSize << " KB." << endl;
    cout << "Partition method: " << g_pre.partMethod <<endl;
    cout << "Partition aggregation strategy: " << g_pre.strategy <<endl;
    cout << "-------------------------" <<endl;
//    exit(0);
    switch (task_type){
        case 1:{//task 1: preprocess without original criterion
            /********* Basic criterion Data Generation *********/
            cout << "--------- Basic criterion data generation ----------" << endl;
            g_pre.BaseCriGenerate(n_num);//read, generate and write, memory usage: 2*|E|
            cout << "------------------ Graph Partitioning --------------------" << endl;
            g_pre.GraphPartitioning();
            //generate random OD pairs
            cout << "--------- OD pairs generation ----------" << endl;
            g_pre.ODpairGenerate(10000);
            //generate OD pairs according to distance
            g_pre.ODpairGenerate_Dis2(100, false);
            cout << "-------------------------------------" << endl;
            break;
        }
        case 2:{//task 2: EM_Dijk preprocess
            /********* Data Preprocessing *********/
//            g_pre.IDMinusOne("/Users/zhouxj/Documents/1-Research/Datasets/PopeElection2013/PopeElection2013_MCEdgesS.txt");
            cout << "------------------ Graph Partitioning --------------------" << endl;
            g_pre.GraphPartitioning();
            //generate random OD pairs
            cout << "--------- OD pairs generation ----------" << endl;
            g_pre.ODpairGenerate(10000);
            //generate OD pairs according to distance
            g_pre.ODpairGenerate_Dis2(100, false);
            cout << "-------------------------------------" << endl;
//            cout << "--------- Preprocessing of EM_Dijk ----------" << endl;
//            g_pre.EMDIJk_preprocess(false);
//            cout << "---------------------------------------------" << endl;
       //    //generate bin file for MCEdges_io
       //    g_pre.text_to_bin();
            break;
        }
        case 3:{//queries generation
            cout << "--------- OD pairs generation ----------" << endl;
//            g_pre.ODpairGenerate(10000);
            g_pre.ODpairGenerate_Dis2(100, true);
            cout << "-------------------------------------" << endl;
            break;
        }
        case 4:{// Preprocessing of other partitioning methods
            cout<<"--------- Preprocessing for partition methods ----------"<<endl;
            g_pre.PartitionPreprocess(g_pre.partMethod);
            cout << "-------------------------------------" << endl;
            break;
        }
        case 5:{// Postprocessing of other partitioning methods
            cout<<"--------- Postprocessing for partition methods ----------"<<endl;
            g_pre.PartitionPostprocess(g_pre.partMethod,partitionNum);
            cout << "-------------------------------------" << endl;
            break;
        }
        default:{
            cout << "Task type is wrong!" << endl;
            //    g_pre.Read(string(DataPath)+g_pre.dataset+"/"+g_pre.dataset+"_original.txt",20);
            break;
        }

    }

    tt.stop();
    cout << "\nOverall time cost:" << tt.GetRuntime() << " seconds" << endl;

    PeakMemory();
    exit(0);
    return 0;
}
