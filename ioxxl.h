/*
 * Filename:    ioxxl.h
 * Description: head file of basic data structures and functions based on stxxl
 * Created:     05 March 2022
 * Authors:     Xinjie ZHOU
 */

#ifndef MCSPS_IOXXL_H
#define MCSPS_IOXXL_H

#include "head.h"
#include <stxxl/vector>
#include <stxxl/priority_queue>
#include <stxxl/sorter>
#include <stxxl/stack>
#include <stxxl/stats>
#include <stxxl/queue>

namespace ioxxl{
    using namespace std;

    /*** Sorter ***/
    typedef stxxl::sorter<EdgePair,EdgePairComparator> EdgeSorter;
    typedef stxxl::sorter<pair<NodeId,uint>,PairComparator> RankSorter;
    typedef stxxl::sorter<EdgePairW,EdgePairWComparator> EdgeWSorter;
    /*** External vector ***/
    typedef stxxl::VECTOR_GENERATOR<uint>::result VectorUint;//262144bytes=256KB,32*256KB=8MB ,1,16,Block_SZ
//    typedef stxxl::VECTOR_GENERATOR<pair<NodeId, Edge>,BlockPerPage,128,Block_SZ>::result VectorEdgesEM;//1*128*256KB=32MB
//    typedef stxxl::VECTOR_GENERATOR<pair<NodeId, MCEdge>,BlockPerPage,128,Block_SZ>::result VectorMCEdgesEM;//4*8*2MB=64MB
    ///for algorithms with io optimization
    //for EM_Dijk
    typedef stxxl::VECTOR_GENERATOR<EdgePairW,BlockPerPage, PageNumberEMDijk,Block_SZ>::result VectorEdges_EMDijk;   //4*32*256KB=32MB
    //for OneHop, MultiHops
    typedef stxxl::VECTOR_GENERATOR<MCEdgeT,BlockPerPage,PageNumberIO,Block_SZ>::result VectorMCEdgesEMTuple_IO; //4*32*256KB=32MB
    //for BiMultiHops
    typedef stxxl::VECTOR_GENERATOR<MCEdgeT,BlockPerPage,PageNumberIO_Bi,Block_SZ>::result VectorMCEdgesEMTuple_IO_Bi; //4*32*256KB=32MB
    ///for algorithms without io optimization
    //for Dijkstra_IO
    typedef stxxl::VECTOR_GENERATOR<MCEdgeT,BlockPerPage,PageNumber_DijkIO,Block_SZ>::result VectorMCEdgesEMTuple_DijkIO;
    //for BiDijkstra_IO
    typedef stxxl::VECTOR_GENERATOR<MCEdgeT,BlockPerPage,PageNumber_BiDijkIO,Block_SZ>::result VectorMCEdgesEMTuple_BiDijkIO;
    //for OneHop_NoIO and MultiHops_NoIO
    typedef stxxl::VECTOR_GENERATOR<MCEdgeT,BlockPerPage,PageNumber_NoIO,Block_SZ>::result VectorMCEdgesEMTuple_NoIO;
    //for BiMultiHops_NoIO
    typedef stxxl::VECTOR_GENERATOR<MCEdgeT,BlockPerPage,PageNumber_NoIO_Bi,Block_SZ>::result VectorMCEdgesEMTuple_NoIO_Bi;     //4*PageNumber*256KB = BASIC_SZ for algorithms without io optimization
//    typedef stxxl::VECTOR_GENERATOR<EdgePair>::result VectorEdgePairs;
    /*** Priority Queue ***/
    //overall memory cost for one priority queue is about 50MB in practice
    typedef stxxl::PRIORITY_QUEUE_GENERATOR<VertexCost, PQCompareLess, MEMORY_FOR_PRIORITY_QUEUE, PRIORITY_QUEUE_MAX_SIZE/sizeof(VertexCost)>::result PriorityQueue;
    stxxl::read_write_pool<PriorityQueue::block_type> PQ_Pool((mem_for_pools / 2) / PriorityQueue::block_type::raw_size, (mem_for_pools / 2) / PriorityQueue::block_type::raw_size);
    typedef stxxl::PRIORITY_QUEUE_GENERATOR<pair<uint,EdgePair>, PQEdgePairCompareLess, MEMORY_FOR_PRIORITY_QUEUE, PRIORITY_QUEUE_MAX_SIZE/sizeof(pair<uint,EdgePair>)>::result PriorityQueueEdgePair;
    stxxl::read_write_pool<PriorityQueueEdgePair::block_type> PQEdgePair_Pool((mem_for_pools / 2) / PriorityQueueEdgePair::block_type::raw_size, (mem_for_pools / 2) / PriorityQueueEdgePair::block_type::raw_size);
    /*** Queue ***/
    typedef stxxl::queue<VertexCost> em_queue;//memory cost: about 25MB per queue, 8MB each
    //Function for configuring the disk parameters
    void DiskConfigure() {
        // get uninitialized config singleton
        stxxl::config* cfg = stxxl::config::get_instance();
        // create a disk_config structure.
        cfg->add_disk(stxxl::disk_config(disk_conf));
    }

    class Stats {
    public:
        Stats(double sec, double readTime, double writeTime, double pioTime, double ioWaitTime, unsigned noOfReads, unsigned noOfWrites, long long int readVolume, long long int writeVolume, unsigned aveBlockSizeR, unsigned aveBlockSizeW) :m_sec(sec), m_readTime(readTime), m_writeTime(writeTime), m_pioTime(pioTime), m_ioWaitTime(ioWaitTime), m_noOfReads(noOfReads), m_noOfWrites(noOfWrites), m_readVolume(readVolume), m_writeVolume(writeVolume), m_aveBlockSizeR(aveBlockSizeR), m_aveBlockSizeW(aveBlockSizeW) {};
        Stats() :m_sec(0), m_readTime(0), m_writeTime(0), m_pioTime(0), m_ioWaitTime(0), m_noOfReads(0), m_noOfWrites(0), m_readVolume(0), m_writeVolume(0), m_aveBlockSizeR(0), m_aveBlockSizeW(0) {};
        Stats& operator -=(const Stats& a) {
            m_sec -= a.m_sec;
            m_readTime -= a.m_readTime;
            m_writeTime -= a.m_writeTime;
            m_pioTime -= a.m_pioTime;
            m_ioWaitTime -= a.m_ioWaitTime;
            m_noOfReads -= a.m_noOfReads;
            m_noOfWrites -= a.m_noOfWrites;
            m_readVolume -= a.m_readVolume;
            m_writeVolume -= a.m_writeVolume;
            return *this;
        };
        Stats& operator +=(const Stats& a) {
            m_sec += a.m_sec;
            m_readTime += a.m_readTime;
            m_writeTime += a.m_writeTime;
            m_pioTime += a.m_pioTime;
            m_ioWaitTime += a.m_ioWaitTime;
            m_noOfReads += a.m_noOfReads;
            m_noOfWrites += a.m_noOfWrites;
            m_readVolume += a.m_readVolume;
            m_writeVolume += a.m_writeVolume;
            return *this;
        };
        void start();
        Stats get_stats();

        Timer TT_IO;
        double m_sec//overall time
        , m_readTime//read io time
        , m_writeTime//write io time
        , m_pioTime//read+write io time, second
        , m_ioWaitTime;//io wait time
        unsigned m_noOfReads//number of read io
        , m_noOfWrites//number of write io
        , m_aveBlockSizeR//average size of write block
        , m_aveBlockSizeW;//average size of write block
        unsigned long long int m_readVolume//volume of read bytes
        , m_writeVolume;//volume of write bytes
    private:
        //        friend void get_current_stats(Stats&);//recording current statistics
        friend Stats operator + (const Stats&, const Stats&);
        friend Stats operator - (const Stats&, const Stats&);
        friend std::ostream& operator<<(std::ostream&, const Stats&);
    };
    void Stats::start() {
        stxxl::stats* stxxlStats = stxxl::stats::get_instance();
        m_readTime = stxxlStats->get_read_time();   //get read io time, second
        m_writeTime = stxxlStats->get_write_time(); //get write io time, second
        m_pioTime = stxxlStats->get_pio_time();     //get overall io time, second
        m_ioWaitTime = stxxlStats->get_io_wait_time();//get io wait time, second
        m_noOfReads = stxxlStats->get_reads();      //get the number of read io
        m_noOfWrites = stxxlStats->get_writes();    //get the number of write io
        m_readVolume = stxxlStats->get_read_volume();//get the volume of bytes read from disk
        m_writeVolume = stxxlStats->get_written_volume();//get the volume of bytes read from disk
        TT_IO.start();
    }
    Stats Stats::get_stats() {
        Stats SoFar;
        stxxl::stats* stxxlStats = stxxl::stats::get_instance();
        SoFar.m_readTime = stxxlStats->get_read_time() - m_readTime;//get read io time
        SoFar.m_writeTime = stxxlStats->get_write_time() - m_writeTime;//get write io time
        SoFar.m_pioTime = stxxlStats->get_pio_time() - m_pioTime;//get overall io time
        SoFar.m_ioWaitTime = stxxlStats->get_io_wait_time() - m_ioWaitTime;//get io wait time
        SoFar.m_noOfReads = stxxlStats->get_reads() - m_noOfReads;//get the number of read io
        SoFar.m_noOfWrites = stxxlStats->get_writes() - m_noOfWrites;//get the number of write io
        SoFar.m_readVolume = stxxlStats->get_read_volume() - m_readVolume;//get the volume of bytes read from disk
        SoFar.m_writeVolume = stxxlStats->get_written_volume() - m_writeVolume;//get the volume of bytes read from disk
        SoFar.m_aveBlockSizeR = (SoFar.m_noOfReads ? SoFar.m_readVolume / SoFar.m_noOfReads : 0);
        SoFar.m_aveBlockSizeW = (SoFar.m_noOfWrites ? SoFar.m_writeVolume / SoFar.m_noOfWrites : 0);
        TT_IO.stop();
        SoFar.m_sec = TT_IO.GetRuntime();
        return SoFar;
    }
    Stats operator + (const Stats& a, const Stats& b) {
        Stats r = a;
        r += b;
        return r;
    };
    Stats operator- (const Stats& a, const Stats& b) {
        Stats r = a;
        r -= b;
        return r;
    };
    std::ostream& operator<<(std::ostream& os, const Stats& stats) {
        os << std::left << std::setw(30) << "Total number of reads: " << stats.m_noOfReads << std::endl;
        os << std::left << std::setw(30) << "Total number of writes: " << stats.m_noOfWrites << std::endl;
        os << std::left << std::setw(30) << "Average block size (read): " << stats.m_aveBlockSizeR << " (" << (double)stats.m_aveBlockSizeR / 1024 << " KiB)" << std::endl;
        os << std::left << std::setw(30) << "Average block size (write): " << stats.m_aveBlockSizeW << " (" << (double)stats.m_aveBlockSizeW / 1024 << " KiB)" << std::endl;
        os << std::left << std::setw(30) << "Total read volume: " << stats.m_readVolume << " (" << (double)stats.m_readVolume / 1048576.0 << " MiB)" << std::endl;
        os << std::left << std::setw(30) << "Total written volume: " << stats.m_writeVolume << " (" << (double)stats.m_writeVolume / 1048576.0 << " MiB)" << std::endl;
        os << std::left << std::setw(30) << "Read time: " << stats.m_readTime << " s @ " << ((double)stats.m_readVolume / 1048576.0 / stats.m_readTime) << " MiB/s" << std::endl;
        os << std::left << std::setw(30) << "Write time: " << stats.m_writeTime << " s @ " << ((double)stats.m_writeVolume / 1048576.0 / stats.m_writeTime) << " MiB/s" << std::endl;
        os << std::left << std::setw(30) << "Total I/O wait time: " << stats.m_ioWaitTime << " s" << std::endl;
        os << std::left << std::setw(30) << "Total I/O time: " << stats.m_pioTime << " s" << std::endl;
        os << std::left << std::setw(30) << "Total time: " << stats.m_sec << " seconds." << std::endl;
        return os;
    }
    void get_current_stats(Stats& stats) {
        //        stats.m_sec = stats.TT_IO.GetRuntime();
        stxxl::stats* stxxlStats = stxxl::stats::get_instance();
        stats.m_readTime = stxxlStats->get_read_time();//get read io time
        stats.m_writeTime = stxxlStats->get_write_time();//get write io time
        stats.m_pioTime = stxxlStats->get_pio_time();//get overall io time
        stats.m_ioWaitTime = stxxlStats->get_io_wait_time();//get io wait time
        stats.m_noOfReads = stxxlStats->get_reads();//get the number of read io
        stats.m_noOfWrites = stxxlStats->get_writes();//get the number of write io
        stats.m_readVolume = stxxlStats->get_read_volume();//get the volume of bytes read from disk
        stats.m_writeVolume = stxxlStats->get_written_volume();//get the volume of bytes read from disk
    }
}


#endif //MCSPS_IOXXL_H
