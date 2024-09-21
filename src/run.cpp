#include "run.hpp"
void run_experiment (const char* input, const char* outPath, MPIIO &mpiIo, int slaveNum, Method_t method, int memSize, int repeat, int bufLen, double tolerance)
{

    int seed = 0;

    struct timeval diff, startTV, endTV;

	if (mpiIo.isMaster())
	{
		struct stat sb;
		if (stat(outPath, &sb) == 0)
		{
			if (S_ISDIR(sb.st_mode))
				;
			else if (S_ISREG(sb.st_mode))
				;
			else
				;
		} 
		else 
		{
			mkdir(outPath, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
		}
	}

    for(int i =0 ; i < repeat; i++) {

        if (mpiIo.isMaster()) {

            gettimeofday(&startTV, NULL);

            std::vector<float> nodeToCount;

            double srcCompCost = 0;
            double slaveCompCostMax = 0;
            double slaveCompCostSum = 0;

            double globalCount = run_experiment_mpi(input, mpiIo, slaveNum, method, memSize, bufLen, tolerance, seed + repeat * slaveNum * i, nodeToCount, srcCompCost, slaveCompCostMax, slaveCompCostSum);

            gettimeofday(&endTV, NULL);

            timersub(&endTV, &startTV, &diff);

            double elapsedTime = diff.tv_sec * 1000 + diff.tv_usec / 1000 ;

            cout<<"globalCount:"<<std::setprecision(std::numeric_limits<double>::max_digits10)<<globalCount<<endl;
            cout<<"elapsedTime:"<<std::setprecision(std::numeric_limits<double>::max_digits10)<<elapsedTime / 1000<<endl;

        } else {

            double srcCompCost = 0;
            double slaveCompCostMax = 0;
            double slaveCompCostSum = 0;
            std::vector<float> nodeToCount;
            run_experiment_mpi(input, mpiIo, slaveNum, method, memSize, bufLen, tolerance, seed + repeat * slaveNum * i, nodeToCount, srcCompCost, slaveCompCostMax, slaveCompCostSum);
        }
    }
}

double run_experiment_mpi(const char* filename, MPIIO &mpiIo, int slaveNum, const Method_t method, int memSize, int lenBuf, double tolerance, unsigned int seed, std::vector<float> & oLocalCount, double &srcCompCost, double &slaveCompCostMax, double &slaveCompCostSum)
{

    clock_t begin = clock();

    mpiIo.init(lenBuf, slaveNum);

    string filenameStr(filename);
    string filenameFront    =   getFilename(filenameStr, "_true_front_91");
    string filenameBack     =   getFilename(filename, "_true_back_91");

    if (mpiIo.isMaster())
    {
        Source source(slaveNum, method, tolerance);
        Edge edge;
        MID dst1(0);
        MID dst2(0);
        bool isBroadCast;

        int count = 0;
        EdgeParser parserFront(filenameFront.c_str());
        int sendCount = 0;
        while (parserFront.getEdge(edge))
        {
            if(edge.src != edge.dst) {
                isBroadCast = source.processEdge(edge, dst1, dst2);

                if (!isBroadCast)
                {
                    mpiIo.sendEdge(edge, edge.srcSlaveId);
                } else
                {
                    mpiIo.bCastEdge(edge);
                }
            }
        }
        parserFront.close();
        mpiIo.sendEndSignal();

        clock_t geoTriBegin = clock();

        EdgeParser parserBack(filenameBack.c_str());
        while (parserBack.getEdge(edge))
        {
            source.processBatchGraphEdge(edge);
        }
        parserBack.close();

        sendCount += mpiIo.getSendCount();
        mpiIo.reInitEBuf();

        source.partionNewNode();

        std::vector<Edge>   w2wInfoVec;
        source.calculateEdgeSendDirection(w2wInfoVec);
        std::vector<Edge>* sendEdgeVec    = source.getSendEdgeVec();
        for(auto &edge : *sendEdgeVec)
        {
            if(edge.srcSlaveId == edge.dstSlaveId) {
                mpiIo.sendEdge(edge, edge.srcSlaveId);
            } else {
                mpiIo.sendEdge(edge, edge.srcSlaveId);
                mpiIo.sendEdge(edge, edge.dstSlaveId);
            }
           
        }
        mpiIo.sendEndSignal();

        MID* nodeToSlave = source.getNodeToSlave();
        for(Edge &edge : w2wInfoVec) {
            if(edge.srcSlaveId == edge.computeSlaveId) {
                mpiIo.sendSlaveToSlaveEdge(edge.dst, nodeToSlave[edge.dst], nodeToSlave[edge.src]);
            } else {
                mpiIo.sendSlaveToSlaveEdge(edge.src, nodeToSlave[edge.src], nodeToSlave[edge.dst]);
            }
        }
        mpiIo.sendNeighbourEndSignal();

        int receiveNeighbourNum;
        mpiIo.recvRecvNeighbourNum(receiveNeighbourNum);

        double batchGraphCount = 0;
        mpiIo.recvBatchGraphCount(batchGraphCount);

        double globalCount = 0;

        mpiIo.receiveCount(source.getMaxVId(), globalCount, oLocalCount);

        mpiIo.receiveTime(slaveCompCostMax, slaveCompCostSum);

        if(method == Method_t::NAIVE)
        {
            globalCount = globalCount / slaveNum;
            for(auto it = oLocalCount.begin(); it != oLocalCount.end(); ++it)
            {
                *it  = *it / slaveNum;
            }
        }

        return globalCount;
    }
    else
    {
        int sendCount = 0;
        srand(mpiIo.getSlaveId() + (unsigned)time(NULL));
        Slave  slave(memSize, seed + mpiIo.getSlaveId(), mpiIo.getSlaveId());
        Edge edge;
        while(mpiIo.receiveEdges(edge))
        {
            if(slave.doStore(edge) || method == Method_t::NAIVE)
            {
                slave.processEdge(edge);
            } else
            {
                slave.processEdgeWithoutSampling(edge);
            }
        }

        mpiIo.reInitEBuf();

        Edge newEdge;
        while(mpiIo.receiveEdges(newEdge))
        {
            slave.doStore(newEdge);
            slave.processBatchEdge(newEdge);
        }

        Node node;
        Node nodeMaster;
        Node sendInfoNode;
        std::unordered_map<VID, std::vector<std::set<VID>>>* slaveNodeToNeighbours = slave.getNodeToNeighbors();
        int receiveNeighbourNum = 0;
        vector<Node> infoNodeVec;

        while(mpiIo.receiveNeighbour(sendInfoNode))
        {
            infoNodeVec.push_back(sendInfoNode);
        }


        for(int i = 0;i < mpiIo.slaveNum;i++) {
            if(mpiIo.getSlaveId() == i) {
                for(Node &node : infoNodeVec) {
                    mpiIo.sendNeighbour(node.slaveID, node.index, (*slaveNodeToNeighbours)[node.index], mpiIo.TAG_W2W_NEIGHBOUR);
                }
                mpiIo.sendW2WNeighbourEnd();
            } else {
                while(mpiIo.receiveNeighbour(node, mpiIo.TAG_W2W_NEIGHBOUR, i + 1)) {
                    if(node.slaveID >= 0) {
                        nodeMaster = node;
                    } else {
                        slave.processSlaveNeighbour(nodeMaster, node);
                        receiveNeighbourNum++;
                    }
                }
            }
            sendCount += mpiIo.getSendCount();
            mpiIo.reInitNBuf();
        }

        std::vector<Edge>* countEdge = slave.getCountEdge();
        for(auto &oEdge : *countEdge) {
            slave.tirangleCount(oEdge);
        }

        mpiIo.sendRecvNeighbourNum(receiveNeighbourNum);

        mpiIo.sendBatchGraphCount(slave.batchGraphCount);

        mpiIo.sendCount(slave.getGlobalCount(), slave.getLocalCount());

        double slaveCompCost = (double(clock() - begin) - mpiIo.getCPUIOTime()) / CLOCKS_PER_SEC;

        mpiIo.sendTime(slaveCompCost);


        return 0;
    }
}

void print_count(const char* outPath, double globalCount, const std::vector<float> &localCount, int id)
{

	std::ostringstream gCountFileName;
	gCountFileName << outPath << "/global" << id << ".txt";
	std::fstream	gfp;
	gfp.open(gCountFileName.str(), std::fstream::out | std::fstream::trunc);
	gfp << std::setprecision(std::numeric_limits<double>::max_digits10) <<  globalCount << endl;
	gfp.close();

	std::ostringstream lCountFileName;
	lCountFileName << outPath << "/local" << id << ".txt";
	std::fstream	lfp;
	lfp.open(lCountFileName.str(), std::fstream::out | std::fstream::trunc);

	for (int nid = 0; nid < localCount.size(); nid++)
	{
		lfp << nid << "\t"  << localCount[nid] << endl;
	}
	lfp.close();
}