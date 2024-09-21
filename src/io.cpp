#include "io.hpp"

MPI_Datatype 		MPIIO::MPI_TYPE_EDGE;
MPI_Datatype 		MPIIO::MPI_TYPE_ELEMCNT;
MPI_Datatype 		MPIIO::MPI_TYPE_NODE;
const Edge 			MPIIO::END_STREAM(INVALID_VID, INVALID_VID, -1, -1, -1);
const Node 			MPIIO::END_STREAM_NODE(INVALID_VID, -3);

unsigned short		MPIIO::lenBuf;

MPIIO::MPIIO(int &argc, char** &argv)
{
	MPI_Init(&argc, &argv);

	MPI_Comm_size(MPI_COMM_WORLD, &szProc);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	int          lenAttr[Edge::szAttr] = {1, 1, 1, 1, 1};
	MPI_Datatype arrType[Edge::szAttr] = {MPI_UNSIGNED, MPI_UNSIGNED, MPI_SHORT, MPI_SHORT, MPI_SHORT};

	MPI_Aint     offsets[Edge::szAttr];
	offsets[0] = offsetof(Edge, src);
	offsets[1] = offsetof(Edge, dst);
	offsets[2] = offsetof(Edge, srcSlaveId);
	offsets[3] = offsetof(Edge, dstSlaveId);
	offsets[4] = offsetof(Edge, computeSlaveId);
	MPI_Type_create_struct(Edge::szAttr, lenAttr, offsets, arrType, &MPI_TYPE_EDGE);
	MPI_Type_commit(&MPI_TYPE_EDGE);

	arrType[0] = MPI_UNSIGNED;
	arrType[1] = MPI_DOUBLE;
	offsets[0] = offsetof(ElemCount, vid);
	offsets[1] = offsetof(ElemCount, count);
	MPI_Type_create_struct(ElemCount::szAttr, lenAttr, offsets, arrType, &MPI_TYPE_ELEMCNT);
	MPI_Type_commit(&MPI_TYPE_ELEMCNT);


	int          lenAttrNode[Node::szAttrNode] = {1, 1};
	MPI_Datatype arrTypeNode[Node::szAttrNode] = {MPI_UNSIGNED, MPI_SHORT};
	MPI_Aint     offsetsNode[Node::szAttrNode];

	offsetsNode[0] = offsetof(Node, index);
	offsetsNode[1] = offsetof(Node, slaveID);
	MPI_Type_create_struct(Node::szAttrNode, lenAttrNode, offsetsNode, arrTypeNode, &MPI_TYPE_NODE);
	MPI_Type_commit(&MPI_TYPE_NODE);
}

void MPIIO::init(int lenBuf, int slaveNum)
{
    CPUIOTime 			= 0;
	sendCount			= 0;
    commCostGather 		= 0;
	commCostDistribute 	= 0;

    eBuf.clear();
	nBuf.clear();
    MPIIO::lenBuf = lenBuf;
    MPIIO::slaveNum = slaveNum;
	MPIIO::w2wEndSignalCount = slaveNum - 1;
	MPIIO::w2wEndRequest = MPI_REQUEST_NULL;

	if(isMaster()) {
		eBuf.resize(slaveNum);
		nBuf.resize(slaveNum);
		for (int i = 0; i < slaveNum; i++)
		{
			eBuf[i].init(i);
			nBuf[i].init(i, getSlaveId());
		}
	} else {
		eBuf.resize(slaveNum + 1);
		nBuf.resize(slaveNum + 1);
		for (int i = 0; i < slaveNum + 1; i++)
		{
			eBuf[i].init(i);
			nBuf[i].init(i, getSlaveId());
		}
	}
}

void MPIIO::reInitNBuf(){
	nBuf.clear();
	nBuf.resize(slaveNum + 1);
	for (int i = 0; i < slaveNum + 1; i++)
	{
		nBuf[i].reInit(i, getSlaveId());
	}
}

void MPIIO::reInitEBuf(){
	eBuf.clear();
	eBuf.resize(slaveNum);
	for (int i = 0; i < slaveNum; i++)
	{
		eBuf[i].reInit(i);
	}
}

bool MPIIO::isMaster()
{
    return rank == MPI_MASTER;
}

MID MPIIO::getSlaveId()
{
    return (MID)(rank - 1);
}

long MPIIO::getCommCostDistribute()
{
    return commCostDistribute;
}

long MPIIO::getCommCostGather()
{
    return commCostGather;

}

void MPIIO::cleanup()
{
	MPI_Finalize();
}

bool MPIIO::sendEdge(const Edge &iEdge, MID dst)
{
	Edge tmpEdge = iEdge;
	commCostDistribute++;
	eBuf[dst].putNextNode(tmpEdge);

	return true;
}

bool MPIIO::sendNeighbour(MID dst, VID nodeID, std::vector<std::set<VID>>& neighbours, const int tag){
	Node node(nodeID, dst);
	int bufIndex = nBuf[dst].putNextNode(node, tag);
	commCostDistributeNode++;
	for(unsigned int i = 0; i < neighbours.size();i++) {
		int flag = i == 0 ? -1 : -2;
		for(auto it = neighbours[i].begin(); it != neighbours[i].end();it++){
			if(bufIndex == 0) {
				bufIndex = nBuf[dst].putNextNode(node, tag);
				commCostDistributeNode++;
			}
			Node neighbourNode(*it, flag);
			bufIndex = nBuf[dst].putNextNode(neighbourNode, tag);
			commCostDistributeNode++;
		}
	}
	return true;
}

bool MPIIO::sendSlaveToSlaveEdge(VID nodeID, MID src, MID dst){
	Node sendInfoNode(nodeID, dst);
	nBuf[src].putNextNode(sendInfoNode);
	commCostDistributeNode++;
	
	return true;
}

bool MPIIO::bCastEdge(Edge &iEdge)
{
	Edge tmpEdge(iEdge);
	for (int mit = 0; mit < slaveNum; mit++)
	{
		eBuf[mit].putNextNode(tmpEdge);
	}
	commCostDistribute += slaveNum;

	return true;
}

bool MPIIO::IreceiveEdges(Edge *buf, MPI_Request &iReq)
{
	return (MPI_SUCCESS == MPI_Irecv(buf, lenBuf, MPI_TYPE_EDGE, MPI_MASTER, TAG_STREAM, MPI_COMM_WORLD, &iReq));
}

bool MPIIO::IsendEdge(Edge *buf, int mid, MPI_Request &iReq){
	MPI_Isend(buf, lenBuf, MPI_TYPE_EDGE, mid + 1, TAG_STREAM, MPI_COMM_WORLD, &iReq);
	return true;
}

bool MPIIO::IrecvNode(Node *nodeBuf, MPI_Request &iReq, const int tag, const int source){
	return (MPI_SUCCESS == MPI_Irecv(nodeBuf, lenBuf, MPI_TYPE_NODE, source, tag, MPI_COMM_WORLD, &iReq));
};

bool MPIIO::IsendNode(Node *nodeBuf, int mid, MPI_Request &iReq, const int tag){
	MPI_Isend(nodeBuf, lenBuf, MPI_TYPE_NODE, mid + 1, tag, MPI_COMM_WORLD, &iReq);
	return true;
};

bool MPIIO::receiveEdges(Edge &oEdge)
{
	eBuf[0].getNext(oEdge);
	if (oEdge == END_STREAM)
	{
		eBuf[0].cleanup();
	}
	return (oEdge != END_STREAM);
}

bool MPIIO::receiveNeighbour(Node &oNode, const int tag, const int source)
{
	nBuf[getSlaveId()].getNext(oNode, tag, source);
	if (oNode == END_STREAM_NODE)
	{
		nBuf[getSlaveId()].cleanup();
	}
	return (oNode != END_STREAM_NODE);
}

bool MPIIO::sendEndSignal()
{
	Edge signal(END_STREAM);
	for (int mit = 0; mit < slaveNum; mit++)
	{
		eBuf[mit].putNextNode(signal);
		eBuf[mit].flushCacheSend();
	}
	return true;
}

bool MPIIO::sendNeighbourEndSignal()
{
	Node signal(END_STREAM_NODE);
	for (int mit = 0; mit < slaveNum; mit++)
	{
		nBuf[mit].putNextNode(signal);
		nBuf[mit].flushCacheSend();
	}
	return true;
}

bool MPIIO::sendW2WNeighbourEnd()
{
	Node signal(END_STREAM_NODE);
	for (int mit = 0; mit < slaveNum; mit++)
	{
		if(mit == getSlaveId()) {
			continue;
		}
		nBuf[mit].putNextNode(signal, TAG_W2W_NEIGHBOUR);
		nBuf[mit].slaveFlushSend();
	}
	return true;
}

bool MPIIO::sendCount(double gCount, unordered_map<VID, float> &lCount)
{
    clock_t begin = clock();
	MPI_Reduce(&gCount, nullptr, 1, MPI_DOUBLE, MPI_SUM, MPI_MASTER, MPI_COMM_WORLD);
    VID maxVId;
    MPI_Bcast(&maxVId, 1, MPI_UNSIGNED, MPI_MASTER, MPI_COMM_WORLD);
    CPUIOTime += double(clock() - begin);


    float* lCountArr = new float[maxVId+1];
    std::fill_n(lCountArr, maxVId+1, 0.0);
    unordered_map<VID, float>::const_iterator it;
    for (it = lCount.begin(); it != lCount.end(); it++) {
        lCountArr[it->first] = it->second;
    }

    begin = clock();
    MPI_Reduce(lCountArr, nullptr, maxVId+1, MPI_FLOAT, MPI_SUM, MPI_MASTER, MPI_COMM_WORLD);
    CPUIOTime += double(clock() - begin);

    delete lCountArr;
    return true;
}

bool MPIIO::receiveCount(VID maxVId, double &gCount, std::vector<float> &lCount)
{
	double empty = 0;
    clock_t begin = clock();
	MPI_Reduce(&empty, &gCount, 1, MPI_DOUBLE, MPI_SUM, MPI_MASTER, MPI_COMM_WORLD);
    MPI_Bcast(&maxVId, 1, MPI_UNSIGNED, MPI_MASTER, MPI_COMM_WORLD);
    CPUIOTime += double(clock() - begin);

    float* lCountArr = new float[maxVId+1];
    std::fill_n(lCountArr, maxVId+1, 0.0);

    begin = clock();
    MPI_Reduce(MPI_IN_PLACE, lCountArr, maxVId+1, MPI_FLOAT, MPI_SUM, MPI_MASTER, MPI_COMM_WORLD);
    CPUIOTime += double(clock() - begin);

    commCostGather = (maxVId + 1) * (getSzProc()-1);

    lCount.insert(lCount.end(), &lCountArr[0], &lCountArr[maxVId]);
    delete lCountArr;

	return true;
}

double MPIIO::getCPUIOTime()
{
	if (rank == MPI_MASTER)
	{
		double totalIOCPUTime = CPUIOTime;
		for (int i = 0; i < slaveNum; i++)
		{
			totalIOCPUTime += eBuf[i].CPUIOTime;
		}
		return totalIOCPUTime;
	}
	else
	{
		return eBuf[0].CPUIOTime;
	}
}

long MPIIO::getSendCount()
{
	if (rank == MPI_MASTER)
	{
		long totalSendCount = 0;
		for (int i = 0; i < slaveNum; i++)
		{
			totalSendCount += eBuf[i].sendCount;
			totalSendCount += nBuf[i].sendCount;
		}
		return totalSendCount;
	}
	else
	{
		long totalSendCount = 0;
		for (int i = 0; i < slaveNum + 1; i++)
		{
			totalSendCount += eBuf[i].sendCount;
			totalSendCount += nBuf[i].sendCount;
		}
		return totalSendCount;
	}
}

bool MPIIO::sendTime(double compTime)
{
    MPI_Reduce(&compTime, nullptr, 1, MPI_DOUBLE, MPI_MAX, MPI_MASTER, MPI_COMM_WORLD);
    MPI_Reduce(&compTime, nullptr, 1, MPI_DOUBLE, MPI_SUM, MPI_MASTER, MPI_COMM_WORLD);
    return true;
}

bool MPIIO::receiveTime(double &compTimeMax, double &compTimeSum)
{
	double empty = 0;
	MPI_Reduce(&empty, &compTimeMax, 1, MPI_DOUBLE, MPI_MAX, MPI_MASTER, MPI_COMM_WORLD);
    MPI_Reduce(&empty, &compTimeSum, 1, MPI_DOUBLE, MPI_SUM, MPI_MASTER, MPI_COMM_WORLD);
	return true;

}

bool MPIIO::sendRecvNeighbourNum(int neighbourNum)
{
    MPI_Reduce(&neighbourNum, nullptr, 1, MPI_INT, MPI_SUM, MPI_MASTER, MPI_COMM_WORLD);
    return true;
}

bool MPIIO::recvRecvNeighbourNum(int& neighbourNum)
{
	double empty = 0;
    MPI_Reduce(&empty, &neighbourNum, 1, MPI_INT, MPI_SUM, MPI_MASTER, MPI_COMM_WORLD);
	return true;
}

bool MPIIO::sendBatchGraphCount(double count){
	MPI_Reduce(&count, nullptr, 1, MPI_DOUBLE, MPI_SUM, MPI_MASTER, MPI_COMM_WORLD);
    return true;
}
bool MPIIO::recvBatchGraphCount(double& count){
	double empty = 0;
    MPI_Reduce(&empty, &count, 1, MPI_DOUBLE, MPI_SUM, MPI_MASTER, MPI_COMM_WORLD);
	return true;
}

#ifdef _TEST_

#endif