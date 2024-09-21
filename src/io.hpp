#ifndef _IO_HPP_
#define _IO_HPP_

#include "struct.hpp"
#include "split.cpp"
#include <mpi.h>
#include <cstdlib>
#include <vector>
#include <set>
#include <unordered_map>
#include <iostream>
#include <stddef.h>
#include <ctime>
#include <unistd.h>
#include <sys/time.h>

class MPIIO {
_PRIVATE:
	class EdgeContainer {
	_PRIVATE:
		MID							mid;
		unsigned int				bit;
		unsigned int				qit;
		bool						isEmpty;
		double						CPUIOTime;
		long						sendCount;

	public:
		Edge						*buf[2];
		MPI_Request 				req[2];

		EdgeContainer(): bit(0), qit(0), isEmpty(true), CPUIOTime(0), sendCount(0),
				buf{nullptr, nullptr}, req{MPI_REQUEST_NULL, MPI_REQUEST_NULL}{}

		~EdgeContainer()
		{
			if (buf[0] != nullptr)
			{
				delete[] buf[0];
				buf[0] = nullptr;
			}

			if (buf[1] != nullptr)
			{
				delete[] buf[1];
				buf[1] = nullptr;
			}

			waitIOCompletion(req[0]);
			waitIOCompletion(req[1]);
		}

		inline bool init(MID iMid)
		{
			mid = iMid;
			buf[0] = new Edge[lenBuf];
			buf[1] = new Edge[lenBuf];
			return true;
		}

		inline bool reInit(MID iMid)
		{
			mid = iMid;
			if (buf[0] != nullptr)
			{
				delete[] buf[0];
			}
			buf[0] = new Edge[lenBuf];

			if (buf[1] != nullptr)
			{
				delete[] buf[1];
			}
			buf[1] = new Edge[lenBuf];

			waitIOCompletion(req[0]);
			waitIOCompletion(req[1]);
			return true;
		}

		inline void getNext(Edge &oEdge)
		{
			if (isEmpty)
			{
				clock_t begin = clock();

				MPIIO::IreceiveEdges(buf[qit], req[qit]);
				waitIOCompletion(req[qit]);
				MPIIO::IreceiveEdges(buf[(qit+1)%2], req[(qit+1)%2]);

				CPUIOTime += double(clock() - begin);

				isEmpty = false;
			}
			else if (bit == lenBuf)
			{
				clock_t begin = clock();

				bit = 0;
				MPIIO::IreceiveEdges(buf[qit], req[qit]);
				qit = (qit + 1) % 2;
				waitIOCompletion(req[qit]);

				CPUIOTime += double(clock() - begin);
			}

			oEdge = buf[qit][bit++];

			return;
		}

		void cleanup(){
			int flag(0);
			MPI_Status st;
			for (int i = 0; i < 2; i++)
			{
				if (req[i] != MPI_REQUEST_NULL)
				{
					MPI_Test(&req[i], &flag, &st);
					if (flag == false)
					{
						MPI_Cancel(&req[i]);
					}
					req[i] = MPI_REQUEST_NULL;
				}
			}
		}

		inline bool putNextNode(const Edge &iEdge)
		{

			buf[qit][bit++] = iEdge;
			if (bit == lenBuf)
			{
				clock_t begin = clock();

				bit = 0;
				MPIIO::IsendEdge(buf[qit], mid, req[qit]);
				qit = (qit + 1) % 2;
				waitIOCompletion(req[qit]);

				CPUIOTime += double(clock() - begin);
				sendCount++;
			}

			return true;
		}

		void flushCacheSend()
		{

			clock_t begin = clock();

			if (bit != 0)
			{
				MPIIO::IsendEdge(buf[qit], mid, req[qit]);
			}

			for (int i = 0; i < 2; i++)
			{
				waitIOCompletion(req[i]);
			}

			CPUIOTime += double(clock() - begin);
		}
	};


	class NodeContainer {
	_PRIVATE:
		MID							mid;
		MID							slaveId;
		unsigned int				bit;
		unsigned int				qit;
		bool						isEmpty;
		double						CPUIOTime;
		long						sendCount;

	public:
		Node						*nodeBuf[2];
		MPI_Request 				nodeReq[2];

		NodeContainer(): bit(0), qit(0), isEmpty(true), CPUIOTime(0), sendCount(0),
						 nodeBuf{nullptr, nullptr}, nodeReq{MPI_REQUEST_NULL, MPI_REQUEST_NULL}{}

		~NodeContainer()
		{
			if (nodeBuf[0] != nullptr)
			{
				delete[] nodeBuf[0];
				nodeBuf[0] = nullptr;
			}

			if (nodeBuf[1] != nullptr)
			{
				delete[] nodeBuf[1];
				nodeBuf[1] = nullptr;
			}

			waitIOCompletion(nodeReq[0]);
			waitIOCompletion(nodeReq[1]);
		}

		inline bool init(MID iMid, MID iSlaveId)
		{
			mid = iMid;
			slaveId = iSlaveId;
			nodeBuf[0] = new Node[lenBuf];
			nodeBuf[1] = new Node[lenBuf];
			return true;
		}

		inline bool reInit(MID iMid, MID iSlaveId)
		{
			mid = iMid;
			slaveId = iSlaveId;
			if (nodeBuf[0] != nullptr)
			{
				delete[] nodeBuf[0];
			}
			nodeBuf[0] = new Node[lenBuf];

			if (nodeBuf[1] != nullptr)
			{
				delete[] nodeBuf[1];
			}
			nodeBuf[1] = new Node[lenBuf];

			waitIOCompletion(nodeReq[0]);
			waitIOCompletion(nodeReq[1]);
			return true;
		}

		inline void getNext(Node &oNode, const int tag, const int source)
		{
			if (isEmpty)
			{
				clock_t begin = clock();

				MPIIO::IrecvNode(nodeBuf[qit], nodeReq[qit], tag, source);
				waitIOCompletion(nodeReq[qit]);
				MPIIO::IrecvNode(nodeBuf[(qit+1)%2], nodeReq[(qit+1)%2], tag, source);

				CPUIOTime += double(clock() - begin);

				isEmpty = false;
			}
			else if (bit == lenBuf)
			{
				clock_t begin = clock();

				bit = 0;
				MPIIO::IrecvNode(nodeBuf[qit], nodeReq[qit], tag, source);
				qit = (qit + 1) % 2;
				waitIOCompletion(nodeReq[qit]);

				CPUIOTime += double(clock() - begin);
			}

			oNode = nodeBuf[qit][bit++];

			return;
		}

		void cleanup(){
			int flag(0);
			MPI_Status st;
			for (int i = 0; i < 2; i++)
			{
				if (nodeReq[i] != MPI_REQUEST_NULL)
				{
					MPI_Test(&nodeReq[i], &flag, &st);
					if (flag == false)
					{
						MPI_Cancel(&nodeReq[i]);
					}
					nodeReq[i] = MPI_REQUEST_NULL;
				}
			}
		}

		inline unsigned int putNextNode(const Node &iNode)
		{
			nodeBuf[qit][bit++] = iNode;
			if (bit == lenBuf)
			{
				clock_t begin = clock();

				bit = 0;
				MPIIO::IsendNode(nodeBuf[qit], mid, nodeReq[qit]);
				qit = (qit + 1) % 2;
				waitIOCompletion(nodeReq[qit]);

				CPUIOTime += double(clock() - begin);
				sendCount++;
			}

			return bit;
		}

		inline unsigned int putNextNode(const Node &iNode, const int tag)
		{
			nodeBuf[qit][bit++] = iNode;
			if (bit == lenBuf)
			{
				clock_t begin = clock();

				bit = 0;
				MPIIO::IsendNode(nodeBuf[qit], mid, nodeReq[qit], tag);
				qit = (qit + 1) % 2;
				waitIOCompletion(nodeReq[qit]);

				CPUIOTime += double(clock() - begin);
				sendCount++;
			}

			return bit;
		}

		void flushCacheSend()
		{

			clock_t begin = clock();

			if (bit != 0)
			{
				MPIIO::IsendNode(nodeBuf[qit], mid, nodeReq[qit]);
			}

			for (int i = 0; i < 2; i++)
			{
				waitIOCompletion(nodeReq[i]);
			}

			CPUIOTime += double(clock() - begin);
		}

		void slaveFlushSend()
		{

			clock_t begin = clock();

			if (bit != 0)
			{
				MPIIO::IsendNode(nodeBuf[qit], mid, nodeReq[qit], TAG_W2W_NEIGHBOUR);
			}

			for (int i = 0; i < 2; i++)
			{
				waitIOCompletion(nodeReq[i]);
			}

			CPUIOTime += double(clock() - begin);
		}
	};

	static const int 		TAG_STREAM 				= 0;
	static const int 		TAG_RET					= 1;
	static const int 		TAG_W2W_NEIGHBOUR		= 2;
	static const int 		TAG_END_W2W_NEIGHBOUR	= 3;
	static MPI_Datatype 	MPI_TYPE_EDGE;
	static MPI_Datatype 	MPI_TYPE_ELEMCNT;
	static MPI_Datatype 	MPI_TYPE_NODE;
	static const Edge 		END_STREAM;
	static const Node 		END_STREAM_NODE;
	static unsigned short	lenBuf;

	inline static void waitIOCompletion(MPI_Request &iReq)
	{
		MPI_Status 	st;
		(MPI_REQUEST_NULL != iReq) && MPI_Wait(&iReq, &st);
		iReq = MPI_REQUEST_NULL;
		return;
	}

	static bool IreceiveEdges(Edge *buf, MPI_Request &iReq);
	static bool IsendEdge(Edge *buf, int dst, MPI_Request &iReq);

	static bool IrecvNode(Node *nodeBuf, MPI_Request &iReq, const int tag = TAG_STREAM, const int source = MPI_MASTER);
	static bool IsendNode(Node *nodeBuf, int dst, MPI_Request &iReq, const int tag = TAG_STREAM);

	MPI_Request					req;
	int 						rank;
	int 						szProc;
	int 						slaveNum;
	int 						w2wEndSignalCount;
	vector<EdgeContainer>		eBuf;
	vector<NodeContainer>		nBuf;
	long						commCostDistribute;
	long						commCostDistributeNode;
    long						commCostGather;
	long						CPUIOTime;
	long						sendCount;
	MPI_Request					w2wEndRequest;
	MPI_Status 					endP2PMS;


public:

	MPIIO(int &argc, char** &argv);
	~MPIIO(){}

	void init(int buffersize, int slaveNum);
	void cleanup();

	void reInitNBuf();
	void reInitEBuf();

	bool isMaster();

	MID getSlaveId();

	long getCommCostDistribute();
	long getCommCostGather();

	bool sendEdge(const Edge &iEdge, MID dst);
	bool bCastEdge(Edge& iEdge);
	bool receiveEdges(Edge &oEdge);

	bool sendNeighbour(MID dst, VID nodeID, std::vector<std::set<VID>>& neighbours, const int tag = TAG_STREAM);
	bool receiveNeighbour(Node &oNode, const int tag = TAG_STREAM, const int source = MPI_MASTER);

	bool sendSlaveToSlaveEdge(VID nodeID, MID src, MID dst);

	bool sendCount(double gCount, unordered_map<VID, float> &lCount);
    bool receiveCount(VID maxVId, double &gCount, std::vector<float> &lCount);

	double getCPUIOTime();
	long getSendCount();

	bool sendTime(double compTime);
	bool receiveTime(double &compTimeMax, double &compTimeSum);

	bool sendRecvNeighbourNum(int neighbourNum);
	bool recvRecvNeighbourNum(int& neighbourNum);

	bool sendBatchGraphCount(double count);
	bool recvBatchGraphCount(double& count);

	bool sendEndSignal();
	bool sendNeighbourEndSignal();
	bool recvW2WEndSignal();
	bool sendW2WNeighbourEnd();

	int getRank(){ return rank;}
	int getSzProc(){ return szProc;}

};

#endif
