#ifndef _SLAVE_HPP_
#define _SLAVE_HPP_

#include "source.hpp"
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <algorithm>
#include <stdlib.h>     /* srand, rand */
#include <random>
#include <vector>
#include "time.h"

class Slave {
_PRIVATE:

	const int k;
	std::vector<Edge> samples;
	long n; 
	double globalCount;
	VID batchGraphCount;
	MID slaveId;
	std::unordered_set<VID> storeNode;
	std::vector<Edge> countEdge;
	std::unordered_map<VID, float> nodeToCount;
	std::unordered_map<VID, std::vector<std::set<VID>>> nodeToNeighbors;
	std::unordered_map<VID, std::vector<std::set<VID>>> slaveNeighbors;

	static const double OAOWEIGHT;
	static const double OANWEIGHT;
	static const double NANWEIGHT;


	std::default_random_engine generator;
	std::uniform_real_distribution<double> distribution;

public:

	Slave(int k, unsigned int seed, MID slaveId);
	
	int oneNewEdgeTriangle = 0;
	int twoNewEdgesTriangle = 0;
	int threeNewEdgesTriangle = 0;

	void updateCount(const Edge &iEdge);

	int deleteEdge();

	void processEdge(const Edge &iEdge);
	void processEdgeWithoutSampling(const Edge &iEdge);
	void processBatchEdge(const Edge &iEdge);
	void processSlaveNeighbour(Node nodeMaster, Node iNode);
	void tirangleCount(const Edge &iEdge);
	void mergeCount(std::set<VID> &firstSet, std::set<VID> &secondSet, VID srcId, VID dstId, int flag);
	void weightTirangleCount(const Edge &iEdge);
	void weightMergeCount(std::set<VID> &firstSet, std::set<VID> &secondSet, double sumWeight);
	std::unordered_map<VID, std::vector<std::set<VID>>>* getNodeToNeighbors();
	std::unordered_map<VID, std::vector<std::set<VID>>>* getSlaveNeighbors();
	std::vector<Edge>* getCountEdge();
	double getGlobalCount();
	double getBatchGraphCount();
	std::unordered_map<VID, float> & getLocalCount();
	bool doStore(Edge &iEdge);
	int computeTimes = 0;
};



#endif
