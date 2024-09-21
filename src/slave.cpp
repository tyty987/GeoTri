#include "slave.hpp"

const double Slave::OAOWEIGHT	= 	1.0 / 1;
const double Slave::OANWEIGHT	= 	1.0 / 2;
const double Slave::NANWEIGHT	= 	1.0 / 3;

Slave::Slave(int k, unsigned int seed, MID slaveId): k(k), n(0), globalCount(0), batchGraphCount(0), generator(seed), distribution(0.0, 1.0), slaveId(slaveId) {
	srand(seed+time(NULL));
	samples.reserve(k);
};

void Slave::updateCount(const Edge &iEdge){

	VID src = iEdge.src;
	VID dst = iEdge.dst;

	if(nodeToNeighbors.find(src) == nodeToNeighbors.end() || nodeToNeighbors.find(dst) == nodeToNeighbors.end()) {
		return;
	}

	std::set<VID> &srcMap = nodeToNeighbors[src][0];
	std::set<VID> &dstMap = nodeToNeighbors[dst][0];

	double countSum = 0;
	std::set<VID>::iterator srcIt;
	for (srcIt = srcMap.begin(); srcIt != srcMap.end(); srcIt++) {
		VID neighbor = *srcIt;
		if (dstMap.find(neighbor) != dstMap.end()) {
			double curSampleNum = k >= n ? n : k;
			double prob = (curSampleNum / n * (curSampleNum - 1) / (n - 1));
			double count = 1.0 / prob;
			countSum += count;
			if (nodeToCount.find(neighbor) == nodeToCount.end()) {
				nodeToCount[neighbor] = (float)count;
			} else {
				nodeToCount[neighbor] += (float)count;
			}
		}
	}

	if(countSum > 0) {
		if (nodeToCount.find(src) == nodeToCount.end()) {
			nodeToCount[src] = (float)countSum;
		} else {
			nodeToCount[src] += (float)countSum;
		}

		if (nodeToCount.find(dst) == nodeToCount.end()) {
			nodeToCount[dst] = (float)countSum;
		} else {
			nodeToCount[dst] += (float)countSum;
		}

		globalCount += countSum;
	}


	return;
}

void Slave::weightMergeCount(std::set<VID> &firstSet, std::set<VID> &secondSet, double sumWeight){
	auto firstIt = firstSet.begin();
	auto secondIt = secondSet.begin();
	while(firstIt != firstSet.end() && secondIt != secondSet.end()) {
		if(*firstIt == *secondIt) {
			batchGraphCount += sumWeight;
			firstIt++;
			secondIt++;
			if(sumWeight == OAOWEIGHT) {
				oneNewEdgeTriangle++;
			}else if(sumWeight == NANWEIGHT){
				threeNewEdgesTriangle++;
			} else {
				twoNewEdgesTriangle++; 
			}
		} else if(*firstIt > *secondIt) {
			secondIt++;
		} else {
			firstIt++;
		}
	}
}

void Slave::weightTirangleCount(const Edge &iEdge){
	if(iEdge.srcSlaveId == iEdge.dstSlaveId) {
		weightMergeCount(nodeToNeighbors[iEdge.src][0], nodeToNeighbors[iEdge.dst][0], OAOWEIGHT);
		weightMergeCount(nodeToNeighbors[iEdge.src][1], nodeToNeighbors[iEdge.dst][0], OANWEIGHT);
		weightMergeCount(nodeToNeighbors[iEdge.src][0], nodeToNeighbors[iEdge.dst][1], OANWEIGHT);
		weightMergeCount(nodeToNeighbors[iEdge.src][1], nodeToNeighbors[iEdge.dst][1], NANWEIGHT);
		computeTimes++;
	} else if(iEdge.srcSlaveId == iEdge.computeSlaveId) {
		if(slaveNeighbors.find(iEdge.dst) == slaveNeighbors.end()) {
			return;
		}
		weightMergeCount(nodeToNeighbors[iEdge.src][0], slaveNeighbors[iEdge.dst][0], OAOWEIGHT);
		weightMergeCount(nodeToNeighbors[iEdge.src][1], slaveNeighbors[iEdge.dst][0], OANWEIGHT);
		weightMergeCount(nodeToNeighbors[iEdge.src][0], slaveNeighbors[iEdge.dst][1], OANWEIGHT);
		weightMergeCount(nodeToNeighbors[iEdge.src][1], slaveNeighbors[iEdge.dst][1], NANWEIGHT);
		computeTimes++;
	} else if(iEdge.dstSlaveId == iEdge.computeSlaveId) {
		if(slaveNeighbors.find(iEdge.src) == slaveNeighbors.end()) {
			return;
		}
		weightMergeCount(slaveNeighbors[iEdge.src][0], nodeToNeighbors[iEdge.dst][0], OAOWEIGHT);
		weightMergeCount(slaveNeighbors[iEdge.src][1], nodeToNeighbors[iEdge.dst][0], OANWEIGHT);
		weightMergeCount(slaveNeighbors[iEdge.src][0], nodeToNeighbors[iEdge.dst][1], OANWEIGHT);
		weightMergeCount(slaveNeighbors[iEdge.src][1], nodeToNeighbors[iEdge.dst][1], NANWEIGHT);
		computeTimes++;
	}
}

int Slave::deleteEdge() {
	int index = rand() % k;
	Edge removedEdge = samples[index];
	nodeToNeighbors[removedEdge.src][0].erase(removedEdge.dst);
	nodeToNeighbors[removedEdge.dst][0].erase(removedEdge.src);
	return index;
}

void Slave::tirangleCount(const Edge &iEdge){
	if(iEdge.srcSlaveId == iEdge.dstSlaveId) {
		mergeCount(nodeToNeighbors[iEdge.src][0], nodeToNeighbors[iEdge.dst][0], iEdge.src, iEdge.dst, 0);
		mergeCount(nodeToNeighbors[iEdge.src][1], nodeToNeighbors[iEdge.dst][0], iEdge.src, iEdge.dst, 1);
		mergeCount(nodeToNeighbors[iEdge.src][0], nodeToNeighbors[iEdge.dst][1], iEdge.src, iEdge.dst, 2);
		mergeCount(nodeToNeighbors[iEdge.src][1], nodeToNeighbors[iEdge.dst][1], iEdge.src, iEdge.dst, 3);
		computeTimes++;
	} else if(iEdge.srcSlaveId == iEdge.computeSlaveId) {
		if(slaveNeighbors.find(iEdge.dst) == slaveNeighbors.end()) {
			return;
		}
		mergeCount(nodeToNeighbors[iEdge.src][0], slaveNeighbors[iEdge.dst][0], iEdge.src, iEdge.dst, 0);
		mergeCount(nodeToNeighbors[iEdge.src][1], slaveNeighbors[iEdge.dst][0], iEdge.src, iEdge.dst, 1);
		mergeCount(nodeToNeighbors[iEdge.src][0], slaveNeighbors[iEdge.dst][1], iEdge.src, iEdge.dst, 2);
		mergeCount(nodeToNeighbors[iEdge.src][1], slaveNeighbors[iEdge.dst][1], iEdge.src, iEdge.dst, 3);
		computeTimes++;
	} else if(iEdge.dstSlaveId == iEdge.computeSlaveId) {
		if(slaveNeighbors.find(iEdge.src) == slaveNeighbors.end()) {
			return;
		}
		mergeCount(slaveNeighbors[iEdge.src][0], nodeToNeighbors[iEdge.dst][0], iEdge.src, iEdge.dst, 0);
		mergeCount(slaveNeighbors[iEdge.src][1], nodeToNeighbors[iEdge.dst][0], iEdge.src, iEdge.dst, 1);
		mergeCount(slaveNeighbors[iEdge.src][0], nodeToNeighbors[iEdge.dst][1], iEdge.src, iEdge.dst, 2);
		mergeCount(slaveNeighbors[iEdge.src][1], nodeToNeighbors[iEdge.dst][1], iEdge.src, iEdge.dst, 3);
		computeTimes++;
	}
}

void Slave::mergeCount(std::set<VID> &firstSet, std::set<VID> &secondSet, VID srcId, VID dstId, int flag){
	std::set<VID>::iterator firstIt = firstSet.begin();
	std::set<VID>::iterator secondIt = secondSet.begin();
	if(flag == 1) {
		auto pair2 = upper_bound(secondSet.begin(), secondSet.end(), dstId);
		secondIt = pair2;
	}else if(flag == 2) {
		auto pair1 = upper_bound(firstSet.begin(), firstSet.end(), srcId);
		firstIt = pair1;
	}else if(flag == 3){
		auto pair1 = upper_bound(firstSet.begin(), firstSet.end(), srcId);
		firstIt = pair1;
		auto pair2 = upper_bound(secondSet.begin(), secondSet.end(), dstId);
		secondIt = pair2;
	}

	while(firstIt != firstSet.end() && secondIt != secondSet.end()) {
		if(*firstIt == *secondIt) {
			batchGraphCount += 1;
			firstIt++;
			secondIt++;
		} else if(*firstIt > *secondIt) {
			secondIt++;
		} else {
			firstIt++;
		}
	}
}

void Slave::processBatchEdge(const Edge &iEdge){
	if(iEdge.computeSlaveId == slaveId) {
		Edge edge(iEdge);
		countEdge.push_back(edge);
	}

	VID src = iEdge.src;
	VID dst = iEdge.dst;

	if(src == dst) {
		return;
	}

	if(nodeToNeighbors.find(src) == nodeToNeighbors.end()) {
		std::vector<std::set<VID>> tempVec(2, std::set<VID>());
		nodeToNeighbors[src] = tempVec;
	}
	nodeToNeighbors[src][1].insert(dst);

	if(nodeToNeighbors.find(dst) == nodeToNeighbors.end()) {
		std::vector<std::set<VID>> tempVec(2, std::set<VID>());
		nodeToNeighbors[dst] = tempVec;
	}
	nodeToNeighbors[dst][1].insert(src);

	return;
}

void Slave::processEdge(const Edge &iEdge){

	VID src = iEdge.src;
	VID dst = iEdge.dst;

	if(src == dst) {
		return;
	}

	updateCount(iEdge);

	bool isSampled = false;
	if(n < k) {
		isSampled = true;
	}
	else {

		if(distribution(generator) < k / (1.0+n)) {
			isSampled = true;
		}
	}

	if(isSampled) {
		if(nodeToNeighbors.find(src)==nodeToNeighbors.end()) {
			std::vector<std::set<VID>> tempVec(2, std::set<VID>());
			nodeToNeighbors[src] = tempVec;
		}
		nodeToNeighbors[src][0].insert(dst);

		if(nodeToNeighbors.find(dst)==nodeToNeighbors.end()) {
			std::vector<std::set<VID>> tempVec(2, std::set<VID>());
			nodeToNeighbors[dst] = tempVec;
		}
		nodeToNeighbors[dst][0].insert(src);
	}

	n++;

	return;

}

void Slave::processSlaveNeighbour(Node nodeMaster, Node iNode)
{
	if(slaveNeighbors.find(nodeMaster.index) == slaveNeighbors.end()) {
		std::vector<std::set<VID>> tempVec(2, std::set<VID>());
		slaveNeighbors[nodeMaster.index] = tempVec;
	}

	int newOrOldIndex = iNode.slaveID == -1 ? 0 : 1;
	slaveNeighbors[nodeMaster.index][newOrOldIndex].insert(iNode.index);

	return;
}

double Slave::getGlobalCount()
{
	return globalCount;
}

void Slave::processEdgeWithoutSampling(const Edge &iEdge)
{
	if(iEdge.src == iEdge.dst) {
		return;
	}

	updateCount(iEdge);
}

std::unordered_map<VID, std::vector<std::set<VID>>>* Slave::getSlaveNeighbors(){
	return &slaveNeighbors;
}

double Slave::getBatchGraphCount(){
	return batchGraphCount;
}

std::unordered_map<VID, float> & Slave::getLocalCount()
{
	return nodeToCount;
}

std::unordered_map<VID, std::vector<std::set<VID>>>* Slave::getNodeToNeighbors(){
	return &nodeToNeighbors;
}



std::vector<Edge>* Slave::getCountEdge(){
	return &countEdge;
}

bool Slave::doStore(Edge &iEdge){
	VID src = iEdge.src;
	VID dst = iEdge.dst;
	bool flag = false;
	if(iEdge.srcSlaveId == slaveId) {
		storeNode.insert(src);
		flag = true;
	}
	if(iEdge.dstSlaveId == slaveId) {
		storeNode.insert(dst);
		flag = true;
		
	}
	return flag;
}

#ifdef _TEST_

#endif
