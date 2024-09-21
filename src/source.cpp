#include "source.hpp"

Source::Source(int slaveNum, Method_t method, double tolerance):
        method(method), slaveNum(slaveNum), tolerancePlusOne(tolerance+1.0), tolerancePlusOnePointTwo(tolerancePlusOne + 0.4), slaveToLoad(nullptr), minLoad(0), threshold(0.0), minLoadSlave(0),
        nodeToSlave(nullptr), nodeNeighbourNum(nullptr), maxVId(128), missingMId(-1)
{

    if(method == Method_t::OPT)
    {
        capacity = maxVId;
        slaveToLoad = new long[slaveNum];
        for (MID i = 0; i < slaveNum; i++)
        {
            slaveToLoad[i] = 0;
        }

        nodeToSlave = new MID[maxVId];
        nodeNeighbourNum = new VID[maxVId];
        for (VID i = 0; i < maxVId; i++)
        {
            nodeToSlave[i] = missingMId;
            nodeNeighbourNum[i] = 0;
        }
    }

}

Source::~Source()
{
    if(slaveToLoad != nullptr) {
        delete[] slaveToLoad;
    }

    if(nodeToSlave != nullptr) {
        delete[] nodeToSlave;
    }

    if(nodeNeighbourNum != nullptr) {
        delete[] nodeNeighbourNum;
    }

}

VID Source::getMaxVId()
{
    return maxVId;
}


/**
 *
 * @param iEdge input edge
 * @param oDstMID1 destination machine 1
 * @param oDstMID2 destination machine 2
 * @return whether to broadcast
 */
bool Source::processEdge(Edge &iEdge, MID &oDstMID1, MID &oDstMID2)
{
    VID src = iEdge.src;
    VID dst = iEdge.dst;
    
    if (src > maxVId) {
        maxVId = src;
    }
    if (dst > maxVId) {
        maxVId = dst;
    }
	if(method == Method_t::NAIVE)
    {
        return true;
	}
	else if(method == Method_t::HASH)
    {
        oDstMID1 = src % slaveNum;
        oDstMID2 = dst % slaveNum;

        return oDstMID1 != oDstMID2;
	}
	else 
	{
        if(maxVId >= capacity) {

            VID newCapacity = 2 * maxVId;
            MID * newNodeToSlave = new MID[newCapacity];
            memcpy(newNodeToSlave, nodeToSlave, capacity * sizeof(MID));
            memset(newNodeToSlave + capacity, missingMId, (newCapacity - capacity) * sizeof(MID));
            delete[] nodeToSlave;
            nodeToSlave = newNodeToSlave;

            VID * newNodeNeighbourNum = new VID[newCapacity];
            memcpy(newNodeNeighbourNum, nodeNeighbourNum, capacity * sizeof(VID));
            memset(newNodeNeighbourNum + capacity, 0, (newCapacity - capacity) * sizeof(VID));
            delete[] nodeNeighbourNum;
            nodeNeighbourNum = newNodeNeighbourNum;

            capacity = newCapacity;
        }

        if(nodeToSlave[src] == missingMId)
        {
            if(nodeToSlave[dst] == missingMId)
            {
                oDstMID1 = minLoadSlave;
                oDstMID2 = minLoadSlave;
                nodeToSlave[src] = oDstMID1;
                nodeToSlave[dst] = oDstMID2;
            }
            else 
            {
                oDstMID2 = nodeToSlave[dst];
                if(slaveToLoad[oDstMID2] <= threshold)
                {
                    oDstMID1 = oDstMID2;
                }
                else
                {
                    oDstMID1 = minLoadSlave;
                }
                nodeToSlave[src] = oDstMID1;
            }
        }
        else
        {
            if(nodeToSlave[dst] == missingMId) 
            {
                oDstMID1 = nodeToSlave[src];
                if(slaveToLoad[oDstMID1] <= threshold)
                {
                    oDstMID2 = oDstMID1;
                }
                else
                {
                    oDstMID2 = minLoadSlave;
                }
                nodeToSlave[dst] = oDstMID2;
            }
            else 
            {
                oDstMID1 = nodeToSlave[src];
                oDstMID2 = nodeToSlave[dst];
            }
        }

        if(oDstMID1 == oDstMID2) {
            slaveToLoad[oDstMID1] += 1;
        }
        else {
            slaveToLoad[oDstMID1] += 1;
            slaveToLoad[oDstMID2] += 1;
        }

       
        if(oDstMID1 == minLoadSlave || oDstMID2 == minLoadSlave)
        {
            minLoad = slaveToLoad[minLoadSlave];
            for(MID i = 0; i < slaveNum; i++) {
                if (slaveToLoad[i] < minLoad) {
                    minLoad = slaveToLoad[i];
                    minLoadSlave = i;
                }
            }
            threshold = minLoad * tolerancePlusOnePointTwo;
        }

        iEdge.srcSlaveId = nodeToSlave[src];
        iEdge.dstSlaveId = nodeToSlave[dst];

        nodeNeighbourNum[src]++;
        nodeNeighbourNum[dst]++;

        return oDstMID1 != oDstMID2;
	}
}

bool Source::processBatchGraphEdge(Edge &iEdge)
{
    
    VID src = iEdge.src;
    VID dst = iEdge.dst;

    if (src > maxVId) {
        maxVId = src;
    }
    if (dst > maxVId) {
        maxVId = dst;
    }
    if(maxVId >= capacity) {

        VID newCapacity = 2 * maxVId;
        
        MID * newNodeToSlave = new MID[newCapacity];
        memcpy(newNodeToSlave, nodeToSlave, capacity * sizeof(MID));
        memset(newNodeToSlave + capacity, missingMId, (newCapacity - capacity) * sizeof(MID));
        delete[] nodeToSlave;
        nodeToSlave = newNodeToSlave;

        VID * newNodeNeighbourNum = new VID[newCapacity];
        memcpy(newNodeNeighbourNum, nodeNeighbourNum, capacity * sizeof(VID));
        memset(newNodeNeighbourNum + capacity, 0, (newCapacity - capacity) * sizeof(VID));
        delete[] nodeNeighbourNum;
        nodeNeighbourNum = newNodeNeighbourNum;

        capacity = newCapacity;
    }

    batchGraph[src].insert(dst);
    batchGraph[dst].insert(src);

    
    if(nodeToSlave[src] != missingMId && nodeToSlave[dst] != missingMId) {
        if(nodeToSlave[src] == nodeToSlave[dst]) {
            slaveToLoad[nodeToSlave[src]] += 1;
        } else {
            slaveToLoad[nodeToSlave[src]] += 1;
            slaveToLoad[nodeToSlave[dst]] += 1;
        }

    } else if(nodeToSlave[src] == missingMId && nodeToSlave[dst] != missingMId) {
        newNode.insert(src);
        slaveToLoad[nodeToSlave[dst]] += 1;
    } else if(nodeToSlave[src] != missingMId && nodeToSlave[dst] == missingMId) {
        newNode.insert(dst);
        slaveToLoad[nodeToSlave[src]] += 1;
    } else {
        newNode.insert(src);
        newNode.insert(dst);
    }

    iEdge.srcSlaveId = nodeToSlave[src];
    iEdge.dstSlaveId = nodeToSlave[dst];

    nodeNeighbourNum[src]++;
    nodeNeighbourNum[dst]++;

    Edge e(iEdge);
    batchGraphEdgeVec.push_back(e);

    return true;
}

void Source::partionNewNode(){
    int goodChioce = 0;
    int minLoadChioce = 0;
    int length = slaveNum;
    vector<pair<MID, double>> valueVec;
    for(auto it = newNode.begin(); it != newNode.end();it++) {
        valueVec.clear();
        VID nodeId = *it;
        double value = numeric_limits<double>::max();
        MID slaveId = nodeId % slaveNum;
        for(int i = 0;i < length;i++) {
            double tempValue = countSlaveValue(nodeId, (slaveId + i) % slaveNum);
            if(tempValue != 0.0) {
                valueVec.push_back(make_pair((slaveId + i) % slaveNum, tempValue));
            }
        }

        sort(valueVec.begin(), valueVec.end(), [](const pair<MID, double> &x, const pair<MID, double> &y) -> int {
            return x.second < y.second;
        });

        bool flag = false;
        for(int i = 0;i < valueVec.size();i++) {
            if(slaveToLoad[valueVec[i].first] + nodeNeighbourNum[nodeId] <= threshold) {
                nodeToSlave[nodeId] = valueVec[i].first;
                slaveToLoad[valueVec[i].first] += nodeNeighbourNum[nodeId];
                goodChioce++;
                flag = true;
                break;
            }
        }

        if(!flag) {
            nodeToSlave[nodeId] = minLoadSlave;
            slaveToLoad[minLoadSlave] += nodeNeighbourNum[nodeId];
            minLoad = slaveToLoad[minLoadSlave];
            minLoadChioce++;
            for(MID i = 0; i < slaveNum; i++) {
                if (slaveToLoad[i] < minLoad) {
                    minLoad = slaveToLoad[i];
                    minLoadSlave = i;
                }
            }
            threshold = minLoad * tolerancePlusOne;
        }
    }
}

double Source::countSlaveValue(VID nodeId, MID slaveId){
    double value = 0;
    double tempValue = 0;
    for(auto neighbourId : batchGraph[nodeId]) {
        if(nodeToSlave[neighbourId] == missingMId || nodeToSlave[neighbourId] == slaveId) {
            continue;
        }
        countEdgeValue(nodeId, slaveId, neighbourId, nodeToSlave[neighbourId], tempValue);
        value += tempValue;
    }
    return value;
}

bool Source::countEdgeValue(VID srcId, MID srcSlaveId, VID dstId, MID dstSlaveId, double& value){
    int srcNeighbourNum = nodeNeighbourNum[srcId];
    int dstNeighbourNum = nodeNeighbourNum[dstId];
    double srcToDstTransTime = max(srcNeighbourNum / bandWidth[srcSlaveId][0], srcNeighbourNum / bandWidth[dstSlaveId][1]);
    double dstToSrcTransTime = max(dstNeighbourNum / bandWidth[srcSlaveId][0], dstNeighbourNum / bandWidth[dstSlaveId][1]);
    double srcToDstTransCost = srcNeighbourNum * bandCost[srcSlaveId][0];
    double dstToSrcTransCost = dstNeighbourNum * bandCost[dstSlaveId][0];
    double tempValue = (srcToDstTransCost - dstToSrcTransCost) * (1 - theta) + (srcToDstTransTime - dstToSrcTransTime) * theta;
    if(tempValue <= 0) {
        value = srcToDstTransCost * (1 - theta) + srcToDstTransTime * theta;
        return true;
    } else {
        value = dstToSrcTransCost * (1 - theta) + dstToSrcTransTime * theta;
        return false;
    }
}

void Source::countEdgeNodeWeight(Edge& edge, double& srcWeight, double& dstWeight){
    int srcNeighbourNum = nodeNeighbourNum[edge.src];
    int dstNeighbourNum = nodeNeighbourNum[edge.dst];
    double srcToDstTransTime = max(srcNeighbourNum / bandWidth[edge.srcSlaveId][0], srcNeighbourNum / bandWidth[edge.dstSlaveId][1]);
    double dstToSrcTransTime = max(dstNeighbourNum / bandWidth[edge.srcSlaveId][0], dstNeighbourNum / bandWidth[edge.dstSlaveId][1]);
    double srcToDstTransCost = srcNeighbourNum * bandCost[edge.srcSlaveId][0];
    double dstToSrcTransCost = dstNeighbourNum * bandCost[edge.dstSlaveId][0];
    srcWeight = srcToDstTransCost * (1 - theta) + srcToDstTransTime * theta;
    dstWeight = dstToSrcTransCost * (1 - theta) + dstToSrcTransTime * theta;
    return ;
}

void Source::calculateEdgeSendDirection(std::vector<Edge>& resultVector){
    std::vector<std::vector<std::vector<Edge>>> edgeThreeDimesionVector(slaveNum, std::vector<std::vector<Edge>>(slaveNum, std::vector<Edge>()));
    std::vector<std::vector<std::map<VID, double>>> weightThreeDimesionMap(slaveNum, std::vector<std::map<VID, double>>(slaveNum, std::map<VID, double>()));

    double srcWeight = 0;
    double dstWeight = 0;
    for(Edge& edge : batchGraphEdgeVec) {
        if(edge.srcSlaveId == missingMId) {
            edge.srcSlaveId = nodeToSlave[edge.src];
        }
        if(edge.dstSlaveId == missingMId) {
            edge.dstSlaveId = nodeToSlave[edge.dst];
        }
        if(edge.srcSlaveId == edge.dstSlaveId) {
            edge.computeSlaveId = edge.srcSlaveId;
            sendEdgeVec.push_back(edge);
        } 
        else {
            edgeThreeDimesionVector[min(edge.srcSlaveId, edge.dstSlaveId)][max(edge.srcSlaveId, edge.dstSlaveId)].push_back(edge);

            countEdgeNodeWeight(edge, srcWeight, dstWeight);

            weightThreeDimesionMap[min(edge.srcSlaveId, edge.dstSlaveId)][max(edge.srcSlaveId, edge.dstSlaveId)][edge.src] = srcWeight;
            weightThreeDimesionMap[min(edge.srcSlaveId, edge.dstSlaveId)][max(edge.srcSlaveId, edge.dstSlaveId)][edge.src] = dstWeight;
        }
    }

    for(int i = 0;i < slaveNum;i++) {
        for(int j = i + 1; j < slaveNum;j++){
            bar(edgeThreeDimesionVector[i][j], weightThreeDimesionMap[i][j], resultVector);
        }
    }
}

void Source::bar(std::vector<Edge>& cutEdgeVector, std::map<VID, double>& edgeNodeWeight, std::vector<Edge>& resultVector){
    std::set<VID> sendNodeIdSet;
    
    for(Edge& edge : cutEdgeVector) {
        if(sendNodeIdSet.find(edge.src) == sendNodeIdSet.end() && sendNodeIdSet.find(edge.dst) == sendNodeIdSet.end() ) {
            if(edgeNodeWeight[edge.src] > edgeNodeWeight[edge.dst]) {
                edgeNodeWeight[edge.src] = edgeNodeWeight[edge.src] - edgeNodeWeight[edge.dst];
                edgeNodeWeight[edge.dst] = 0;
                sendNodeIdSet.insert(edge.dst);
                edge.computeSlaveId = edge.srcSlaveId;
            } else {
                edgeNodeWeight[edge.dst] = edgeNodeWeight[edge.dst] - edgeNodeWeight[edge.src];
                edgeNodeWeight[edge.src] = 0;
                sendNodeIdSet.insert(edge.src);
                edge.computeSlaveId = edge.dstSlaveId;
            }
            resultVector.push_back(edge);
        }
        if(sendNodeIdSet.find(edge.src) != sendNodeIdSet.end()) {
            edge.computeSlaveId = edge.dstSlaveId;
        }
        if(sendNodeIdSet.find(edge.dst) != sendNodeIdSet.end()) {
            edge.computeSlaveId = edge.srcSlaveId;
        }
        sendEdgeVec.push_back(edge);
    }
}

bool Source::isNewNode(VID &nodeId){
    return nodeToSlave[nodeId] == -1 ? true : false;
}

MID* Source::getNodeToSlave(){
    return nodeToSlave;
};

VID* Source::getNodeNeighbourNum(){
    return nodeNeighbourNum;
}

std::set<VID>*  Source::getNewNode(){
    return &newNode;
};

std::unordered_map<VID, std::set<VID>>* Source::getBatchGraph(){
    return &batchGraph;
};

std::vector<Edge>* Source::getBatchGraphEdgeVec(){
    return &batchGraphEdgeVec;
}

std::vector<Edge>* Source::getSendEdgeVec(){
    return &sendEdgeVec;
}


#ifdef _TEST_

#endif
