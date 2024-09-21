#ifndef _SOURCE_HPP_
#define _SOURCE_HPP_

#include "struct.hpp"
#include "io.hpp"

#include <mpi.h>
#include <set>
#include <unordered_map>
#include <vector>
#include <iostream>
#include <algorithm>

class Source {
_PRIVATE:

    const Method_t method;

    const int slaveNum;

    double tolerancePlusOne;

    double tolerancePlusOnePointTwo;

    long* slaveToLoad;

    long minLoad;

    double threshold

    long minLoadSlave;

    MID* nodeToSlave;

    VID* nodeNeighbourNum;

    std::set<VID> newNode;

    std::unordered_map<VID, std::set<VID>> batchGraph;

    std::vector<Edge> batchGraphEdgeVec;

    std::vector<Edge> sendEdgeVec;

    VID maxVId;

    VID capacity;

    const MID missingMId;

    const double theta = 1.0;

    const int bandWidth[16][2] = {{8, 8}, {8, 8}, {8, 8}, {4, 4},
                                  {4, 4}, {1, 1}, {1, 1}, {1, 1},
                                  {4, 8}, {4, 4}, {4, 4}, {4, 4},
                                  {1, 1}, {1, 1}, {1, 1}, {1, 1}};
    const int bandCost[16][2] = {{8, 0}, {8, 0}, {8, 0}, {4, 0},
                                  {4, 0}, {4, 0}, {4, 0}, {8, 0},
                                  {4, 0}, {4, 0}, {4, 0}, {4, 0},
                                  {1, 0}, {1, 0}, {1, 0}, {1, 0}};

    Source(): method(Method_t::OPT), slaveNum(0), missingMId(-1) {;}

public:

    Source(int slaveNum, Method_t mode, double tolerance);

    ~Source();

    bool processEdge(Edge &iEdge, MID &oDstMID1, MID &oDstMID2);

    bool processBatchGraphEdge(Edge &iEdge);

    bool isNewNode(VID &nodeId);

    void partionNewNode();

    MID* getNodeToSlave();

    VID* getNodeNeighbourNum();

    std::set<VID>* getNewNode();

    std::unordered_map<VID, std::set<VID>>* getBatchGraph();

    std::vector<Edge>* getBatchGraphEdgeVec();

    VID getMaxVId();

    double countSlaveValue(VID nodeId, MID slaveId);

    bool countEdgeValue(VID srcId, MID srcSlaveId,  VID dstId, MID dstSlaveId, double& value);

    std::vector<Edge>* getSendEdgeVec();

    void calculateEdgeSendDirection(std::vector<Edge>& resultVector);
    void bar(std::vector<Edge>& cutEdgeVector, std::map<VID, double>& edgeNodeWeight, std::vector<Edge>& resultVector);
    void countEdgeNodeWeight(Edge& edge, double& srcWeight, double& dstWeight);
};

#endif
