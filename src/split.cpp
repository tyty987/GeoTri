#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <set>
#include <vector>
#include <queue>
#include <iomanip>

using namespace std;
static int findDot(string filename, string symbol = ".") {
    int index = filename.find(symbol, 0);
    return index < filename.length() ? index : -1;
}

static string getFilename(string filename, string suffix) {
    string resultFilename(filename);
    resultFilename.insert(findDot(filename), suffix);
    return resultFilename;
}

static long getFileSize(string fileName) {
    ifstream file(fileName, ios::in);
    long size = 0;
    string str;
    while(getline(file, str)) {
        size++;
    }
    file.close();
    return size;
}

static void split_graph(string filename, double percent) {
    long edgeCount = getFileSize(filename);
    long dividePoint = edgeCount * percent;

    fstream file(filename, std::fstream::in);
    file.seekg(0, file.beg);
    
    string suffix = to_string((int)(10 * percent)) + to_string((int)(10 - 10 * percent));
    string tempFrontStr = "_true_front_" + suffix;
    string tempBackStr = "_true_back_" + suffix;

    string outFilenameFront = getFilename(filename, tempFrontStr);
    string outFilenameBack = getFilename(filename, tempBackStr);

    std::ofstream ofpFront(outFilenameFront, std::ios::out | std::ios::trunc);
    std::ofstream ofpBack(outFilenameBack, std::ios::out | std::ios::trunc);

    int src = 0, dst = 0;
    int count = 0;
    if (file.is_open()) {
        while (!file.eof()) {
            file >> src >> dst;
            count++;
            if(count <= dividePoint) {
                ofpFront << src << " " << dst << endl;
            } else {
                ofpBack << src << " " << dst << endl;
            }
        }
        file.close();
        ofpFront.close();
        ofpBack.close();
    }
    else {
        cout << "Unable to open file";
    }

    return;
}

static void split_graph_random(string filename, double percent) {
    long edgeCount = getFileSize(filename);
    long dividePoint = edgeCount * percent;

    fstream file(filename, std::fstream::in);
    file.seekg(0, file.beg);

    string suffix = to_string((int)(10 * percent)) + to_string((int)(10 - 10 * percent));
    string tempFrontStr = "_random_front_" + suffix;
    string tempBackStr = "_random_back_" + suffix;

    string outFilenameFront = getFilename(filename, tempFrontStr);
    string outFilenameBack = getFilename(filename, tempBackStr);


    std::ofstream ofpFront(outFilenameFront, std::ios::out | std::ios::trunc);
    std::ofstream ofpBack(outFilenameBack, std::ios::out | std::ios::trunc);

    int src = 0, dst = 0;
    int num = 0;
    std::map<int, std::vector<int>> graph;
    if (file.is_open()) {
        while (!file.eof()) {
            file >> src >> dst;
            graph[src].push_back(dst);
            graph[dst].push_back(src);
        }
        file.close();
    }

    srand((unsigned)time(NULL));
    int beginNode = 1;
    for(auto it = graph.begin();it != graph.end();it++) {
        if(rand() % 1000 < 10) {
            beginNode = it->first;
            break;
        }
    }

    std::set<int> nodeSet;
    nodeSet.insert(beginNode);
    std::set<int> usedNode;
    int count = 0;
    while(!nodeSet.empty()) {
        int nodeId = *(nodeSet.begin());
        nodeSet.erase(nodeSet.begin());
        usedNode.insert(nodeId);
        for(auto iter = graph[nodeId].begin();iter != graph[nodeId].end();iter++) {
            int src = nodeId;
            int dst = *iter;
            if(usedNode.find(dst) != usedNode.end()) {
                continue;
            }
            count++;
            if(count <= dividePoint) {
                ofpFront << src << " " << dst << endl;
            } else {
                ofpBack << src << " " << dst << endl;
            }
            if(usedNode.find(dst) == usedNode.end()) {
                nodeSet.insert(dst);
            }
        }
    }

    ofpFront.close();
    ofpBack.close();

    return;
}

static void count_new_old_node(string filename, double percent) {
    string suffix = to_string((int)(10 * percent)) + to_string((int)(10 - 10 * percent));
    string tempFrontStr = "_random_front_" + suffix;
    string tempBackStr = "_random_back_" + suffix;

    string outFilenameFront = getFilename(filename, tempFrontStr);
    string outFilenameBack = getFilename(filename, tempBackStr);

    fstream fileFront(outFilenameFront, std::fstream::in);
    fstream fileBack(outFilenameBack, std::fstream::in);
    
    std::set<int> oldNode;
    std::set<int> newNode;
    std::set<int> oldBackNode;
    map<int, set<int>> graphNodeNeightbougNum;

    int src = 0, dst = 0;
    if (fileFront.is_open()) {
        while (!fileFront.eof()) {
            if(graphNodeNeightbougNum.find(src) == graphNodeNeightbougNum.end()) {
                graphNodeNeightbougNum[src] = set<int>();
            }
             if(graphNodeNeightbougNum.find(dst) == graphNodeNeightbougNum.end()) {
                graphNodeNeightbougNum[dst] = set<int>();
            }
            fileFront >> src >> dst;
            oldNode.insert(src);
            oldNode.insert(dst);
            graphNodeNeightbougNum[src].insert(dst);
            graphNodeNeightbougNum[dst].insert(src);
        
        }
        fileFront.close();
    }
    else {
        cout << "Unable to open file";
    }

    int oldNertexNeighbourNum = 0;
    if (fileBack.is_open()) {
        while (!fileBack.eof()) {
            fileBack >> src >> dst;
            if(oldNode.find(src) == oldNode.end()) {
                newNode.insert(src);
            } else {
                oldBackNode.insert(src);
            }
            if(oldNode.find(dst) == oldNode.end()) {
                newNode.insert(dst);
            } else {
                oldBackNode.insert(dst);
            }
        }
        fileBack.close();
    }
    else {
        cout << "Unable to open file";
    }

    for(auto item : oldBackNode) {
        for(auto it : graphNodeNeightbougNum[item]) {
            if(oldBackNode.find(it) == oldBackNode.end()) {
                oldNertexNeighbourNum += 1;
            } else {
                oldNertexNeighbourNum += 0.5;
            }
        }
    }
    return;
}