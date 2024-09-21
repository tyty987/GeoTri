#ifndef _STRUCT_HPP
#define _STRUCT_HPP

#include "macro.hpp"
#include "choice.hpp"

#include <fstream>
#include <cstring>
#include <limits>

#include <ctime>
#include <cstdlib>

using namespace std;

typedef unsigned int VID; 
typedef short MID;

const VID INVALID_VID(numeric_limits<VID>::max());
const MID INVALID_MID(numeric_limits<MID>::max());

struct Edge 
{
	static const unsigned short szAttr 	= 5;
	static const bool			TRUE	= true;
	static const bool 			FALSE	= false;

	VID	src;
	VID dst;
	MID srcSlaveId = -1;
	MID dstSlaveId = -1;
	MID computeSlaveId = -1;

	Edge(): src(INVALID_VID), dst(INVALID_VID), srcSlaveId(-1), dstSlaveId(-1), computeSlaveId(-1) {}
	Edge(VID iSrc, VID iDst): src(iSrc), dst(iDst) {}
	Edge(VID iSrc, VID iDst, MID iSrcSlaveId, MID iDstSlaveId, MID computeSlaveId): src(iSrc), dst(iDst), srcSlaveId(iSrcSlaveId), dstSlaveId(iDstSlaveId), computeSlaveId(computeSlaveId) {}

	inline bool operator==(const Edge& iEdge) const
	{
		return (src == iEdge.src) && (dst == iEdge.dst);
	}
	inline bool operator!=(const Edge& iEdge) const
	{	
		return !operator==(iEdge);
	}
	inline bool operator<(const Edge& iEdge) const
	{
		if(src < iEdge.src) {
			return true;
		} else if (src == iEdge.src) {
			return dst < iEdge.dst;
		}
	}

};

struct Node 
{
	static const unsigned short szAttrNode 	= 2;

	VID	index;
	MID slaveID;

	Node(): index(INVALID_VID), slaveID(-4) {}
	Node(VID iIndex, MID iSlaveId): index(iIndex), slaveID(iSlaveId) {}

	inline bool operator==(const Node& iNode) const 
	{
		return index == iNode.index;
	}
	inline bool operator!=(const Node& iNode) const 
	{	
		return !operator==(iNode);
	}
	inline bool operator<(const Node& iNode) const 
	{	
		return index < iNode.index;
	}
};

struct ElemCount
{
	static const unsigned short szAttr = 2;
	VID		vid;
	double	count;
	inline void setValue(VID iVid, double iCount)
	{
		vid = iVid;
		count = iCount;
	}
};

inline ostream& operator<<(ostream &os, const Edge &e){
	os << "(" << e.src << ", " << e.dst << ")" << "\t";
	return os;
}

class EdgeParser
{
_PRIVATE:
	fstream fp;

public:
	EdgeParser(){}
	EdgeParser(const char* filename): fp(filename, std::fstream::in) {}

	inline bool getEdge(Edge &oEdge)
	{
		if (fp.eof())
		{
			return false;
		}

		fp >> oEdge.src >> oEdge.dst;
		if (fp.eof())
		{
			return false;
		}

		return true;
	}

	inline void close() {
		fp.close();
	}

	inline void rewind()
	{
		fp.seekg(0, ios_base::beg);
	}
};



#endif
