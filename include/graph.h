/***************************************************************************
                            A Graph Class
                         -------------------
    last modified   : 21/6/2016
    author       :  2016 by Angus Kenny
                      (based on network.h by Dhananjay Thiruvady)
    libraries		  : .
    description		: contains a data structure for a graph
 ***************************************************************************/

#ifndef Graph_H
#define Graph_H

#include <iostream>
#include <new>
#include <vector>
#include <string>
#include <set>
#include <cstdio>
#include <cstring>
#include <algorithm>
#include <stdint.h>

class Arc;
class Node;

typedef uint32_t ArcID;	// 32 bit only (compared to 64 bit long)
typedef int32_t NodeID;		// 32 bit
const int BACKWARD = 0;
const int FORWARD = 1;
const int NUM_DIRS = 2;

class Graph{
public:
    // allow better control over memory usage by limiting ID type
    Graph() {}
    Graph(int nVertices) { resizeNodes(nVertices);} // ,const std::vector<std::vector<int> > &_pred);
    ~Graph();
    const Node *getNode(NodeID nodeID) const  { return nodes[nodeID];}
    const Arc *getArc(ArcID arcID) const {return arcs[arcID];}
    const Arc *getArcPair(NodeID srcID, NodeID tgtID) const;
    Node *getNode(NodeID nodeID) { return nodes[nodeID];}
    Arc *getArc(ArcID arcID) {return arcs[arcID];}
    NodeID getNumNodes() const { return nodes.size();}
    ArcID getNumArcs() const  { return arcs.size();}
    void resizeNodes(NodeID n);
    void addNode(NodeID nodeID);

    int addArc(NodeID src,NodeID tgt) {
	addArc(arcs.size(),src,tgt); return arcs.size()-1;
    }
    void addArc(ArcID arcID, NodeID src, NodeID tgt);
    void deleteArcs(std::vector<ArcID> &arcIDs);
    void deleteNodes(const std::vector<NodeID> &nodesToRemove);
    //AE: This looks like a really bad idea - why store predecessors
    //    here when we already have them stored as  inArcs for each node?
    //    Seems to increase memory use for little/no benefit
    //OK: Uncommenting until replaced
    const std::vector<int> & getPreds(int vertex)
     {   return pred[vertex]; }

private:
    std::vector<Node*> nodes;
    std::vector<Arc*> arcs;
    //OK: restored
     std::vector<std::vector<int> > pred;
};

class Arc{
  public:
    Arc(ArcID id_, Node *src_, Node *tgt_); //: id(id_),src(src_->getID()), tgt(tgt_->getID()) {}
    Arc(ArcID id_, NodeID src_, NodeID tgt_) : id(id_) {st[0]=src_;st[1]=tgt_;}
    ~Arc() {}
    //const Node *getSrc() const {return src;}
    //const Node *getTgt() const {return tgt;}
    //Node *getSrc() {return src;}
    //Node *getTgt() {return tgt;}
    int getSrcID() const { return st[0];}
    int getTgtID() const { return st[1];}
    ArcID getID() const { return id;}
    void setID(ArcID id_) {id=id_;}

  protected:
    NodeID st[2];		// force placing 2 x 32 bit into 1 64 bit word
    ArcID id;
    //NodeID src,tgt;
    //Node *src,*tgt;
};


class Node{
  public:
    Node(Graph *g_,NodeID id_) : g(g_),id(id_) {}
    ~Node() {}
    void insertInArc(Arc *in){ inArcs.push_back(in->getID());}
    void insertOutArc(Arc *out){ outArcs.push_back(out->getID());}
    NodeID getID() const {  return id;}
    void setID(NodeID id){this->id = id;}
    int getInDegree() const {return inArcs.size();}
    int getOutDegree() const {return outArcs.size();}
    int getDegree() const {return this->getInDegree() + this->getOutDegree();}
    const Arc * getInArc(int index) const { return g->getArc(inArcs[index]);}
    const Arc * getOutArc(int index) const { return g->getArc(outArcs[index]);}
    Arc * getInArc(int index)  { return g->getArc(inArcs[index]);}
    Arc * getOutArc(int index) { return g->getArc(outArcs[index]);}
    const int getArcs (const int direction, std::vector<Arc*> &arcs);
    const int getConnected (const int direction, std::vector<Node*> &connected);
    void deleteArcs() { inArcs.clear(); outArcs.clear();}
    bool isChain() const{
	return (inArcs.size() == 1 && outArcs.size() == 1);
    }

  private:
    Graph *g;
    NodeID id;
    std::vector<ArcID> inArcs;
    std::vector<ArcID> outArcs;

};

#endif
