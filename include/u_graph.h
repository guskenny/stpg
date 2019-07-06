/***************************************************************************
                            An Undirected Weighted Graph Class for the STPG
                         -------------------
    last modified   : 17/10/2018
    author       :  2018 by Angus Kenny
    libraries		  : .
    description		: contains a data structure for an undirected weighted graph
 ***************************************************************************/

#ifndef U_Graph_H
#define U_Graph_H

#include <iostream>
#include <new>
#include <vector>
#include <string>
#include <set>
#include <cstdio>
#include <cstring>
#include <algorithm>
#include <stdint.h>
#include <debug.h>

class Edge;
class Node;

typedef uint32_t EdgeID;	// 32 bit only (compared to 64 bit long)
typedef int32_t NodeID;		// 32 bit

class U_Graph{
public:
    // allow better control over memory usage by limiting ID type
    U_Graph() {}
    U_Graph(int nVertices) { resizeNodes(nVertices);} // ,const std::vector<std::vector<int> > &_pred);
    ~U_Graph();
    const Node *getNode(NodeID nodeID) const  { return nodes[nodeID];}
    const Edge *getEdge(EdgeID edgeID) const {return edges[edgeID];}
    Edge *getEdgePair(NodeID n1, NodeID n2);
    Node *getNode(NodeID nodeID) { return nodes[nodeID];}
    Edge *getEdge(EdgeID edgeID) {return edges[edgeID];}
    NodeID getNumNodes() const { return nodes.size();}
    EdgeID getNumEdges() const  { return edges.size();}
    void resizeNodes(NodeID n);
    void addNode(NodeID nodeID);
    void insertNode(int pos);
    int addEdge(NodeID n1,NodeID n2,EdgeID wt=1) {
	addEdge(edges.size(),n1,n2,wt); return edges.size()-1;
    }
    void addEdge(EdgeID edgeID, NodeID n1, NodeID n2, EdgeID wt);
    void deleteEdges(std::vector<EdgeID> &edgeIDs);
    void deleteNodes(std::vector<NodeID> &nodeIDs);
    //AE: This looks like a really bad idea - why store predecessors
    //    here when we already have them stored as  inEdges for each node?
    //    Seems to increase memory use for little/no benefit
    //OK: Uncommenting until replaced
    const std::vector<int> & getPreds(int vertex)
     {   return pred[vertex]; }

private:
    std::vector<Node*> nodes;
    std::vector<Edge*> edges;
    //OK: restored
     std::vector<std::vector<int> > pred;
};



class Node{
  public:
    Node(U_Graph *g_,NodeID id_) : g(g_),id(id_),orig_id(id_),terminal(false) {}
    ~Node() {}


    void insertInEdge(Edge *in){ inEdges.push_back(in);}
    void insertOutEdge(Edge *out){ outEdges.push_back(out);}
    NodeID getID() const {  return id;}
    NodeID getOrigID() const {  return orig_id;}
    void setID(NodeID id){this->id = id;}
    int getInDegree() const {return inEdges.size();}
    int getOutDegree() const {return outEdges.size();}
    int getDegree() const {return this->getInDegree() + this->getOutDegree();}
    const Edge * getInEdge(int index) const { return inEdges[index];}
    const Edge * getOutEdge(int index) const { return outEdges[index];}
    Edge * getInEdge(int index)  { return inEdges[index];}
    Edge * getOutEdge(int index) { return outEdges[index];}
    const int getEdges (std::vector<Edge*> &edges);
    const int getConnected (std::vector<Node*> &connected);
    const int getConnectedIDs (std::vector<int> &connected);
    void deleteInEdge(int edge_idx);
    void deleteOutEdge(int edge_idx);
    void deleteEdges() { inEdges.clear(); outEdges.clear();}
    bool isChain() const{
    return ((inEdges.size() ==1) && (outEdges.size() == 1));
    }
    void setTerm(bool val){terminal = val;}
    void setKey(bool val){key = val;}
    bool isTerm()const{return terminal;}
    bool isKey()const{return key;}

  private:
    bool terminal;
    bool key;
    U_Graph *g;
    NodeID id;
    NodeID orig_id;
    std::vector<Edge*> inEdges;
    std::vector<Edge*> outEdges;

};

// undirected edges are arcs where the "source" node is the lowest ID end point,
// this way it is easier to retrieve a specific edge knowing only the two end points
class Edge{
  public:
    Edge(EdgeID id_, Node *src_, Node *tgt_, EdgeID wt_): id(id_), wt(wt_),src(src_),tgt(tgt_){} //: id(id_),src(src_->getID()), tgt(tgt_->getID()) {}
    // Edge(EdgeID id_, NodeID src_, NodeID tgt_, EdgeID wt_) : id(id_),wt(wt_),src(src_),tgt(tgt_){}
    ~Edge() {}
    //const Node *getSrc() const {return src;}
    //const Node *getTgt() const {return tgt;}
    //Node *getSrc() {return src;}
    //Node *getTgt() {return tgt;}
    int getSrcID()  { return src->getID();}
    int getTgtID()  { return tgt->getID();}
    Node* getSrc() const { return src;}
    Node* getTgt() const { return tgt;}
    int getWt() const {return wt;}
    EdgeID getID() const { return id;}
    void setID(EdgeID id_) {id=id_;}
    void setWt(EdgeID new_wt) {wt = new_wt;}
    void setSrc(Node* new_src);
    void setTgt(Node* new_tgt);

  protected:
    Node *src;
    Node *tgt;
    EdgeID id;
    EdgeID wt;          // edge weight
    //NodeID src,tgt;
    //Node *src,*tgt;
};

#endif
