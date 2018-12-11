/***************************************************************************
                          Implementation of Graph class
                         -------------------------------
    last modified   : 21/6/2016
    author          : 2016 by Angus Kenny
                      (based on network.cpp by Dhananjay Thiruvady)
    libraries		    : .
    description		  : implementation of the graph data structure
 ***************************************************************************/


#include "../include/graph.h"

/*
  Graph related function implementations
*/
// The graph destructor
Graph::~Graph(){
  for(int i=0;i<nodes.size();i++)
    if(nodes[i]!=NULL){
		delete nodes[i];
		nodes[i] = NULL;
  }
  for(int i=0;i<arcs.size();i++)
    if(arcs[i]!=NULL){
		delete arcs[i];
		arcs[i] = NULL;
  }
}
/*
const Node * Graph::getNode(int nodeID) const { return nodes[nodeID];}
const Arc * Graph::getArc(long arcID) const {return arcs[arcID];}
Node * Graph::getNode(int nodeID) { return nodes[nodeID];}
Arc * Graph::getArc(long arcID)   {return arcs[arcID];}
int Graph::getNumNodes() const { return nodes.size();}
int Graph::getNumArcs() const { return arcs.size();}
int Graph::addArc(int src,int tgt) {
	const int a=arcs.size();
		addArc(a,src,tgt);
	return a;
}
*/
// Method to add node to graph
void Graph::addNode(NodeID nodeID){
    Node *newNode = new Node(this,nodeID);
    nodes.push_back(newNode);
}


// Method to add arc to graph
void Graph::addArc(ArcID arcID, NodeID src, NodeID tgt) {
  Arc * newArc = new Arc(arcID, nodes[src], nodes[tgt]);
  arcs.push_back(newArc);
  nodes[src]->insertOutArc(newArc);
  nodes[tgt]->insertInArc(newArc);
}

// Method to delete a batch of arcs and reset arc IDs
// takes vector of type long with arc IDs to remove
void Graph::deleteArcs(std::vector<ArcID> &arcIDs){
  // sort in descending order so removing arc doesn't affect subsequent IDs
  std::sort(arcIDs.rbegin(), arcIDs.rend());

  // remove arcs from arc vector
  for (int i = 0; i < arcIDs.size(); i++){
    arcs.erase(arcs.begin() + arcIDs[i]);
  }

  // single pass to reset all of the arcIDs
  for (int i = 0; i < arcs.size(); i++){
    arcs[i]->setID(i);
  }
}

// Method to delete a batch of nodes and reset node IDs
// takes {1,0} vector of int with length |N|, if 1 then remove node
void Graph::deleteNodes(const std::vector<NodeID> &nodesToRemove){

  // remove nodes where vector equals 1 in reverse order
  // so removing node doesnt affect subsequent IDs
  for (int i = nodesToRemove.size()-1; i >= 0; i--){
    if (nodesToRemove[i] == 1)
      nodes.erase(nodes.begin() + i);
  }

  // single pass to reset all of the nodeIDs
  for (int i = 0; i < nodes.size(); i++){
    nodes[i]->setID(i);
  }
}

// Method to get specific arc in graph (returns NULL pointer if no arc exists)
const Arc *Graph::getArcPair(int srcID, int tgtID) const {
  for (int i=0; i < nodes[srcID]->getOutDegree(); i++)
    if (nodes[srcID]->getOutArc(i)->getTgtID() == tgtID)
      return nodes[srcID]->getOutArc(i);

  return NULL;
}

void Graph::resizeNodes(int n)
{
  if( n < (int)nodes.size()){
    for(int i=n;i<(int)nodes.size();++i) delete nodes[i];
    nodes.resize(n);
    return;
  }
  nodes.reserve(n);
  for(int i=(int) nodes.size();i < n;++i)
      nodes.push_back(new Node(this,i));
}


/*const std::vector<int> & Graph::getPreds(int vertex){
  return pred[vertex];
  }*/

/*
  Node related function implementations
*/


const int Node::getConnected (const int direction, std::vector<Node*> &connected){
  connected.clear();

  if (direction == BACKWARD){
    for (size_t prev_arc = 0; prev_arc < this->getInDegree(); ++prev_arc){
      int prev_idx = this->getInArc(prev_arc)->getSrcID();
      connected.push_back(g->getNode(prev_idx));
    }
  }
  else{
    for (size_t next_arc = 0; next_arc < this->getOutDegree(); ++next_arc){
      int next_idx = this->getOutArc(next_arc)->getTgtID();
      connected.push_back(g->getNode(next_idx));
    }
  }
  return connected.size();
}

const int Node::getArcs(const int direction, std::vector<Arc*> &arcs){

  arcs.clear();

  if (direction == BACKWARD){
    for (int i = 0; i < getInDegree(); ++i){
      arcs.push_back(getInArc(i));
    }
  }
  else{
    for (int i = 0; i < getOutDegree(); ++i){
      arcs.push_back(getOutArc(i));
    }
  }
  return arcs.size();
}

/*
// Node class constructor
Node::Node(int id){
  this->id = id;
}
Node::~Node(){}
void Node::insertInArc(Arc *in){ inArcs.push_back(in);}
void Node::insertOutArc(Arc *out){ outArcs.push_back(out);}
int Node::getID() const {  return id;}
void Node::setID(int id){this->id = id;}
int Node::getInDegree() const {return inArcs.size();}
int Node::getOutDegree() const {return outArcs.size();}
int Node::getDegree() const {return this->getInDegree() + this->getOutDegree();}
const Arc * Node::getInArc(int index) const { return inArcs[index];}
const Arc * Node::getOutArc(int index) const { return outArcs[index];}
Arc * Node::getInArc(int index)  { return inArcs[index];}
Arc * Node::getOutArc(int index) { return outArcs[index];}
void Node::deleteArcs() { inArcs.clear(); outArcs.clear();}
bool Node::isChain() const{
  if (inArcs.size() == 1 && outArcs.size() == 1)
    return true;
  else
    return false;
}
*/
/*
  Arc related function implementations
*/
// Arc class constructor
Arc::Arc(ArcID id_, Node *src_, Node *tgt_)
    : id(id_) {st[0]=src_->getID();st[1]=tgt_->getID();}//,src(src_->getID()), tgt(tgt_->getID()) {}
 /*
Arc::~Arc(){}
const Node *Arc::getSrc() const{ return src;}
const Node *Arc::getTgt() const { return tgt;}
Node *Arc::getSrc() { return src;}
Node *Arc::getTgt() { return tgt;}
int Arc::getSrcID() const { return src->getID();}
int Arc::getTgtID() const { return tgt->getID();}
int Arc::getID() const { return id;}
void Arc::setID(long id){this->id = id;}
*/
