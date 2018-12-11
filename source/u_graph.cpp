/***************************************************************************
                          Implementation of Graph class
                         -------------------------------
    last modified   : 21/6/2016
    author          : 2016 by Angus Kenny
                      (based on network.cpp by Dhananjay Thiruvady)
    libraries		    : .
    description		  : implementation of the graph data structure
 ***************************************************************************/


#include "u_graph.h"

/*
  Graph related function implementations
*/
// The graph destructor
U_Graph::~U_Graph(){
  for(int i=0;i<edges.size();i++)
    if(edges[i]!=NULL){
    delete edges[i];
    edges[i] = NULL;
  }
  for(int i=0;i<nodes.size();i++)
    if(nodes[i]!=NULL){
		delete nodes[i];
		nodes[i] = NULL;
  }
}
/*
const Node * Graph::getNode(int nodeID) const { return nodes[nodeID];}
const Edge * Graph::getEdge(long edgeID) const {return edges[edgeID];}
Node * Graph::getNode(int nodeID) { return nodes[nodeID];}
Edge * Graph::getEdge(long edgeID)   {return edges[edgeID];}
int Graph::getNumNodes() const { return nodes.size();}
int Graph::getNumEdges() const { return edges.size();}
int Graph::addEdge(int src,int tgt) {
	const int a=edges.size();
		addEdge(a,src,tgt);
	return a;
}
*/
// Method to add node to graph
void U_Graph::addNode(NodeID nodeID){
    Node *newNode = new Node(this,nodeID);
    nodes.push_back(newNode);
}

void U_Graph::insertNode(int pos){
  Node *newNode = new Node(this,pos);
  nodes.insert(nodes.begin()+pos, newNode);

  for (int i = pos; i < nodes.size(); i++){
    nodes[i]->setID(i);
  }
}

// Method to add edge to graph
void U_Graph::addEdge(EdgeID edgeID, NodeID n1, NodeID n2, EdgeID wt) {
  // order source and target
  NodeID srcID = std::min(n1,n2);
  NodeID tgtID = std::max(n1,n2);

  Edge * newEdge = new Edge(edgeID, nodes[srcID], nodes[tgtID], wt);
  edges.push_back(newEdge);
  nodes[srcID]->insertOutEdge(newEdge);
  nodes[tgtID]->insertInEdge(newEdge);
}

// Method to delete a batch of edges and reset edge IDs
// takes vector of type long with edge IDs to remove
void U_Graph::deleteEdges(std::vector<EdgeID> &edgeIDs){
  // sort in descending order so removing edge doesn't affect subsequent IDs
  std::sort(edgeIDs.rbegin(), edgeIDs.rend());

  // remove edges from edge vector
  for (int i = 0; i < edgeIDs.size(); i++){
    int src = edges[edgeIDs[i]]->getSrc()->getID();
    int tgt = edges[edgeIDs[i]]->getTgt()->getID();

    for (int e = 0; e < nodes[src]->getOutDegree(); ++e){
      if (nodes[src]->getOutEdge(e)->getID() == edgeIDs[i]){
        nodes[src]->deleteOutEdge(e);
        break;
      }
    }

    for (int e = 0; e < nodes[tgt]->getInDegree(); ++e){
      if (nodes[tgt]->getInEdge(e)->getID() == edgeIDs[i]){
        nodes[tgt]->deleteInEdge(e);
        break;
      }
    }
    edges[edgeIDs[i]] = edges.back();
    edges.pop_back();
    // edges.erase(edges.begin() + edgeIDs[i]);
  }

  // single pass to reset all of the edgeIDs
  for (int i = 0; i < edges.size(); i++){
    edges[i]->setID(i);
  }
}

// Method to delete a batch of nodes and reset node IDs
// takes {1,0} vector of int with length |N|, if 1 then remove node
void U_Graph::deleteNodes(std::vector<NodeID> &nodeIDs){

  // // remove nodes where vector equals 1 in reverse order
  // // so removing node doesnt affect subsequent IDs
  // for (int i = nodesToRemove.size()-1; i >= 0; i--){
  //   if (nodesToRemove[i] == 1)
  //     nodes.erase(nodes.begin() + i);
  // }

  // sort in descending order so removing node doesn't affect subsequent IDs
  std::sort(nodeIDs.rbegin(), nodeIDs.rend());

  for (int i = 0; i < nodeIDs.size(); ++i){
    nodes.erase(nodes.begin() + nodeIDs[i]);
    // nodes[nodeIDs[i]] = nodes.back();
    // nodes.pop_back();
  }

  // single pass to reset all of the nodeIDs
  for (int i = 0; i < nodes.size(); i++){
    nodes[i]->setID(i);
  }
}

// Method to get specific edge in graph (returns NULL pointer if no edge exists)
Edge *U_Graph::getEdgePair(int n1, int n2) {
  NodeID srcID = std::min(n1,n2);
  NodeID tgtID = std::max(n1,n2);
  PF1("-> getting edge pair.. ")
  for (int i=0; i < nodes[srcID]->getOutDegree(); ++i)
    if (nodes[srcID]->getOutEdge(i)->getTgt()->getID() == tgtID){
      PF1("found! ")
      return nodes[srcID]->getOutEdge(i);
    }
  PF1("edge (" << srcID << ", "<<tgtID <<") not found! ")
  return NULL;
}

void U_Graph::resizeNodes(int n)
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


const int Node::getConnected (std::vector<Node*> &connected){
  connected.clear();
  
    for (size_t prev_edge = 0; prev_edge < this->getInDegree(); ++prev_edge){
      int prev_idx = this->getInEdge(prev_edge)->getSrc()->getID();
      connected.push_back(g->getNode(prev_idx));
    }
    for (size_t next_edge = 0; next_edge < this->getOutDegree(); ++next_edge){
      int next_idx = this->getOutEdge(next_edge)->getTgt()->getID();
      connected.push_back(g->getNode(next_idx));
    }

  return connected.size();
}

const int Node::getEdges(std::vector<Edge*> &edges){

  edges.clear();

  for (int i = 0; i < getInDegree(); ++i){
    edges.push_back(getInEdge(i));
  }
  for (int i = 0; i < getOutDegree(); ++i){
    edges.push_back(getOutEdge(i));
  }
  return edges.size();
}

void Node::deleteInEdge(int edge_idx){
  inEdges[edge_idx] = inEdges.back();
  inEdges.pop_back();
}

void Node::deleteOutEdge(int edge_idx){
  outEdges[edge_idx] = outEdges.back();
  outEdges.pop_back();
}

/*
// Node class constructor
Node::Node(int id){
  this->id = id;
}
Node::~Node(){}
void Node::insertInEdge(Edge *in){ inEdges.push_back(in);}
void Node::insertOutEdge(Edge *out){ outEdges.push_back(out);}
int Node::getID() const {  return id;}
void Node::setID(int id){this->id = id;}
int Node::getInDegree() const {return inEdges.size();}
int Node::getOutDegree() const {return outEdges.size();}
int Node::getDegree() const {return this->getInDegree() + this->getOutDegree();}
const Edge * Node::getInEdge(int index) const { return inEdges[index];}
const Edge * Node::getOutEdge(int index) const { return outEdges[index];}
Edge * Node::getInEdge(int index)  { return inEdges[index];}
Edge * Node::getOutEdge(int index) { return outEdges[index];}
void Node::deleteEdges() { inEdges.clear(); outEdges.clear();}
bool Node::isChain() const{
  if (inEdges.size() == 1 && outEdges.size() == 1)
    return true;
  else
    return false;
}
*/
/*
  Edge related function implementations
*/
 /*
Edge::~Edge(){}
const Node *Edge::getSrc() const{ return src;}
const Node *Edge::getTgt() const { return tgt;}
Node *Edge::getSrc() { return src;}
Node *Edge::getTgt() { return tgt;}
int Edge::getSrcID() const { return src->getID();}
int Edge::getTgtID() const { return tgt->getID();}
int Edge::getID() const { return id;}
void Edge::setID(long id){this->id = id;}
*/
