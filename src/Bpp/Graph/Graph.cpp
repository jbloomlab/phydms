
#include <string>
#include <vector>
#include <sstream>
#include <algorithm>

#include "GraphObserver.h"
#include "Graph.h"
#include "../Exceptions.h"

using namespace bpp;
using namespace std;

SimpleGraph::SimpleGraph(bool directed_p):
  directed_(directed_p),
  observers_(set<GraphObserver*>()),
  numberOfNodes_(0),
  highestNodeID_(0),
  highestEdgeID_(0),
  nodeStructure_(nodeStructureType()),
  edgeStructure_(edgeStructureType()),
  
  root_(0)
{
 
}


void SimpleGraph::registerObserver(GraphObserver* observer)
{
  if(!observers_.insert(observer).second)
    throw(Exception("This GraphObserver was already an observer of this Graph"));;
}

void SimpleGraph::unregisterObserver(GraphObserver* observer)
{
  if(!observers_.erase(observer))
    throw(Exception("This GraphObserver was not an observer of this Graph"));
}

const SimpleGraph::Edge SimpleGraph::getEdge(SimpleGraph::Node nodeA, SimpleGraph::Node nodeB) const
{
   nodeStructureType::const_iterator firstNodeFound = nodeStructure_.find(nodeA);
   if(firstNodeFound == nodeStructure_.end())
     throw(Exception("The fist node was not the origin of an edge."));
   map<Node,Edge>::const_iterator secondNodeFound = firstNodeFound->second.first.find(nodeB);
   if(secondNodeFound == firstNodeFound->second.first.end())
     throw(Exception("The second node was not in a relation with the first one."));
   return(secondNodeFound->second);
}

const SimpleGraph::Edge SimpleGraph::getAnyEdge(SimpleGraph::Node nodeA, SimpleGraph::Node nodeB) const
{
   try
   {
     // trying in the given order A->B
     return getEdge(nodeA, nodeB);
   }
   catch(Exception e)
   {
     // didn’t work, hence trying in the opposite order B->A
     return getEdge(nodeB, nodeA);
   }
}

const SimpleGraph::Edge SimpleGraph::link(SimpleGraph::Node nodeA, SimpleGraph::Node nodeB)
{
  // the nodes must exist
  nodeMustExist_(nodeA, "first node");
  nodeMustExist_(nodeB, "second node");
  
  // which ID is available?
  unsigned int edgeID = ++highestEdgeID_;
  
  // writing the new relation to the structure
  linkInNodeStructure_(nodeA, nodeB, edgeID);
  if(!directed_){
    linkInNodeStructure_(nodeB, nodeA, edgeID);
  }
  linkInEdgeStructure_(nodeA, nodeB, edgeID);

  return edgeID;
}

void SimpleGraph::nodeMustExist_(Node node, string name) const
{
  if(nodeStructure_.find(node) == nodeStructure_.end())
  {
    ostringstream errMessage;
    errMessage << "This node must exist: " << node << " as " << name << ".";
    throw(Exception(errMessage.str()));
  }
}

void SimpleGraph::edgeMustExist_(SimpleGraph::Edge edge, string name) const
{
  if(edgeStructure_.find(edge) != edgeStructure_.end())
  {
    ostringstream errMessage;
    errMessage << "This edge must exist: " << edge << " as " << name << ".";
    throw(Exception(errMessage.str()));
  }
}


std::vector<SimpleGraph::Edge> SimpleGraph::unlink(const Node nodeA, const Node nodeB)
{
  // the nodes must exist
  nodeMustExist_(nodeA, "first node");
  nodeMustExist_(nodeB, "second node");
  
  // unlinking in the structure
  vector<Edge> deletedEdges; //what edges ID are affected by this unlinking
  deletedEdges.push_back(unlinkInNodeStructure_(nodeA, nodeB));
  for(vector<Edge>::iterator currEdgeToDelete = deletedEdges.begin(); currEdgeToDelete != deletedEdges.end();currEdgeToDelete++)
  {
    unlinkInEdgeStructure_(*currEdgeToDelete);
  }
  
  // telling the observers
  notifyDeletedEdges(deletedEdges);
  
  this->topologyHasChanged_();
  return deletedEdges;
}

void SimpleGraph::unlinkInEdgeStructure_(Edge edge)
{
  edgeStructureType::iterator foundEdge = edgeStructure_.find(edge);
  edgeStructure_.erase(foundEdge);
  this->topologyHasChanged_();
}

void SimpleGraph::linkInEdgeStructure_(SimpleGraph::Node nodeA, SimpleGraph::Node nodeB, SimpleGraph::Edge edge)
{
  edgeStructure_[edge] = pair<Node,Node>(nodeA,nodeB);
  this->topologyHasChanged_();
}



SimpleGraph::Edge SimpleGraph::unlinkInNodeStructure_(SimpleGraph::Node nodeA, SimpleGraph::Node nodeB)
{
  // Forward
  nodeStructureType::iterator nodeARow = nodeStructure_.find(nodeA);
  map<Node,Edge>::iterator foundForwardRelation = nodeARow->second.first.find(nodeB);
  Edge foundEdge = foundForwardRelation->second;
  nodeARow->second.first.erase(foundForwardRelation);
  // Backwards
  nodeStructureType::iterator nodeBRow = nodeStructure_.find(nodeB);
  map<Node,Edge>::iterator foundBackwardsRelation = nodeBRow->second.second.find(nodeA);
  nodeBRow->second.second.erase(foundBackwardsRelation);
  
  this->topologyHasChanged_();
  return foundEdge;
}

void SimpleGraph::linkInNodeStructure_(SimpleGraph::Node nodeA, SimpleGraph::Node nodeB, SimpleGraph::Edge edge)
{
  nodeStructure_.find(nodeA)->second.first.insert( pair<SimpleGraph::Node,SimpleGraph::Edge>(nodeB,edge));
  nodeStructure_.find(nodeB)->second.second.insert( pair<SimpleGraph::Node,SimpleGraph::Edge>(nodeA,edge));
  this->topologyHasChanged_();
}

const SimpleGraph::Node SimpleGraph::createNode()
{
  Node newNode = highestNodeID_++;
  nodeStructure_[newNode];
  numberOfNodes_++;
  this->topologyHasChanged_();
  return newNode;
}

const SimpleGraph::Node SimpleGraph::createNodeFromNode(SimpleGraph::Node origin)
{
  //origin must be an existing node
  nodeMustExist_(origin,"origin node");
  
  Node newNode = createNode();
  link(origin,newNode);
  this->topologyHasChanged_();
  return newNode;
}

const SimpleGraph::Node SimpleGraph::createNodeOnEdge(SimpleGraph::Edge edge)
{
  //origin must be an existing edge
  edgeMustExist_(edge,"");
  
  Node newNode = createNode();
  
  // determining the nodes on the border of the edge
  pair<Node,Node> nodes = edgeStructure_[edge];
  Node nodeA = nodes.first;
  Node nodeB = nodes.second;
  
  unlink(nodeA,nodeB);
  link(nodeA,newNode);
  link(newNode,nodeB);
  this->topologyHasChanged_();
  return newNode;
}



const SimpleGraph::Node SimpleGraph::createNodeFromEdge(SimpleGraph::Edge origin)
{
  //origin must be an existing edge
  edgeMustExist_(origin,"origin edge");
  
  // splitting the edge
  Node anchor = createNodeOnEdge(origin);
  
  Node newNode = createNodeFromNode(anchor);
  this->topologyHasChanged_();  
  return newNode;
}


void SimpleGraph::notifyDeletedEdges(vector< SimpleGraph::Edge > edgesToDelete) const
{
  for(set<GraphObserver*>::iterator currObserver = observers_.begin(); currObserver != observers_.end(); currObserver++)
  {
    (*currObserver)->deletedEdgesUpdate(edgesToDelete);
  }
}

void SimpleGraph::notifyDeletedNodes(vector< SimpleGraph::Node > nodesToDelete) const
{
  for(set<GraphObserver*>::iterator currObserver = observers_.begin(); currObserver != observers_.end(); currObserver++)
  {
    (*currObserver)->deletedNodesUpdate(nodesToDelete);
  }
}

std::vector< SimpleGraph::Node > SimpleGraph::getNeighbors_(SimpleGraph::Node node, bool outgoing) const
{
  nodeMustExist_(node,"");
  nodeStructureType::const_iterator foundNode = nodeStructure_.find(node);
  if(foundNode == nodeStructure_.end())
    throw(Exception("The requested node is not in the structure."));
  const std::map<Node,Edge> &forOrBack = (outgoing?foundNode->second.first:foundNode->second.second);
  vector<Node> result;
  for(map<Node,Edge>::const_iterator currNeighbor = forOrBack.begin(); currNeighbor!= forOrBack.end(); currNeighbor++)
  {
    result.push_back(currNeighbor->first);
  }
  
  return result;
}

vector< SimpleGraph::Node > SimpleGraph::getIncomingNeighbors(SimpleGraph::Node node) const
{
  return getNeighbors_(node,false);
}

vector< SimpleGraph::Node > SimpleGraph::getOutgoingNeighbors(SimpleGraph::Node node) const
{
  return getNeighbors_(node,true);
}

vector< SimpleGraph::Node > SimpleGraph::getNeighbors(SimpleGraph::Node node) const
{
  vector<SimpleGraph::Node> result;
  vector<SimpleGraph::Node> neighborsToInsert;
  neighborsToInsert = getNeighbors_(node,false);
  result.insert(result.end(),neighborsToInsert.begin(),neighborsToInsert.end());
  neighborsToInsert = getNeighbors_(node,true);
  result.insert(result.end(),neighborsToInsert.begin(),neighborsToInsert.end());
  return(result);
}

std::pair<SimpleGraph::Node, SimpleGraph::Node> SimpleGraph::getNodes(SimpleGraph::Edge edge) const
{
  edgeMustExist_(edge);
  edgeStructureType::const_iterator found = edgeStructure_.find(edge);
  //TODO Except if not found
  return found->second;
}


void SimpleGraph::deleteNode(SimpleGraph::Node node)
{
  //checking the node
  nodeMustExist_(node,"node to delete");
  isolate_(node);
  nodeStructureType::iterator found = nodeStructure_.find(node);
  nodeStructure_.erase(found);
  numberOfNodes_--;
  this->topologyHasChanged_();
}

void SimpleGraph::isolate_(SimpleGraph::Node node)
{
  vector<Node> oneighbors = getOutgoingNeighbors(node);
  for(vector<Node>::iterator currNeighbor = oneighbors.begin(); currNeighbor != oneighbors.end(); currNeighbor++){
    unlink(node,*currNeighbor);
  }
  vector<Node> ineighbors = getIncomingNeighbors(node);
  for(vector<Node>::iterator currNeighbor = ineighbors.begin(); currNeighbor != ineighbors.end(); currNeighbor++){
    unlink(*currNeighbor,node);
  }
}

unsigned int SimpleGraph::getHighestNodeID() const
{
  return highestNodeID_;
}


unsigned int SimpleGraph::getHighestEdgeID() const
{
  return highestEdgeID_;
}

vector<SimpleGraph::Node> SimpleGraph::getAllLeaves() const
{
  vector<Node> listOfLeaves;
  fillListOfLeaves_(root_,listOfLeaves,root_);
  return listOfLeaves;
}

void SimpleGraph::fillListOfLeaves_(Node startingNode, vector<Node>& foundLeaves, Node originNode, bool limitedRecursions, unsigned int maxRecursions) const
{
  const vector<Node> neighbors = getNeighbors(startingNode);
  if (neighbors.size() > 1)
  {
    if(!limitedRecursions || maxRecursions > 0)
      for(vector<Node>::const_iterator currNeighbor = neighbors.begin(); currNeighbor != neighbors.end(); currNeighbor++)
      {
        if (*currNeighbor != originNode)
          fillListOfLeaves_(*currNeighbor, foundLeaves, startingNode, limitedRecursions, maxRecursions-1);
      }
  }
  else
  {
    foundLeaves.push_back(startingNode);
  }
}


std::vector<SimpleGraph::Node> SimpleGraph::getLeavesFromNode(SimpleGraph::Node node,unsigned int maxDepth) const
{
  vector<Node> listOfLeaves;
  fillListOfLeaves_(node,listOfLeaves,node,(maxDepth!=0),maxDepth);
  return listOfLeaves;
}

void SimpleGraph::nodeToDot_(SimpleGraph::Node node, ostream& out,  std::set<std::pair<Node,Node> > &alreadyFigured) const
{
  bool theEnd = true;
  const std::map<Node,Edge> &children = nodeStructure_.at(node).first;
  for(map<Node,Edge>::const_iterator currChild = children.begin();currChild != children.end();currChild++)
  {
     if(alreadyFigured.find(pair<Node,Node>(node,currChild->first))!=alreadyFigured.end() || (!directed_ &&alreadyFigured.find(pair<Node,Node>(currChild->first,node))!=alreadyFigured.end()))
      continue;
    alreadyFigured.insert(pair<Node,Node>(node,currChild->first));
    theEnd = false;
    out << node << (directed_? " -> ":" -- ");
    nodeToDot_(currChild->first,out,alreadyFigured);
  }
  if(theEnd)
    out << node << ";\n    " ;
}

void SimpleGraph::outputToDot(ostream& out, std::string name) const
{
  out << (directed_?"digraph":"graph") << " "<< name <<" {\n    ";
  set<pair<Node,Node> > alreadyFigured;
  nodeToDot_(root_,out,alreadyFigured);
  out << "\r}" << endl;
}

bool SimpleGraph::isTree() const
{
  set<Graph::Node> metNodes;
  bool nodesAreMetOnlyOnce = nodesAreMetOnlyOnce_(root_,metNodes,root_);
  if(!nodesAreMetOnlyOnce)
    return false;
  // now they have only been met once, they have to be met at least once
  bool noNodeMissing = true;
  for(nodeStructureType::const_iterator currNode = nodeStructure_.begin(); noNodeMissing && currNode != nodeStructure_.end(); currNode++)
    noNodeMissing = (metNodes.find(currNode->first) != metNodes.end());
  return noNodeMissing;
}

bool SimpleGraph::nodesAreMetOnlyOnce_(Graph::Node node, set< Graph::Node >& metNodes, Graph::Node originNode) const
{
  //insert().second <=> not yet in the set
  bool neverMetANodeMoreThanOnce = metNodes.insert(node).second;
  vector<Graph::Node> neighbors = getOutgoingNeighbors(node);
  for(vector<Graph::Node>::iterator currNeighbor = neighbors.begin(); neverMetANodeMoreThanOnce && currNeighbor != neighbors.end(); currNeighbor++)
  {
    if(*currNeighbor==originNode)
      continue;
    neverMetANodeMoreThanOnce = nodesAreMetOnlyOnce_(*currNeighbor,metNodes,node);
  }
  return neverMetANodeMoreThanOnce;
}

void SimpleGraph::setRoot(Graph::Node newRoot)
{
  nodeMustExist_(newRoot,"new root");
  root_ = newRoot;
}

Graph::Node SimpleGraph::getRoot() const
{
  return(root_);
}


bool SimpleGraph::isDirected() const
{
  return(directed_);
}

void SimpleGraph::makeDirected()
{
  if(directed_)
    return;
  // save and clean the undirectedStructure
  nodeStructureType undirectedStructure = nodeStructure_;
  nodeStructure_.clear();
  // copy each relation once, without the reciprocal link
  // (first met, first kept)
  // eg: A - B in undirected is represented as A->B and B->A
  //     in directed, becomes A->B only
  std::set<pair<Node,Node> > alreadyConvertedRelations;
  for(nodeStructureType::iterator currNodeRow = undirectedStructure.begin(); currNodeRow != undirectedStructure.end(); currNodeRow++){
    Node nodeA = currNodeRow->first;
    for(map<Node,Edge>::iterator currRelation = currNodeRow->second.first.begin(); currRelation != currNodeRow->second.first.end(); currRelation++)
    {
      Node nodeB = currRelation->first;
      Edge edge = currRelation->second;
      if( !alreadyConvertedRelations.insert(pair<Node,Node>(min(nodeA,nodeB),max(nodeA,nodeB))).second)
        linkInNodeStructure_(nodeA,nodeB,edge);
    }
  }
  directed_=true;
  this->topologyHasChanged_();
}

void SimpleGraph::makeUndirected()
{
  if(!directed_)
    return;
  if(containsReciprocalRelations())
    throw Exception("Cannot make an undirected graph from a directed one containing reciprocal relations.");
  // save and clean the undirectedStructure
  nodeStructureType directedStructure = nodeStructure_;
  nodeStructure_.clear();
  // copy each relation twice, making the reciprocal link
  // eg: A - B in directed is represented as A->B
  //     in undirected, becomes A->B and B->A
  for(nodeStructureType::iterator currNodeRow = directedStructure.begin(); currNodeRow != directedStructure.end(); currNodeRow++){
    Node nodeA = currNodeRow->first;
    for(map<Node,Edge>::iterator currRelation = currNodeRow->second.first.begin(); currRelation != currNodeRow->second.first.end(); currRelation++)
    {
      Node nodeB = currRelation->first;
      Edge edge = currRelation->second;
      linkInNodeStructure_(nodeA,nodeB,edge);
      linkInNodeStructure_(nodeB,nodeA,edge);
    }
  }
  directed_=false;
  this->topologyHasChanged_();
}

bool SimpleGraph::containsReciprocalRelations() const
{
  if(!directed_)
    throw Exception("Cannot state reciprocal link in an undirected graph.");
  std::set<pair<Node,Node> > alreadyMetRelations;
  for(nodeStructureType::const_iterator currNodeRow = nodeStructure_.begin(); currNodeRow != nodeStructure_.end(); currNodeRow++){
    Node nodeA = currNodeRow->first;
    for(map<Node,Edge>::const_iterator currRelation = currNodeRow->second.first.begin(); currRelation != currNodeRow->second.first.end(); currRelation++)
    {
      Node nodeB = currRelation->first;
      if(!alreadyMetRelations.insert(pair<Node,Node>(min(nodeA,nodeB),max(nodeA,nodeB))).second)
        return true;
    }
  }
  return false;
}

void SimpleGraph::topologyHasChanged_() const
{
  // do nothing: a Graph does not care to be modified
}