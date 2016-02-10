#ifndef _ASSOCIATIONGRAPHOBSERVER_HPP_
#define _ASSOCIATIONGRAPHOBSERVER_HPP_

#include "GraphObserver.h"

#include "../Exceptions.h"
#include "../Clonable.h"

#include <vector>
#include <map>
#include <iostream>
#include <ostream>


namespace bpp
{
    
    /**
   * @brief Defines a Graph Associator. It is a template which follows
   * (subscribed to) a Graph.
   *
   * @author Thomas Bigot
   */
  
    // interface
    
    template <class N, class E>
    class AssociationGraphObserver:
    public virtual GraphObserver
    {
      
    private:
      typedef Graph::Node NodeGraphid;
      typedef Graph::Edge EdgeGraphid;
      
    public:    
      typedef unsigned int NodeIndex;
      typedef unsigned int EdgeIndex;
      
      
      /** @name Graph Relations Management
      *  Modificating the structure of the graph.
      */
      ///@{
      
      /**
      * Creates an orphaned node from a NodeClass object.
      * @param newNodeObject the N object associated to the node in the graph.
      * 
      */
      virtual void createNode(N* newNodeObject) = 0;
      
      
      /**
      * Creates an node linked to an existing node. Order of parameters match
      * the link method.
      * @param newNodeObject the N object associated to the node in the graph.
      * @param objectOriginNode existing node. In a directed graph: origin -> newNode.
      */
      virtual void createNode(N* objectOriginNode, N* newNodeObject) = 0;
      
      /**
      * Creates a link between two existing nodes.
      * If directed graph: nodeA -> nodeB.
      * @param nodeObjectA source node (or first node if undirected)
      * @param nodeObjectB target node (or second node if undirected)
      */
      virtual void link(N* nodeObjectA, N* nodeObjectB, E* edgeObject = 00) = 0;
      
      /**
      * Creates a link between two existing nodes.
      * If directed graph: nodeA -> nodeB.
      * @param nodeObjectA source node (or first node if undirected)
      * @param nodeObjectB target node (or second node if undirected)
      */
      virtual void unlink(N* nodeObjectA, N* nodeObjectB) = 0;
      
      /**
      * Deletes a node
      * @param nodeObject node to be deleted
      */
      virtual void deleteNode(N* nodeObject) = 0;
      
      ///@}
      
      /** @name Object Association
      *  Associate or dissociate N and E objects to pre-existing Graph Nodes and Graph Edges
      */
      ///@{
      
      /**
      * Associate a N or E object to a node or an edge in the graph.
      * @param nodeObject object to associate
      * @param node/edge existing node/edge to be associated
      */
      virtual void associateNode(N* nodeObject, NodeGraphid node) = 0;
      virtual void associateEdge(E* edgeObject, EdgeGraphid edge) = 0;
      
      /**
      * Dissociate a N or E object to a node or an edge in the graph.
      * @param nodeObject object to dissociate
      */
      virtual void forgetNode(N* nodeObject) = 0;
      virtual void forgetEdge(E* edgeObject) = 0;
      
      
      /**
      * Return the associated Node ID
      * @param nodeObject object which to return the node ID
      * @return a node ID
      */
      virtual NodeGraphid getNodeGraphid(const N* nodeObject) const = 0;
      
      /**
      * Return the associated Node ID
      * @param edgeObject object which to return the node ID
      * @return a node ID
      */
      virtual EdgeGraphid getEdgeGraphid(const E* edgeObject) const = 0;
      
      ///@}
      
      
       /** @name Object Indexation
      *  Get or set indexes to nodes and edges
      */
      ///@{
      
     /**
      * Return the associated Node index
      * @param nodeObject object which to return the node index
      * @return a node index
      */
      virtual NodeIndex getNodeIndex(const N* nodeObject) const = 0;
      virtual std::vector<NodeIndex> getNodeIndexes(std::vector<N*> nodeObjects) const = 0;
      
            
      /**
      * Return the associated Node index
      * @param edgeObject object which to return the node index
      * @return a node index
      */
      virtual EdgeIndex getEdgeIndex(const E* edgeObject) const = 0;
     
      /**
      * Set an index associated to a node
      * @param nodeObject object to which one want to set the index
      * @param index intex to be given, 0 to get the first free index
      * @return the given index
      */
      virtual NodeIndex setNodeIndex(const N* nodeObject, NodeIndex index = 0) = 0;
      
      /**
      * Set an index associated to an edge
      * @param edgeObject object to which one want to set the index
      * @param index intex to be given, 0 to get the first free index
      * @return the given index
      */
      virtual EdgeIndex setEdgeIndex(const E* edgeObject, EdgeIndex index = 0) = 0;
      
      /**
      * Return the associated Node, querying with an index
      * @param nodeIndex the index of the wanted node
      * @return N, a node object
      */
      virtual N* getNode(NodeIndex nodeIndex) const = 0;
      
            
      /**
      * Return the associated Node index
      * @param edgeIndex the index of the wanted edge
      * @return E, an edge object
      */
      virtual E* getEdge(EdgeIndex edgeIndex) const = 0;
      
      ///@}    
      
      /** @name Topology exploration
      *  These methodes of the graph concern the topology exploration.
      */
      ///@{
      
      /**
      * Get all the neighbors of a node in the graph.
      * @param node the node one wants to get its neighbors
      * @return a vector containing the neighbors
      */
      virtual std::vector<N*> getNeighbors(N* node) const = 0;
      virtual std::vector<NodeIndex> getNeighbors(NodeIndex node) const = 0;
      
      /**
      * In an directed graph, get all the neighbors which
      * are leaving a node in the graph.
      * @param node the node one wants to get its neighbors
      * @return a vector containing the outgoing neighbors
      */
      virtual std::vector<N*> getOutgoingNeighbors(N* node) const = 0;
      virtual std::vector<NodeIndex> getOutgoingNeighbors(NodeIndex node) const = 0;

      
      /**
      * In an directed graph, get all the neighbors which
      * are coming to a node in the graph.
      * @param node the node one wants to get its neighbors
      * @return a vector containing the incoming neighbors
      */
      virtual std::vector<N*> getIncomingNeighbors(N* node) const = 0;
      virtual std::vector<NodeIndex> getIncomingNeighbors(NodeIndex node) const = 0;

      
      /**
      * Get the leaves of a graph, ie, nodes with only one neighbor,
      * starting from a peculiar node.
      * @param node the starting node
      * @param maxDepth the maximum number of allowed depth, 0 means no max.
      * @return a vector containing the leaves
      */
      virtual std::vector<N*> getLeavesFromNode(N* node, unsigned int maxDepth) const = 0;
      
      /**
      * Get all the leaves objects of a graph, ie, nodes with only one neighbor,
      * @return a vector containing the leaves
      */
      virtual std::vector<N*> getAllLeaves() const = 0;
      
      /**
      * Get all the defined nodes of a graphO,
      * @return a vector containing the nodesObjects
      */
      virtual std::vector<N*> getAllNodes() const = 0;
      
     /**
      * Returns the Edge between two nodes nodeA -> nodeB
      * @param nodeA source node (if directed)
      * @param nodeB destination node (if directed)
      * @return the edge between these two nodes
      */
      virtual E* getEdgeLinking(N* nodeA, N* nodeB) const = 0;
      
      
      ///@}
      
    };
    
    
    template <class N, class E, class GraphImpl>
    class SimpleAssociationGraphObserver: public virtual AssociationGraphObserver<N,E>
    {

    public:
      typedef typename AssociationGraphObserver<N,E>::NodeIndex NodeIndex;
      typedef typename AssociationGraphObserver<N,E>::EdgeIndex EdgeIndex;
      
      typedef typename Graph::Node NodeGraphid;
      typedef typename Graph::Edge EdgeGraphid;
      
    private:
      
      //index management
      NodeIndex highestNodeIndex_;
      EdgeIndex highestEdgeIndex_;
      
      /**
      * List of nodes, stored at the same ID than the corresponding nodes
      * in the observed graph.
      */
      std::vector<N*> graphidToN_;
      
      /**
      * List of edges, stored at the same ID than the corresponding edges
      * in the bserved graph.
      */
      std::vector<E*> graphidToE_;
      
      /**
      * Can find a Node with the corresponding object.
      */
      std::map<N*,NodeGraphid> NToGraphid_;
      
      /**
      * Can find an Edge with the corresponding object.
      */
      std::map<E*,EdgeGraphid> EToGraphid_;
      
      
      /**
      * List of edges, stored at a given index.
      */
      std::vector<N*> indexToN_;
      
      /**
      * List of nodes, stored at a given index.
      */
      std::vector<E*> indexToE_;
      
      /**
      * Can find a Node index with the corresponding object.
      */
      std::map<N*,NodeIndex> NToIndex_;
      
      /**
      * Can find an Edge index with the corresponding object.
      */
      std::map<E*,EdgeIndex> EToIndex_;
      
      /**
       * defines a type of neighbors : incoming and/or outgoing
       */
      enum neighborType {INCOMING,OUTGOING,BOTH};
      
      /**
      * Get incoming / outgoing neighbors according to the enum type
      */
      std::vector<N*> getNeighbors_(N* nodeObject, neighborType type) const;
      
      
      /**
      * The observed Graph. Anytime this graph is observed,
      * the current object will be warned to take changes into account.
      */
      GraphImpl* subjectGraph_;
      
      /**
      * Set the observed Graph
      * @param subjectGraph the graph which is observed
      */
      void observe_(bpp::SimpleGraph subjectGraph);
      

      
      /**
      * Allows to create the proper observed subject class, thanks to the virtuality.
      * @param directed_p if the subject graph is directed
      */
      void allocateNewSubject(bool directed_p);
      
      
      
    public:
      
      /**
      * Constructor
      * @param directed is the graph directed
      */
      SimpleAssociationGraphObserver(bool directed = false);
      
      /**
      * Constructor
      * @param subjectGraph the graph which is observed
      */
      SimpleAssociationGraphObserver(bpp::SimpleGraph* subjectGraph = 00);
      
      /**
      * Copy Constructor
      * @param graphObserver the graphObserver to be copied
      */
      SimpleAssociationGraphObserver(bpp::SimpleAssociationGraphObserver<N,E,GraphImpl> const& graphObserver);
      
      /**
      * = Operator
      * @param graphObserver the graphObserver we want to copy the values
      */
      SimpleAssociationGraphObserver<N,E,GraphImpl>& operator=(bpp::SimpleAssociationGraphObserver<N,E,GraphImpl> const& graphObserver);
      
      /**
      * Destructor
      * @param graphObserver the graphObserver to be copied
      */
      ~SimpleAssociationGraphObserver();
      
      
      
      /**
      * clone function
      * @param graphObserver the graphObserver to be copied
      */
      #ifdef NO_VIRTUAL_COV
        Clonable*
      #else
        SimpleAssociationGraphObserver<N,E,GraphImpl>*
      #endif
      clone() const { return new SimpleAssociationGraphObserver<N,E,GraphImpl>(*this); };
      
      
      GraphImpl* getGraph() const;
      
      /**
      * This function is called to tell the observer that the subject
      * has changed and hence the observer has to take the changes
      * into account.
      */
      void update();
      
      
      
      /** @name Graph Relations Management
      *  Modificating the structure of the graph.
      */
      ///@{
      
      /**
      * Creates an orphaned node from a NodeClass object.
      * @param newNodeObject the N object associated to the node in the graph.
      * 
      */
      void createNode(N* newNodeObject);
      
      /**
      * Creates an node linked to an existing node. Order of parameters match
      * the link method.
      * @param newNodeObject the N object associated to the node in the graph.
      * @param objectOriginNode existing node. In a directed graph: origin -> newNode.
      */
      void createNode(N* objectOriginNode, N* newNodeObject);
      
      /**
      * Creates a link between two existing nodes.
      * If directed graph: nodeA -> nodeB.
      * @param nodeObjectA source node (or first node if undirected)
      * @param nodeObjectB target node (or second node if undirected)
      */
      void link(N* nodeObjectA, N* nodeObjectB, E* edgeObject = 00);
      
      /**
      * Creates a link between two existing nodes.
      * If directed graph: nodeA -> nodeB.
      * @param nodeObjectA source node (or first node if undirected)
      * @param nodeObjectB target node (or second node if undirected)
      */
      void unlink(N* nodeObjectA, N* nodeObjectB);
      
      /**
      * Deletes a node
      * @param nodeObject node to be deleted. The N node object given in argument is not deleted.
      */
      void deleteNode(N* nodeObject);
      
      
      ///@}
      
      /** @name Object Association
      *  Associate or dissociate N and E objects to pre-existing Graph Nodes and Graph Edges
      */
      ///@{
      
      /**
      * Associate a N or E object to a node or an edge in the graph.
      * @param nodeObject object to associate
      * @param node/edge existing node/edge to be associated
      */
      void associateNode(N* nodeObject, NodeGraphid graphNode);
      void associateEdge(E* edgeObject, EdgeGraphid graphEdge);
      
      /**
      * Dissociate a N or E object to a node or an edge in the graph.
      * @param nodeObject object to dissociate
      */
      void forgetNode(N* nodeObject);
      void forgetEdge(E* edgeObject);
      
      
      /**
      * Return the associated Node ID
      * @param nodeObject object which to return the node ID
      * @return a node ID
      */
      NodeGraphid getNodeGraphid(const N* nodeObject) const;

      /**
      * Return the associated Edge ID
      * @param edgeObject object which to return the node ID
      * @return a edge ID
      */
      EdgeGraphid getEdgeGraphid(const E* edgeObject) const;
      
      
      /**
       * Transforms an (a list of) id(s) into an (a list of) object(s)
       */
      N* getNodeFromGraphid(NodeGraphid) const;
      std::vector<N*> getNodesFromGraphid(std::vector<NodeGraphid>) const;
      E* getEdgeFromGraphid(EdgeGraphid) const;
      std::vector<E*> getEdgesFromGraphid(std::vector<EdgeGraphid>) const;
      
      
      ///@}   
      
      /** @name Object Indexation
      *  Get or set indexes to nodes and edges
      */
      ///@{
      
     /**
      * Return the associated Node index
      * @param nodeObject object which to return the node index
      * @return a node index
      */
      NodeIndex getNodeIndex(const N* nodeObject) const;
      std::vector<NodeIndex> getNodeIndexes(std::vector<N*> nodeObjects) const;
            
      /**
      * Return the associated Node index
      * @param edgeObject object which to return the node index
      * @return a node index
      */
      EdgeIndex getEdgeIndex(const E* edgeObject) const;
      
      /**
      * Set an index associated to a node
      * @param nodeObject object to which one want to set the index
      * @param index intex to be given, 0 to get the first free index
      * @return the given index
      */
      NodeIndex setNodeIndex(const N* nodeObject, NodeIndex index = 0);
      
      /**
      * Set an index associated to an edge
      * @param edgeObject object to which one want to set the index
      * @param index intex to be given, 0 to get the first free index
      * @return the given index
      */
      EdgeIndex setEdgeIndex(const E* edgeObject, EdgeIndex index = 0);
      
      
      /**
      * Return the associated Node index
      * @param nodeObject object which to return the node index
      * @return a node index
      */
      virtual N* getNode(NodeIndex nodeIndex) const;
      
            
      /**
      * Return the associated Node index
      * @param edgeObject object which to return the node index
      * @return a node index
      */
      virtual E* getEdge(EdgeIndex edgeIndex) const;
      ///@}    
      
      /** @name Topology exploration
      *  These methodes of the graph concern the topology exploration.
      */
      ///@{
      /**
      * Get all the neighbors of a node in the graph.
      * @param node the node one wants to get its neighbors
      * @return a vector containing the neighbors
      */
      std::vector<N*> getNeighbors(N* node) const;
      std::vector<NodeIndex> getNeighbors(NodeIndex node) const;

      /**
      * In an directed graph, get all the neighbors which
      * are leaving a node in the graph.
      * @param node the node one wants to get its neighbors
      * @return a vector containing the outgoing neighbors
      */
      std::vector<N*> getOutgoingNeighbors(N* node) const;
      std::vector<NodeIndex> getOutgoingNeighbors(NodeIndex node) const;

      /**
      * In an directed graph, get all the neighbors which
      * are coming to a node in the graph.
      * @param node the node one wants to get its neighbors
      * @return a vector containing the incoming neighbors
      */
      std::vector<N*> getIncomingNeighbors(N* node) const;
      std::vector<NodeIndex> getIncomingNeighbors(NodeIndex node) const;
      /**
      * Get the leaves of a graph, ie, nodes with only one neighbor,
      * starting from a peculiar node.
      * @param node the starting node
      * @param maxDepth the maximum number of allowed depth, 0 means no max.
      * @return a vector containing the leaves
      */
      std::vector<N*> getLeavesFromNode(N* node, unsigned int maxDepth) const;
      /**
      * Get all the leaves objects of a graph, ie, nodes with only one neighbor,
      * @return a vector containing the leaves
      */
      std::vector<N*> getAllLeaves() const;
      
      /**
      * Get all the defined nodes of a graph,
      * @return a vector containing the nodesObjects
      */
      std::vector<N*> getAllNodes() const;
      
      
      /**
    * Get nodes located at the extremities of an edge
    * @param edge an edge
    * @return a pair of the Nodes at each extremity of the edge
    *        example : N1--E1-->N2; getNodes(E1) will return (N1,N2);
    */
      std::pair<N*,N*> getNodes(E* edge) const;
      
            
      /**
      * Returns the Edge between two nodes nodeA -> nodeB
      * @param nodeA source node (if directed)
      * @param nodeB destination node (if directed)
      * @return the edge between these two nodes
      */
      E* getEdgeLinking(N* nodeA, N* nodeB) const;
      
      
      ///@}

      
      /** @name Function called by the subjectGraph
      *  These methodes are called by the subject graph to make this observer so fit the subject graph
      */
      ///@{
      
      /**
      * Delete unused object edges, since they have been deleted in the graph
      * @param edgesToDelete a vector of Edges to delete
      */
      void deletedEdgesUpdate(std::vector< EdgeGraphid >& edgesToDelete);
      
      /**
      * Delete unused object nodes, since they have been deleted in the graph
      * @param nodesToDelete a vector of N to delete
      */
      void deletedNodesUpdate(std::vector< NodeGraphid >& nodesToDelete);
      
      ///@}
      
     /** @name General Info
      *  General information about the graph
      */
      ///@{
      
      /**
      * Return the number of defined nodes, ie nodes that have a corresponding object
      * in this GraphObserver
      * @return the number of nodes
      */
      size_t getNumberOfNodes() const;
      
      /**
      * Return the number of defined leaves, ie leaves that have a corresponding object
      * in this GraphObserver
      * @return the number of leaves
      */
      size_t getNumberOfLeaves() const;
      
      ///@}
      
      
    };
    
//     class GraphObserverTools
//     {
//       static void outputToDot(SimpleAssociationGraphObserver<std::string,void>,std::ostream &out);
//     };
//     
    
template <class N, class E, class GraphImpl>
N* SimpleAssociationGraphObserver<N,E,GraphImpl>::getNodeFromGraphid(NodeGraphid node) const
{
  return graphidToN_.at(node); 
}

template <class N, class E, class GraphImpl>
E* SimpleAssociationGraphObserver<N,E,GraphImpl>::getEdgeFromGraphid(EdgeGraphid edge) const
{
  return graphidToE_.at(edge); 
}

template <class N, class E, class GraphImpl>
std::vector<N*> SimpleAssociationGraphObserver<N,E,GraphImpl>::getNodesFromGraphid(std::vector<NodeGraphid> nodes) const
{
  std::vector<N*> nodeObjects;
  for(typename std::vector<NodeGraphid>::iterator currNode = nodes.begin(); currNode != nodes.end(); currNode++)
  {
    if(*currNode > graphidToN_.size())
      continue;
    N* foundNodeObject = graphidToN_.at(*currNode);
    if(foundNodeObject == 00)
      continue;
    nodeObjects.push_back(foundNodeObject);
  }
  return nodeObjects;
}

template <class N, class E, class GraphImpl>
std::vector<typename AssociationGraphObserver<N,E>::NodeIndex> SimpleAssociationGraphObserver<N,E,GraphImpl>::getNodeIndexes(std::vector<N*> nodes) const
{
  std::vector<NodeIndex> nodeIndexes;
  for(typename std::vector<N*>::iterator currNode = nodes.begin(); currNode != nodes.end(); currNode++)
  {
    nodeIndexes.push_back(getNodeIndex(*currNode));
  }
  return nodeIndexes;
}

template <class N, class E, class GraphImpl>
std::vector<E*> SimpleAssociationGraphObserver<N,E,GraphImpl>::getEdgesFromGraphid(std::vector<EdgeGraphid> edges) const
{
  std::vector<N*> edgeObjects;
  for(typename std::vector<NodeGraphid>::iterator currEdge = edges.begin(); currEdge != edges.end(); currEdge++)
  {
    edgeObjects.push_back(graphidToE_.at(*currEdge));
  }
  return edgeObjects;
}

template <class N, class E, class GraphImpl>
void SimpleAssociationGraphObserver<N,E,GraphImpl>::createNode(N* nodeObject)
{
  NodeGraphid newGraphNode = subjectGraph_->createNode();
  associateNode(nodeObject, newGraphNode);
  
}

template <class N, class E, class GraphImpl>
SimpleAssociationGraphObserver<N,E,GraphImpl>::SimpleAssociationGraphObserver(bool directed_p):
  highestNodeIndex_(0),
  highestEdgeIndex_(0),
  graphidToN_(std::vector<N*>()),
  graphidToE_(std::vector<E*>()),
  NToGraphid_(std::map<N*,NodeGraphid>()),
  EToGraphid_(std::map<E*,NodeGraphid>()),
  indexToN_(std::vector<N*>()),
  indexToE_(std::vector<E*>()),
  NToIndex_(std::map<N*,typename AssociationGraphObserver<N,E>::NodeIndex>()),
  EToIndex_(std::map<E*,typename AssociationGraphObserver<N,E>::EdgeIndex>()),
  subjectGraph_(00)
{
  this->allocateNewSubject(directed_p);
}

template <class N, class E, class GraphImpl>
void SimpleAssociationGraphObserver<N,E,GraphImpl>::allocateNewSubject(bool directed_p)
{
    subjectGraph_ = new GraphImpl(directed_p);
    subjectGraph_->registerObserver(this);
}

template <class N, class E, class GraphImpl>
SimpleAssociationGraphObserver<N,E,GraphImpl>::~SimpleAssociationGraphObserver()
{
  this->subjectGraph_->unregisterObserver(this);
}

template <class N, class E, class GraphImpl>
SimpleAssociationGraphObserver<N,E,GraphImpl>::SimpleAssociationGraphObserver(SimpleAssociationGraphObserver<N,E,GraphImpl> const& graphObserver):
  highestNodeIndex_(0),
  highestEdgeIndex_(0),
  graphidToN_(graphObserver.graphidToN_),
  graphidToE_(graphObserver.graphidToE_),
  NToGraphid_(graphObserver.NToGraphid_),
  EToGraphid_(graphObserver.EToGraphid_),
  indexToN_(std::vector<N*>()),
  indexToE_(std::vector<E*>()),
  NToIndex_(std::map<N*,typename AssociationGraphObserver<N,E>::NodeIndex>()),
  EToIndex_(std::map<E*,typename AssociationGraphObserver<N,E>::EdgeIndex>()),
  subjectGraph_(graphObserver.subjectGraph_)
{
}


template <class N, class E, class GraphImpl>
SimpleAssociationGraphObserver<N,E,GraphImpl>& SimpleAssociationGraphObserver<N,E,GraphImpl>::operator=(SimpleAssociationGraphObserver<N,E,GraphImpl> const& graphObserver)
{
  this->highestEdgeIndex_ = graphObserver.highestEdgeIndex_;
  this->graphidToN_ = graphObserver.graphidToN_;
  this->graphidToE_ = graphObserver.graphidToE_;
  this->EToGraphid_ = graphObserver.EToGraphid_;
  this->NToGraphid_ = graphObserver.NToGraphid_;
  this->indexToN_ = graphObserver.indexToN_;
  this->indexToE_ = graphObserver.indexToE_;
  this->NToIndex_ = graphObserver.NToIndex_;
  this->EToIndex_ = graphObserver.EToIndex_;
  this->subjectGraph_ = graphObserver.subjectGraph_;
  return *this;
}

template <class N, class E, class GraphImpl>
void SimpleAssociationGraphObserver<N,E,GraphImpl>::createNode(N* objectOriginNode,N* newNodeObject)
{
  createNode(newNodeObject);
  link(objectOriginNode,newNodeObject);
}

template <class N, class E, class GraphImpl>
void SimpleAssociationGraphObserver<N,E,GraphImpl>::link(N* nodeObjectA, N* nodeObjectB, E* edgeObject)
{
  // checking the nodes
  typename std::map<N*,NodeGraphid>::iterator foundNodeA, foundNodeB;
  foundNodeA = NToGraphid_.find(nodeObjectA);
  foundNodeB = NToGraphid_.find(nodeObjectB);
  if(foundNodeA == NToGraphid_.end() || foundNodeB == NToGraphid_.end())
    throw Exception("One of the nodes is not in the graph observer.");
  
  //checking if the edge is not already in the GraphObserver
  if(edgeObject != 00 && EToGraphid_.find(edgeObject) != EToGraphid_.end())
    throw Exception("The given edge is already associated to a relation in the subjectGraph.");
  
  std::cout << "Trying to link node " << foundNodeA->second << " -> " << foundNodeB->second << std::endl;
  EdgeGraphid newGraphEdge = subjectGraph_->link(foundNodeA->second,foundNodeB->second);
  
  if(graphidToE_.size() < newGraphEdge+1)
    graphidToE_.resize(newGraphEdge+1);
  graphidToE_.at(newGraphEdge) = edgeObject;
  
  EToGraphid_[edgeObject] = newGraphEdge;
}

template <class N, class E, class GraphImpl>
void SimpleAssociationGraphObserver<N,E,GraphImpl>::unlink(N* nodeObjectA, N* nodeObjectB)
{
  //checking the nodes
  typename std::map<N*,NodeGraphid>::iterator foundNodeA, foundNodeB;
  foundNodeA = NToGraphid_.find(nodeObjectA);
  foundNodeB = NToGraphid_.find(nodeObjectB);
  if(foundNodeA == NToGraphid_.end() || foundNodeB == NToGraphid_.end())
    throw Exception("One of the nodes is not in the graph observer.");
  
  subjectGraph_->unlink(foundNodeA->second,foundNodeB->second);
}

template <class N, class E, class GraphImpl>
void SimpleAssociationGraphObserver<N,E,GraphImpl>::deletedEdgesUpdate(std::vector<EdgeGraphid>& edgesToDelete)
{
  for(typename std::vector<EdgeGraphid>::iterator currEdge = edgesToDelete.begin(); currEdge != edgesToDelete.end(); currEdge++){
    E* edgeObject = graphidToE_.at(*currEdge);
    graphidToE_.at(*currEdge) = 00;
    
    EToGraphid_.erase(edgeObject);
    
  }
}

template <class N, class E, class GraphImpl>
void SimpleAssociationGraphObserver<N,E,GraphImpl>::deletedNodesUpdate(std::vector<NodeGraphid>& nodesToDelete){
  for(typename std::vector<EdgeGraphid>::iterator currNode = nodesToDelete.begin(); currNode != nodesToDelete.end(); currNode++){
    N* nodeObject = graphidToN_.at(*currNode);
    graphidToN_.at(*currNode) = 00;
    
    NToGraphid_.erase(nodeObject);
    
  }
}

template <class N, class E, class GraphImpl>
void SimpleAssociationGraphObserver<N,E,GraphImpl>::associateNode(N* nodeObject, NodeGraphid graphNode)
{
  // nodes vector must be the right size. Eg: to store a node with
  // the ID 3, the vector must be of size 4: {0,1,2,3} (size = 4)
  graphidToN_.resize(subjectGraph_->getHighestNodeID()+1);

  // now storing the node
  graphidToN_.at(graphNode) = nodeObject;
  NToGraphid_[nodeObject] = graphNode;
}

template <class N, class E, class GraphImpl>
void SimpleAssociationGraphObserver<N,E,GraphImpl>::associateEdge(E* edgeObject, EdgeGraphid graphEdge)
{
  // nodes vector must be the right size. Eg: to store an edge with
  // the ID 3, the vector must be of size 4: {0,1,2,3} (size = 4)
  graphidToE_.resize(subjectGraph_->getHighestEdgeID()+1);
  
  // now storing the edge
  graphidToE_.at(graphEdge) = edgeObject;
  EToGraphid_[edgeObject] = graphEdge;
}

template <class N, class E, class GraphImpl>
typename AssociationGraphObserver<N,E>::EdgeIndex SimpleAssociationGraphObserver<N,E,GraphImpl>::setEdgeIndex(const E* edgeObject_p, typename AssociationGraphObserver<N,E>::EdgeIndex index)
{
  E* edgeObject = const_cast<E*>(edgeObject_p);
  //TODO: check if this object has already an index?
  if(index == 0)
    index = ++highestEdgeIndex_;
  // nodes vector must be the right size. Eg: to store an edge with
  // the index 3, the vector must be of size 4: {0,1,2,3} (size = 4)
  if(index > highestEdgeIndex_){
    highestEdgeIndex_ = index;
  }
  indexToE_.resize(highestEdgeIndex_+1);
  
  // now storing the edge
  indexToE_.at(index) = edgeObject;
  EToIndex_[edgeObject] = index;
  return index;
}

template <class N, class E, class GraphImpl>
void SimpleAssociationGraphObserver<N,E,GraphImpl>::forgetNode(N* nodeObject)
{
  typename std::map<N*,NodeGraphid>::iterator nodeToForget = NToGraphid_.find(nodeObject);
  graphidToN_.at(nodeToForget->second) = 00;
  NToGraphid_.erase(nodeToForget);
}

template <class N, class E, class GraphImpl>
typename AssociationGraphObserver<N,E>::NodeIndex SimpleAssociationGraphObserver<N,E,GraphImpl>::setNodeIndex(const N* nodeObject_p, typename AssociationGraphObserver<N,E>::NodeIndex index)
{
  N* nodeObject = const_cast<N*>(nodeObject_p);
  //TODO: check if this object has already an index?
  if(index == 0)
    index = ++highestNodeIndex_;
  // nodes vector must be the right size. Eg: to store a node with
  // the index 3, the vector must be of size 4: {0,1,2,3} (size = 4)
  if(index > highestNodeIndex_){
    highestNodeIndex_ = index;
  }
  indexToN_.resize(highestNodeIndex_+1);
  
  // now storing the node
  indexToN_.at(index) = nodeObject;
  NToIndex_[nodeObject] = index;
  return index;
}

template <class N, class E, class GraphImpl>
void SimpleAssociationGraphObserver<N,E,GraphImpl>::forgetEdge(E* edgeObject)
{
  typename std::map<E*,EdgeGraphid>::iterator edgeToForget = EToGraphid_.find(edgeObject);
  graphidToE_.at(edgeToForget->second) = 00;
  EToGraphid_.erase(edgeToForget);
}


template <class N, class E, class GraphImpl>
std::vector<N*> SimpleAssociationGraphObserver<N,E,GraphImpl>::getAllLeaves() const
{
  std::vector<N*> leavesToReturn;
  // fetching all the graph Leaves
  std::vector<NodeGraphid> graphLeaves = subjectGraph_->getAllLeaves();
  // testing if they are defined in this observer
  for(typename std::vector<NodeGraphid>::iterator currGraphLeave = graphLeaves.begin(); currGraphLeave != graphLeaves.end(); currGraphLeave++)
  {
    N* foundLeafObject = graphidToN_.at(*currGraphLeave);
    if(foundLeafObject != 00)
      leavesToReturn.push_back(foundLeafObject);
  }
  return leavesToReturn;
}

template <class N, class E, class GraphImpl>
std::vector<N*> SimpleAssociationGraphObserver<N,E,GraphImpl>::getAllNodes() const
{
  std::vector<N*> nodesToReturn;
  for(typename std::vector<N*>::const_iterator currNodeObject = graphidToN_.begin(); currNodeObject != graphidToN_.end(); currNodeObject++)
  {
    if(*currNodeObject != 00)
    {
      nodesToReturn.push_back(*currNodeObject);
    }
  }
  return nodesToReturn;
}

template <class N, class E, class GraphImpl>
size_t SimpleAssociationGraphObserver<N,E,GraphImpl>::getNumberOfNodes() const
{
  return NToGraphid_.size();
}
 
template <class N, class E, class GraphImpl>
size_t SimpleAssociationGraphObserver<N,E,GraphImpl>::getNumberOfLeaves() const
{
  return getAllLeaves().size();
}
 
template <class N, class E, class GraphImpl>
GraphImpl* SimpleAssociationGraphObserver<N,E,GraphImpl>::getGraph() const
{
  return subjectGraph_;
}

template <class N, class E, class GraphImpl>
void SimpleAssociationGraphObserver<N,E,GraphImpl>::deleteNode(N* nodeObject)
{
  // first deleting the node in the graph
  subjectGraph_->deleteNode(getNodeGraphid(nodeObject));
  // then forgetting
  forgetNode(nodeObject);
}

template <class N, class E, class GraphImpl>
typename SimpleAssociationGraphObserver<N,E,GraphImpl>::NodeGraphid SimpleAssociationGraphObserver<N,E,GraphImpl>::getNodeGraphid(const N* nodeObject) const
{
  typename std::map<N*,NodeGraphid>::const_iterator found = NToGraphid_.find(const_cast<N*>(nodeObject));
  if(found == NToGraphid_.end())
    throw Exception("Unexisting node object.");
  return found->second;
}

template <class N, class E, class GraphImpl>
typename AssociationGraphObserver<N,E>::NodeIndex SimpleAssociationGraphObserver<N,E,GraphImpl>::getNodeIndex(const N* nodeObject) const
{
  typename std::map<N*,typename AssociationGraphObserver<N,E>::NodeIndex>::const_iterator found = NToIndex_.find(const_cast<N*>(nodeObject));
  if(found == NToIndex_.end())
    throw Exception("Unexisting node object.");
  return found->second;
}

template <class N, class E, class GraphImpl>
typename SimpleAssociationGraphObserver<N,E,GraphImpl>::EdgeGraphid SimpleAssociationGraphObserver<N,E,GraphImpl>::getEdgeGraphid(const E* edgeObject) const
{
  typename std::map<E*,EdgeGraphid>::const_iterator found = EToGraphid_.find(const_cast<E*>(edgeObject));
  if(found == EToGraphid_.end())
    throw Exception("Unexisting edge object.");
  return found->second;
}

template <class N, class E, class GraphImpl>
typename AssociationGraphObserver<N,E>::EdgeIndex SimpleAssociationGraphObserver<N,E,GraphImpl>::getEdgeIndex(const E* edgeObject) const
{
  typename std::map<E*,typename AssociationGraphObserver<N,E>::EdgeIndex>::const_iterator found = EToIndex_.find(const_cast<E*>(edgeObject));
  if(found == EToIndex_.end())
    throw Exception("Unexisting edge object.");
  return found->second;
}

template <class N, class E, class GraphImpl>
std::vector< N* > SimpleAssociationGraphObserver<N,E,GraphImpl>::getNeighbors_(N* nodeObject, neighborType type) const
{
  NodeGraphid node = getNodeGraphid(nodeObject);
  
  // PHASE 1: getting the right neighbors
  std::vector<NodeGraphid> neighbors;
  switch(type){
    case OUTGOING:
      neighbors = subjectGraph_->getOutgoingNeighbors(node);
      break;
    case INCOMING:
      neighbors = subjectGraph_->getIncomingNeighbors(node);
      break;
    case BOTH:
      neighbors = subjectGraph_->getNeighbors(node);
  }
  return getNodesFromGraphid(neighbors);
}

template <class N, class E, class GraphImpl>
std::vector< N* > SimpleAssociationGraphObserver<N,E,GraphImpl>::getIncomingNeighbors(N* node) const
{
  return(getNeighbors_(node,INCOMING));
}

template <class N, class E, class GraphImpl>
std::vector<typename AssociationGraphObserver<N,E>::NodeIndex > SimpleAssociationGraphObserver<N,E,GraphImpl>::getIncomingNeighbors(NodeIndex node) const
{
  return getNodeIndexes(getIncomingNeighbors(getNode(node)));
}

template <class N, class E, class GraphImpl>
std::vector< N* > SimpleAssociationGraphObserver<N,E,GraphImpl>::getOutgoingNeighbors(N* node) const
{
  return getNeighbors_(node,OUTGOING);
}

template <class N, class E, class GraphImpl>
std::vector< typename AssociationGraphObserver<N,E>::NodeIndex > SimpleAssociationGraphObserver<N,E,GraphImpl>::getOutgoingNeighbors(NodeIndex node) const
{
  return getNodeIndexes(getOutgoingNeighbors(getNode(node)));
}

template <class N, class E, class GraphImpl>
std::vector< N* > SimpleAssociationGraphObserver<N,E,GraphImpl>::getNeighbors(N* node) const
{
  return getNeighbors_(node,BOTH);
}

template <class N, class E, class GraphImpl>
std::vector< typename AssociationGraphObserver<N,E>::NodeIndex > SimpleAssociationGraphObserver<N,E,GraphImpl>::getNeighbors(NodeIndex node) const
{
  return getNodeIndexes(getNeighbors(getNode(node)));
}


template <class N, class E, class GraphImpl>
std::vector< N* > SimpleAssociationGraphObserver<N,E,GraphImpl>::getLeavesFromNode(N* node, unsigned int maxDepth) const
{
  return getNodesFromGraphid(subjectGraph_->getLeavesFromNode(getNodeGraphid(node)));
}

template <class N, class E, class GraphImpl>
E* SimpleAssociationGraphObserver<N,E,GraphImpl>::getEdgeLinking(N* nodeA, N* nodeB) const
{
  return getEdge(subjectGraph_->getEdge(getNodeGraphid(nodeA),getNodeGraphid(nodeB)));
}

template <class N, class E, class GraphImpl>
N* SimpleAssociationGraphObserver<N,E,GraphImpl>::getNode(typename AssociationGraphObserver<N,E>::NodeIndex node) const
{
  return indexToN_.at(node); 
}

template <class N, class E, class GraphImpl>
E* SimpleAssociationGraphObserver<N,E,GraphImpl>::getEdge(typename AssociationGraphObserver<N,E>::EdgeIndex edge) const
{
  return indexToE_.at(edge); 
}

template <class N, class E, class GraphImpl>
std::pair<N*,N*> SimpleAssociationGraphObserver<N,E,GraphImpl>::getNodes(E* edge) const{
  std::pair<NodeGraphid,NodeGraphid> nodes = subjectGraph_->getNodes(getEdgeGraphid(edge));
  return std::pair<N*,N*>(getNodeFromGraphid(nodes.first),getNodeFromGraphid(nodes.second));
}
 
}

#else

namespace bpp {class AssociationGraphObserver;}

#endif