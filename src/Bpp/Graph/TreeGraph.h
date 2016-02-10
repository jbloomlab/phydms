#ifndef _TREEGRAPH_H_
#define _TREEGRAPH_H_

#include <vector>

#include "Graph.h"

#include "../Exceptions.h"

namespace bpp
{
  class TreeGraph:
  public virtual Graph
  {
  public:
    
    //TODO Incomplete specification. Waiting the basic implementation is functional
    // to specify the interface.
    
    /**
     * Is the graph a tree? A tree must be acyclic and with no isolated node.
     * @return true if valid tree
     */
    virtual bool isValid() const = 0;
  };
  
  template <class GraphImpl>
  class SimpleTreeGraph:
  public virtual TreeGraph,
  public virtual GraphImpl
  {
  private:
    /**
     * Is the graph a tree? Set to false when structure is modified, true after validation.
     */
    mutable bool isValid_;
    
    // unvalidate the tree
    virtual void topologyHasChanged_() const;
    
    // will throw an exception if the tree is not valid
    void mustBeValid_() const;
    
    // will throw an exception if the tree is not rooted
    void mustBeRooted_() const;
    
    // test the validity of the tree
    bool validate_() const;
    
    /**
     * Reorient all the edges starting from a node:
     * the father node becomes a son, and so on.
     */
    void propagateDirection_(Node node);
    
    
    
    // recursive function for getSubtreeNodes
    void fillSubtreeMetNodes_(std::vector<Graph::Node>& metNodes, Graph::Node localRoot) const;
    
    // recursive function for getSubtreeEdges
    void fillSubtreeMetEdges_(std::vector<Graph::Edge>& metEdges, Graph::Node localRoot) const;
    
    
  public:
      
    SimpleTreeGraph();
    SimpleTreeGraph(bool rooted=true);
    
    /**
     * Is the graph a tree? A tree must be acyclic and with no isolated node.
     * @return true if valid tree
     */
    bool isValid() const;
    
    /**
     * Is the tree rooted?
     * @return true if rooted
     */
    bool isRooted() const;
    
    /**
     * Get the father node of a node in a rooted tree
     * @return the father node
     */
    Graph::Node getFather(Graph::Node node) const;
    
    /**
     * Get the branch leading to the father in a rooted tree
     * @return the branch between a node and its father
     */
    Graph::Edge getBranchToFather(Graph::Node node) const;
    
    /**
     * Get the father node of a node in a rooted tree
     */
    bool hasFather(Graph::Node node) const;
    
    /**
     * Get the father node of a node in a rooted tree
     */
    std::vector<Graph::Node> getSons(Graph::Node node) const;
    
    /**
     * Get the father node of a node in a rooted tree
     */
    void setFather(Graph::Node node, Graph::Node fatherNode);
    
    /**
     * Get the father node of a node in a rooted tree
     */
    void addSon(Graph::Node node, Graph::Node sonNode);
    
    /**
     * Remove all the sons
     */
    std::vector<Graph::Node> removeSons(Graph::Node node);
    
    /**
     * Remove one son
     */
    void removeSon(Graph::Node node, Graph::Node son);
    
    /**
     * Re-root the tree with the new root
     */
    void rootAt(Graph::Node newRoot);
    
    /**
     * Set the tree to its flat unrooted version.
     * As an algorithmical convenience, a root node is kept, but it has
     * no logical significance.
     */
    void unRoot(bool joinRootSons);
    
    /**
     * Set a node as a new outgroup in a rooted tree, will make a root between
     * the given node and its father.
     */
    void setOutGroup(Graph::Node newOutGroup);
    
    /**
     * Get all the nodes of a subtree
     */
    std::vector<Graph::Node> getSubtreeNodes(Graph::Node localRoot) const;
    
    /**
     * Get all the branches of a subtree
     */
    std::vector<Graph::Edge> getSubtreeEdges(Graph::Node localRoot) const;
    
    
    
    /////FROM TREETOOLS & TREETOOLS COMPAT
    
    
    std::vector<Graph::Node> getNodePathBetweenTwoNodes(Graph::Node nodeA, Graph::Node nodeB, bool includeAncestor = true) const;
    std::vector<Graph::Edge> getEdgePathBetweenTwoNodes(Graph::Node nodeA, Graph::Node nodeB) const;
    
  };
  
  template <class GraphImpl>
  SimpleTreeGraph<GraphImpl>::SimpleTreeGraph(bool rooted):
  GraphImpl(rooted),
  isValid_(false)
  {
  }
  
  
  template <class GraphImpl>
  bool SimpleTreeGraph<GraphImpl>::isValid() const
  {
    return (isValid_ || validate_());
  }
  
  template <class GraphImpl>
  Graph::Node SimpleTreeGraph<GraphImpl>::getFather(Graph::Node node) const
  {
    mustBeValid_();
    mustBeRooted_();
    std::vector<Graph::Node> incomers = GraphImpl::getIncomingNeighbors(node);
    if(incomers.size() != 1)
      throw Exception("SimpleTreeGraph<GraphImpl>::getFather: more than one father. Should never happen since validity has been controled. Please report this bug.");
    if(incomers.size() == 0)
      throw Exception("SimpleTreeGraph<GraphImpl>::getFather: node has no father.");
    return *incomers.begin();
  }
  
  template <class GraphImpl>
  Graph::Edge SimpleTreeGraph<GraphImpl>::getBranchToFather(Graph::Node node) const
  {
    Node father = getFather(node);
    return GraphImpl::getBranch(father,node);
  }
  
  template <class GraphImpl>
  bool SimpleTreeGraph<GraphImpl>::hasFather(Graph::Node node) const
  {
    mustBeValid_();
    mustBeRooted_();
    std::vector<Graph::Node> incomers = SimpleGraph::getIncomingNeighbors(node);
    return incomers.size() == 1;
  }
  
  template <class GraphImpl>
  void SimpleTreeGraph<GraphImpl>::mustBeRooted_() const
  {
    if(!isRooted())
      throw Exception("SimpleTreeGraph: The tree must be rooted.");
  }
  
  template <class GraphImpl>
  void SimpleTreeGraph<GraphImpl>::mustBeValid_() const
  {
    if(!isValid())
      throw Exception("SimpleTreeGraph: The tree is not valid.");
  }
  
  template <class GraphImpl>
  bool SimpleTreeGraph<GraphImpl>::isRooted() const
  {
    return(SimpleGraph::isDirected());
  }
  
  template <class GraphImpl>
  bool SimpleTreeGraph<GraphImpl>::validate_() const
  {
    isValid_ = SimpleGraph::isTree();
    return(isValid_);
  }
  
  template <class GraphImpl>
  void SimpleTreeGraph<GraphImpl>::topologyHasChanged_() const
  {
    isValid_ = false;
  }
  
  template <class GraphImpl>
  void SimpleTreeGraph<GraphImpl>::rootAt(Graph::Node newRoot)
  {
    GraphImpl::makeDirected();
    // set the new root on the Graph
    GraphImpl::setRoot(newRoot);
    // change edge direction between the new node and the former one
    propagateDirection_(newRoot);
  }
  
  template <class GraphImpl>
  void SimpleTreeGraph<GraphImpl>::propagateDirection_(Graph::Node node)
  {
    if(hasFather(node)){
      Node father = getFather(node);
      GraphImpl::unlink(father,node);
      GraphImpl::link(node,father);
      propagateDirection_(father);
    }
    
  }
  
  template <class GraphImpl>
  void SimpleTreeGraph<GraphImpl>::setFather(Graph::Node node, Graph::Node fatherNode)
  {
    if(hasFather(node))
      GraphImpl::unlink(getFather(node),node);
    GraphImpl::link(fatherNode,node);
  }
  
  template <class GraphImpl>
  void SimpleTreeGraph<GraphImpl>::addSon(Graph::Node node, Graph::Node sonNode)
  {
    GraphImpl::link(node,sonNode);
  }
  
  template <class GraphImpl>
  void SimpleTreeGraph<GraphImpl>::unRoot(bool joinRootSons)
  {
    if(joinRootSons){
      // the root must have exactly two joinRootSons
      std::vector<Node> sons = getSons(GraphImpl::getRoot());
      if(sons.size() != 2)
        throw Exception("The root must have two sons to join them.");
      GraphImpl::unlink(GraphImpl::getRoot(),sons.at(0));
      GraphImpl::unlink(GraphImpl::getRoot(),sons.at(1));
      GraphImpl::link(sons.at(0),sons.at(1));
      GraphImpl::setRoot(sons.at(0));
    }
    GraphImpl::makeUndirected();
  }
  
  template <class GraphImpl>
  std::vector<Graph::Node> SimpleTreeGraph<GraphImpl>::removeSons(Graph::Node node)
  {
    std::vector <Graph::Node> sons = getSons(node);
    for(std::vector <Graph::Node>::iterator currSon = sons.begin(); currSon != sons.end(); currSon++)
    {
      removeSon(node,*currSon);
    }
    return sons;
  }
  
  template <class GraphImpl>
  void SimpleTreeGraph<GraphImpl>::removeSon(Graph::Node node, Graph::Node son)
  {
    GraphImpl::unlink(node,son);
  }
  
  template <class GraphImpl>
  void SimpleTreeGraph<GraphImpl>::setOutGroup(Graph::Node newOutGroup)
  {
    mustBeRooted_();
    Node newRoot = createNodeFromEdge(getEdge(getFather(newOutGroup),newOutGroup));
    rootAt(newRoot);
  }
  
  template <class GraphImpl>
  std::vector<Graph::Node> SimpleTreeGraph<GraphImpl>::getNodePathBetweenTwoNodes(Graph::Node nodeA, Graph::Node nodeB, bool includeAncestor) const
  {
    GraphImpl::nodeMustExist(nodeA);
    GraphImpl::nodeMustExist(nodeB);
    std::vector<Graph::Node> path;
    std::vector<Graph::Node> pathMatrix1;
    std::vector<Graph::Node> pathMatrix2;
    
    Graph::Node nodeUp = nodeA;
    while (hasFather(nodeUp))
    {
      pathMatrix1.push_back(nodeUp);
      nodeUp = getFather(nodeUp);
    }
    pathMatrix1.push_back(nodeUp); // The root.
    
    nodeUp = nodeB;
    while (hasFather(nodeUp))
    {
      pathMatrix2.push_back(nodeUp);
      nodeUp = getFather(nodeUp);
    }
    pathMatrix2.push_back(nodeUp); // The root.
    // Must check that the two nodes have the same root!!!
    
    size_t tmp1 = pathMatrix1.size();
    size_t tmp2 = pathMatrix2.size();
    
    while ((tmp1 > 0) && (tmp2 > 0))
    {
      if (pathMatrix1[tmp1 - 1] != pathMatrix2[tmp2 - 1])
        break;
      tmp1--; tmp2--;
    }
    // (tmp1 - 1) and (tmp2 - 1) now point toward the first non-common nodes
    
    for (size_t y = 0; y < tmp1; ++y)
    {
      path.push_back(pathMatrix1[y]);
    }
    if (includeAncestor) //FIXME: one of the extremities may be the ancestor!!!
      path.push_back(pathMatrix1[tmp1]);  // pushing once, the Node that was common to both.
    for (size_t j = tmp2; j > 0; --j)
    {
      path.push_back(pathMatrix2[j - 1]);
    }
    return path;
  }

  template <class GraphImpl>
  std::vector<Graph::Edge> SimpleTreeGraph<GraphImpl>::getEdgePathBetweenTwoNodes(Graph::Node nodeA, Graph::Node nodeB) const
  {
    std::vector<Graph::Edge> path;
    std::vector<Graph::Node> pathNodes = getNodePathBetweenTwoNodes(nodeA, nodeB, true);
    for(size_t currNodeNr = 0; currNodeNr+1 < pathNodes.size(); currNodeNr++)
      path.push_back(GraphImpl::getAnyEdge(pathNodes.at(currNodeNr),pathNodes.at(currNodeNr+1)));
    return path;
  }
  
  template <class GraphImpl>
  std::vector<Graph::Node> SimpleTreeGraph<GraphImpl>::getSubtreeNodes(Graph::Node localRoot) const
  {
    mustBeValid_();
    mustBeRooted_();
    std::vector<Graph::Edge> metNodes;
    fillSubtreeMetNodes_(metNodes,localRoot);
    return metNodes;
  }
  
  template <class GraphImpl>
  std::vector<Graph::Edge> SimpleTreeGraph<GraphImpl>::getSubtreeEdges(Graph::Node localRoot) const
  {
    mustBeValid_();
    mustBeRooted_();
    std::vector<Graph::Edge> metEdges;
    fillSubtreeMetEdges_(metEdges,localRoot);
    return metEdges;
  }
  
  
  template <class GraphImpl>
  void SimpleTreeGraph<GraphImpl>::fillSubtreeMetNodes_(std::vector<Graph::Node>& metNodes, Graph::Node localRoot) const
  {
    metNodes.push_back(localRoot);
    std::vector<Graph::Node> sons = GraphImpl::getOutgoingNeighbors(localRoot);
    for(std::vector<Graph::Node>::iterator currSon = sons.begin(); currSon != sons.end(); currSon++)
      fillSubtreeMetNodes_(metNodes,*currSon);
  }
  
  template <class GraphImpl>
  void SimpleTreeGraph<GraphImpl>::fillSubtreeMetEdges_(std::vector<Graph::Edge>& metEdges, Graph::Node localRoot) const
  {
    metEdges.push_back(localRoot);
    std::vector<Graph::Edge> edgesToSons = GraphImpl::getOutgoingEdges(localRoot);
    for(std::vector<Graph::Edge>::iterator currEdgeToSon = edgesToSons.begin(); currEdgeToSon != edgesToSons.end(); currEdgeToSon++)
      fillSubtreeMetNodes_(metEdges,*currEdgeToSon);
  }
  
  
  
}



#else

namespace bpp {
  class TreeGraph;
  
  template <class GraphImpl>
  class SimpleTreeGraph;
}

#endif
