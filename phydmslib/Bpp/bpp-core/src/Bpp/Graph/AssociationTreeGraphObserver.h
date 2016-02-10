#ifndef _ASSOCIATIONTREEGRAPHOBSERVER_HPP_
#define _ASSOCIATIONTREEGRAPHOBSERVER_HPP_

#include "TreeGraph.h"
#include "AssociationGraphObserver.h"

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
    class AssociationTreeGraphObserver:
    public virtual AssociationGraphObserver<N,E>
    {
      
    
    };
    
    
    template <class N, class E, class TreeGraphImpl>
    class SimpleAssociationTreeGraphObserver:
    public virtual AssociationTreeGraphObserver<N,E>,
    public virtual SimpleAssociationGraphObserver<N,E,TreeGraphImpl>
    {
    public:
      typedef typename AssociationGraphObserver<N,E>::NodeIndex NodeIndex;
      typedef typename AssociationGraphObserver<N,E>::EdgeIndex EdgeIndex;
      
      typedef typename Graph::Node NodeGraphid;
      typedef typename Graph::Edge EdgeGraphid;

    private:
      TreeGraphImpl* subjectTreeGraph_;
      

    public:
      
      /**
      * Constructor
      * @param directed is the graph directed
      */
      SimpleAssociationTreeGraphObserver(bool rooted = false);
      
      /**
      * Constructor
      * @param subjectGraph the graph which is observed
      */
      SimpleAssociationTreeGraphObserver(TreeGraphImpl* subjectTreeGraph = 00);
      
      /**
      * Copy Constructor
      * @param graphObserver the treeGraphObserver to be copied
      */
      SimpleAssociationTreeGraphObserver(bpp::SimpleAssociationTreeGraphObserver<N,E,TreeGraphImpl> const& treeGraphObserver);
      
      /**
      * = Operator
      * @param graphObserver the treeGraphObserver we want to copy the values
      */
      SimpleAssociationTreeGraphObserver<N,E,TreeGraphImpl> operator=(bpp::SimpleAssociationTreeGraphObserver<N,E,TreeGraphImpl> const& treeGraphObserver);
      
      /**
      * Destructor
      * @param graphObserver the treeGraphObserver to be copied
      */
      ~SimpleAssociationTreeGraphObserver();
      
      
      
      /**
      * clone function
      * @param graphObserver the graphObserver to be copied
      */
      #ifdef NO_VIRTUAL_COV
        Clonable*
      #else
        SimpleAssociationTreeGraphObserver<N,E,TreeGraphImpl>*
      #endif
      clone() const { return new SimpleAssociationTreeGraphObserver<N,E,TreeGraphImpl>(*this); };
      
      
    /**
     * Is the graph a tree? A tree must be acyclic and with no isolated node.
     * @return true if valid tree
     */
     bool isValid() const;
     
     
     /**
     * Return, in a rooted tree, the branch leading to the father
     * @param nodeObject the concerned node
     * @return an Edge which is the branch to the father
     */
     E* getBranchToFather(const N* nodeObject) const;
      
     /**
      * Is the subject tree rooted?
      */
     bool isRooted() const;
     
     
     /**
     * Return, in a rooted tree, the father node
     * @param nodeObject the concerned node
     * @return the father
     */
     N* getFather(const N* nodeObject) const;
     
     /**
      * Has the node a father?
      */
     bool hasFather(const N* nodeObject) const;
     
     
    /**
     * Return, in a rooted tree, the branch leading to the father
     * @param nodeObject the concerned node
     * @return an Edge which is the branch to the father
     */
     std::vector<N*> getSons(const N* node) const;
     std::vector<NodeIndex> getSons(const NodeIndex node) const;
     
     /**
     * Return, in a rooted tree, the number of sons
     * @param nodeObject the concerned node
     * @return the number of sons
     */
     unsigned int getNumberOfSons(const N* node) const;
     unsigned int getNumberOfSons(NodeIndex node) const;
     
     /**
     * Return, in a rooted tree, the number of leaves under a certain node
     * @param nodeObject the concerned node
     * @return the number of leaves
     */
     unsigned int getNumberOfLeaves(const N* node) const;
     unsigned int getNumberOfLeaves(NodeIndex node) const;
     
     /**
     * Remove the son of a node
     * @return a vector containing the removed nodes
     */
     std::vector<N*> removeSons(N* node);
     
     
     /**
     * Change / set the father of a node
     * @param nodeObject the concerned node
     * @param fatherNodeObject the node to be the father
     */
     void setFather(const N* nodeObject, const N* fatherNodeObject);
     
     /**
     * Add a son to a node
     * @param nodeObject the concerned node
     * @param sonNodeObject the node to be added as a son to the father
     */
     void addSon(const N* nodeObject, const N* sonNodeObject);
     
     
     /// FROM TREETOOLS
     
     
    /**
     * @brief Get a vector of ancestor nodes between to nodes.
     *
     * @param nodeId1 first node.
     * @param nodeId2 second node.
     * @param includeAncestor Tell if the common ancestor must be included in the vector.
     * @return A vector of ancestor nodes ids.
     * @throw PhyloNodeNotFoundException If a node is not found.
     */
     std::vector<N*> getNodePathBetweenTwoNodes(const N* nodeObjectA, const N* nodeObjectB, bool includeAncestor = true) const;
     
     std::vector<E*> getEdgePathBetweenTwoNodes(const N* nodeObjectA, const N* nodeObjectB) const;
     
     std::vector<N*> getSubtreeNodes(const N* localRoot);
     std::vector<E*> getSubtreeEdges(const N* localRoot);
     
    };
    
    
    
    template <class N, class E, class TreeGraphImpl>
    SimpleAssociationTreeGraphObserver<N,E,TreeGraphImpl>::SimpleAssociationTreeGraphObserver(bool rooted_p):
    SimpleAssociationGraphObserver<N,E,TreeGraphImpl>::SimpleAssociationGraphObserver(rooted_p),
    subjectTreeGraph_(00)
    {
      subjectTreeGraph_ = static_cast<TreeGraphImpl*>(SimpleAssociationGraphObserver<N,E,TreeGraphImpl>::getGraph());
    }
    
    
    template <class N, class E, class TreeGraphImpl>
    SimpleAssociationTreeGraphObserver<N,E,TreeGraphImpl>::SimpleAssociationTreeGraphObserver(SimpleAssociationTreeGraphObserver<N,E,TreeGraphImpl> const& treeGraphObserver):
    SimpleAssociationGraphObserver<N,E,TreeGraphImpl>::SimpleAssociationGraphObserver(treeGraphObserver),
    subjectTreeGraph_(00)
    {
    }
    
    
    template <class N, class E, class TreeGraphImpl>
    SimpleAssociationTreeGraphObserver<N,E,TreeGraphImpl> SimpleAssociationTreeGraphObserver<N,E,TreeGraphImpl>::operator=(SimpleAssociationTreeGraphObserver<N,E,TreeGraphImpl> const& treeGraphObserver)
    {
      SimpleAssociationGraphObserver<N,E,TreeGraphImpl>::operator=(treeGraphObserver);
    }
    
    template <class N, class E, class TreeGraphImpl>
    SimpleAssociationTreeGraphObserver<N,E,TreeGraphImpl>::~SimpleAssociationTreeGraphObserver()
    {
      //delete subjectTreeGraph_;
    }
    
    
    template <class N, class E, class TreeGraphImpl>
    bool SimpleAssociationTreeGraphObserver<N,E,TreeGraphImpl>::isValid() const
    {
      return subjectTreeGraph_->isValid();
    }
    
    template <class N, class E, class TreeGraphImpl>
    bool SimpleAssociationTreeGraphObserver<N,E,TreeGraphImpl>::isRooted() const
    {
      return subjectTreeGraph_->isRooted();
    }
    
    template <class N, class E, class TreeGraphImpl>
    bool SimpleAssociationTreeGraphObserver<N,E,TreeGraphImpl>::hasFather(const N* nodeObject) const
    {
      return subjectTreeGraph_->hasFather(getNodeGraphid(nodeObject));
    }
    
    template <class N, class E, class TreeGraphImpl>
    N* SimpleAssociationTreeGraphObserver<N,E,TreeGraphImpl>::getFather(const N* nodeObject) const
    {
      return getNodeFromGraphid(subjectTreeGraph_->getFather(getNodeGraphid(nodeObject)));
    }
    
    template <class N, class E, class TreeGraphImpl>
    E* SimpleAssociationTreeGraphObserver<N,E,TreeGraphImpl>::getBranchToFather(const N* nodeObject) const
    {
      return(getEdgeObject(subjectTreeGraph_->getBranchToFather(getNodeGraphid(nodeObject))));
    }
    
    template <class N, class E, class TreeGraphImpl>
    std::vector<N*> SimpleAssociationTreeGraphObserver<N,E,TreeGraphImpl>::getSons(const N* nodeObject) const
    {
      return getNodesFromGraphid(subjectTreeGraph_->getSons(getNodeGraphid(nodeObject)));
    }
    
    template <class N, class E, class TreeGraphImpl>
    std::vector<typename SimpleAssociationTreeGraphObserver<N,E,TreeGraphImpl>::NodeIndex> SimpleAssociationTreeGraphObserver<N,E,TreeGraphImpl>::getSons(SimpleAssociationTreeGraphObserver<N,E,TreeGraphImpl>::NodeIndex node) const
    {
      return getNodeIndexes(getNodesFromGraphid(subjectTreeGraph_->getSons(getNodeGraphid(getNode(node)))));
    }
    
    template <class N, class E, class TreeGraphImpl>
    unsigned int SimpleAssociationTreeGraphObserver<N,E,TreeGraphImpl>::getNumberOfSons(const N* nodeObject) const
    {
      return getNodesFromGraphid(subjectTreeGraph_->getSons(getNodeGraphid(nodeObject))).size();
    }
    
    template <class N, class E, class TreeGraphImpl>
    unsigned int SimpleAssociationTreeGraphObserver<N,E,TreeGraphImpl>::getNumberOfLeaves(const N* nodeObject) const
    {
      return getNodesFromGraphid(subjectTreeGraph_->getLeavesFromNode(getNodeGraphid(nodeObject))).size();
    }
    
    template <class N, class E, class TreeGraphImpl>
    void SimpleAssociationTreeGraphObserver<N,E,TreeGraphImpl>::setFather(const N* nodeObject, const N* fatherNodeObject)
    {
      subjectTreeGraph_->setFather(getNodeGraphid(nodeObject),getNodeGraphid(fatherNodeObject));
    }
    
    template <class N, class E, class TreeGraphImpl>
    void SimpleAssociationTreeGraphObserver<N,E,TreeGraphImpl>::addSon(const N* nodeObject, const N* sonNodeObject)
    {
      subjectTreeGraph_->addSon(getNodeGraphid(nodeObject),getNodeGraphid(sonNodeObject));
    }
    
    template <class N, class E, class TreeGraphImpl>
    std::vector<N*> SimpleAssociationTreeGraphObserver<N,E,TreeGraphImpl>::removeSons(N* const node)
    {
      return getNodesFromGraphid(subjectTreeGraph_->removeSons(getNodeGraphid(node)));
    }
    
    template <class N, class E, class TreeGraphImpl>
    std::vector<N*> SimpleAssociationTreeGraphObserver<N,E,TreeGraphImpl>::getNodePathBetweenTwoNodes(const N* nodeA, const N* nodeB, bool includeAncestor) const
    {
      return getNodesFromGraphid(subjectTreeGraph_->getNodePathBetweenTwoNodes(getNodeGraphid(nodeA),getNodeGraphid(nodeB),includeAncestor));
    }
    
    template <class N, class E, class TreeGraphImpl>
    std::vector<E*> SimpleAssociationTreeGraphObserver<N,E,TreeGraphImpl>::getEdgePathBetweenTwoNodes(const N* nodeA, const N* nodeB) const
    {
      return getEdges(subjectTreeGraph_->getEdgePathBetweenTwoNodes(getNodeGraphid(nodeA),getNodeGraphid(nodeB)));
    }
    
    template <class N, class E, class TreeGraphImpl>
    std::vector<N*> SimpleAssociationTreeGraphObserver<N,E,TreeGraphImpl>::getSubtreeNodes(const N* localRoot)
    {
      return getNodesFromGraphid(subjectTreeGraph_->getSubtreeNodes(getNodeGraphid(localRoot)));
    }
    
    template <class N, class E, class TreeGraphImpl>
    std::vector<E*> SimpleAssociationTreeGraphObserver<N,E,TreeGraphImpl>::getSubtreeEdges(const N* localRoot)
    {
      return getEdges(subjectTreeGraph_->getSubtreeNodes(getNodeGraphid(localRoot)));
    }
 
 
}

#else

namespace bpp {class AssociationTreeGraphObserver;}

#endif