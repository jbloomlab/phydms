#ifndef _GRAPHOBSERVER_HPP_
#define _GRAPHOBSERVER_HPP_

#include "Graph.h"

#include "../Exceptions.h"
#include "../Clonable.h"

#include <vector>
#include <map>
#include <iostream>
#include <ostream>


namespace bpp
{

  /**
   * @brief Defines a Graph Observer. It is a template which follows
   * (subscribed to) a Graph.
   * The graph and the graph observer communicate to keep them up-to-date
   * each other. The observer is also an actor, since it can change
   * the structure of the observed Graph.
   *
   * @author Thomas Bigot
   */
  
    // interface
    
    class GraphObserver:
    public virtual bpp::Clonable
    {
    
    public:
            
      /** @name Function called by the subjectGraph
      *  Methods called by the subject graph to make this observer so fit the subject graph
      */
      ///@{
      
      /**
      * Delete unused object edges, since they have been deleted in the graph
      * @param edgesToDelete a vector of Edges to delete
      */
      virtual void deletedEdgesUpdate(std::vector< bpp::Graph::Edge >& edgesToDelete) = 0;
      
      /**
      * Delete unused object nodes, since they have been deleted in the graph
      * @param nodesToDelete a vector of N to delete
      */
      virtual void deletedNodesUpdate(std::vector< bpp::Graph::Node >& nodesToDelete) = 0;
      
      ///@}
      
    };

}

#else

namespace bpp {class GraphObserver;}

#endif