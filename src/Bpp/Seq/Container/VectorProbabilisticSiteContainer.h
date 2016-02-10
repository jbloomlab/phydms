//
// File: VectorProbabilisticSiteContainer.h
// Created by: Murray Patterson
// Created on: Mon Oct 19 2015
//

/*
   Copyright or © or Copr. CNRS, (November 17, 2004)

   This software is a computer program whose purpose is to provide classes
   for sequences analysis.

   This software is governed by the CeCILL  license under French law and
   abiding by the rules of distribution of free software.  You can  use,
   modify and/ or redistribute the software under the terms of the CeCILL
   license as circulated by CEA, CNRS and INRIA at the following URL
   "http://www.cecill.info".

   As a counterpart to the access to the source code and  rights to copy,
   modify and redistribute granted by the license, users are provided only
   with a limited warranty  and the software's author,  the holder of the
   economic rights,  and the successive licensors  have only  limited
   liability.

   In this respect, the user's attention is drawn to the risks associated
   with loading,  using,  modifying and/or developing or reproducing the
   software by the user in light of its specific status of free software,
   that may mean  that it is complicated to manipulate,  and  that  also
   therefore means  that it is reserved for developers  and  experienced
   professionals having in-depth computer knowledge. Users are therefore
   encouraged to load and test the software's suitability as regards their
   requirements in conditions enabling the security of their systems and/or
   data to be ensured and,  more generally, to use and operate it in the
   same conditions as regards security.

   The fact that you are presently reading this means that you have had
   knowledge of the CeCILL license and that you accept its terms.
*/

#ifndef _VECTORPROBABILISTICSITECONTAINER_H_
#define _VECTORPROBABILISTICSITECONTAINER_H_

#include "VectorSiteContainer.h"

#include "../ProbabilisticSite.h"
#include "../ProbabilisticSequence.h"

// From the STL :
#include <vector>

namespace bpp
{

/**
 * @brief The VectorProbabilisticSiteContainer class.
 *
 * ProbabilisticSites are stored in a std::vector of pointers.
 * ProbabilisticSite access is hence in \f$O(1)\f$, and sequence
 * access in \f$O(l)\f$, where \f$l\f$ is the number of sites in the
 * container.
 *
 * This is a modified copy of VectorSiteContainer with the minimum
 * changes necessary for it to work with bppML and bppAncestor
 *
 * @see ProbabilisticSequence, ProbabilisticSite, VectorSiteContainer
 */
class VectorProbabilisticSiteContainer :
  public VectorSiteContainer
{

 protected :

  std::vector<ProbabilisticSite *> p_sites_;
  mutable std::vector<ProbabilisticSequence *> p_sequences_;

 public :

  /**
   * @brief Build a new emtpy container.
   *
   * @param alpha The alphabet for this container.
   */
  VectorProbabilisticSiteContainer(const Alphabet * alpha);

  /**
   * @name The Clonable interface.
   *
   * @{
   */
  VectorSiteContainer* clone() const { return new VectorSiteContainer(*this); }

  /**
   * @}
   */

  // class destructor
  virtual ~VectorProbabilisticSiteContainer() { pClear(); }

 public :

  const ProbabilisticSite & getProbabilisticSite(std::size_t i) const throw (IndexOutOfBoundsException);

  void addSite(const ProbabilisticSite & site, bool checkPosition = true) throw (Exception);

  const ProbabilisticSequence & getProbabilisticSequence(std::size_t i) const throw (IndexOutOfBoundsException);

  void addSequence(const ProbabilisticSequence & sequence, bool checkName = true) throw (Exception);

  std::size_t getNumberOfProbabilisticSites() const { return p_sites_.size(); }
  std::size_t getNumberOfProbabilisticSequences() const { return p_sequences_.size(); }

  void pClear();
  void reindexpSites();

 protected :

  // create n void probabilistic sites :
  void pRealloc(std::size_t n);
};

} // end of namespace bpp

#endif  // _VECTORPROBABILISTICSITECONTAINER_H_
