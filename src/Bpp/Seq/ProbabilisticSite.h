//
// File ProbabilisticSite.h
// Created by: Murray Patterson
// Created on: Tue Oct 13 2015
//

/*
Copyright or © or Copr. Bio++ Development Team, (November 17, 2004)

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

#ifndef _PROBABILISTIC_SITE_H_
#define _PROBABILISTIC_SITE_H_

#include "ProbabilisticSymbolList.h"
#include "CoreSite.h"

namespace bpp
{

/**
 * @brief The ProbabilisticSite interface.
 *
 * Define specific attributes and methods for probabilistic sites
 * manipulation.  It is very similar to the ProbabilisticSequence
 * object (a site is a vertical probabilistic sequence!),
 * ProbabilisticSites have a 'position' attribute.
 *
 * This attribute stands for an index in an alignment, and may be used
 * as a unique identifier, in the same manner that names identify
 * probabilistic sequence objects.  But for now, we do not allow to
 * construct a ProbabilisticSite directly from a string.  This should
 * not be a constraint, since you never read probabilistic sites
 * directly from a file.
 */
class ProbabilisticSite :
  public virtual CoreSite,
  public virtual ProbabilisticSymbolList 
{

 public :

  /**
   * @name The Colonable interface
   *
   * @{
   */
  ProbabilisticSite * clone() const = 0;

  /**
   * @}
   */

  // class destructor
  virtual ~ProbabilisticSite() {}

  /**
   * @name Setting/getting the position of the probabilistic site.
   *
   * @{
   */

  /**
   * @brief Get the position of this probabilistic site.
   *
   * @return The position of this probabilistic site.
   */
  virtual int getPosition() const = 0;

  /**
   * @brief Set the position of this probabilistic site.
   *
   * @param position The new position of the probabilistic site.
   */
  virtual void setPosition(int position) = 0;

};

/**
 * @brief A basic implementation of the ProbabilisticSite interface.
 *
 * Define specific attributes and methods for probabilistic sites
 * manipulation.  It is very similar to the ProbabilisticSequence
 * object (a site is a vertical probabilistic sequence!),
 * ProbabilisticSites have a 'position' attribute.
 *
 * This attribute stands for an index in an alignment, and may be used
 * as a unique identifier, in the same manner that names identify
 * probabilistic sequence objects.  But for now, we do not allow to
 * construct a ProbabilisticSite directly from a string.  This should
 * not be a constraint, since you never read probabilistic sites
 * directly from a file.
 */
class BasicProbabilisticSite :
  public virtual ProbabilisticSite,
  public AbstractCoreSite,
  public BasicProbabilisticSymbolList
{

 public:
    
  /**
   * @brief Build a new void BasicProbabilisticSite object with the specified alphabet.
   *
   * @param alpha The alphabet to use.
   */
  BasicProbabilisticSite(const Alphabet* alpha);

  /**
   * @brief Build a new void BasicProbabilisticSite object with the specified alphabet and position.
   *
   * @param alpha    The alphabet to use.
   * @param position The position attribute for this probabilistic site.
   */
  BasicProbabilisticSite(const Alphabet* alpha, int position);

  /**
   * @brief Build a new BasicProbabilisticSite object with the specified alphabet.
   * The content of the site is initialized from a DataTable object.
   *
   * @param site     The content of the probabilist site.
   * @param alpha    The alphabet to use.
   * @throw Exception If the content does not match the specified alphabet or is internally inconsistent.
   */
  BasicProbabilisticSite(const DataTable & site, const Alphabet* alpha) throw (Exception);

  /**
   * @brief Build a new BasicProbabilisticSite object with the specified alphabet and position.
   * The content of the probabilistic site is initialized from a DataTable object.
   *
   * @param site     The content of the probabilistic site.
   * @param alpha    The alphabet to use.
   * @param position The position attribute for this site.
   * @throw Exception If the content does not match the specified alphabet or is internally insconsistent.
   */
  BasicProbabilisticSite(const DataTable & site, const Alphabet* alpha, int position) throw (Exception);

  /**
   * @brief The generic copy constructor.
   */
  BasicProbabilisticSite(const ProbabilisticSite & site);

  /**
   * @brief The copy constructor.
   */
  BasicProbabilisticSite(const BasicProbabilisticSite & site);

  /**
   * @brief The generic assignment operator.
   */
  BasicProbabilisticSite & operator=(const ProbabilisticSite & site);

  /**
   * @brief The assignment operator.
   */
  BasicProbabilisticSite & operator=(const BasicProbabilisticSite & site);

  /**
   * @name The Clonable interface
   *
   * @{
   */
  BasicProbabilisticSite * clone() const { return new BasicProbabilisticSite(*this); }

  /**
   * @}
   */

  // class destructor
  virtual ~BasicProbabilisticSite() {}

 public :

  int getPosition() const { return AbstractCoreSite::getPosition(); }
  
  void setPosition(int position) { AbstractCoreSite::setPosition(position); }

};

} //end of namespace bpp.

#endif  //_PROBABILISTIC_SITE_H_
