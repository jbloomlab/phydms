//
// File: CoreSite.h
// Created by: Murray Patterson
//             Julien Dutheil
// Created on: Mon Oct 12 2015
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

#ifndef _CORESITE_H_
#define _CORESITE_H_

namespace bpp
{

/**
 * @brief The core site interface. 
 *
 * The core interface for sites manipulation.  It is very similar to
 * the CoreSequence interface (a site is a vertical sequence!).  Sites
 * have a 'position' attribute.
 * This attribute stands for an index in an alignment, and may be
 * used as a unique identifier, in the same manner that names identify
 * sequence objects.
 */
class CoreSite :
  public virtual Clonable
{

 public :

  /**
   * @name The Clonable interface
   *
   * @{
   */
  CoreSite * clone() const = 0;

  /**
   * @}
   */

  // class destructor
  virtual ~CoreSite() {}

  /**
   * @name Setting/getting the position of the site.
   *
   * @{
   */

  /**
   * @brief Get the position of this site.
   *
   * @return The position of this site.
   */
  virtual int getPosition() const = 0;

  /**
   * @brief Set the position of this site.
   *
   * @param position The new position of the site.
   */
  virtual void setPosition(int position) = 0;

  /**
   * @}
   */

};

/**
 * @brief A partial implementation of the CoreSite interface. 
 */
class AbstractCoreSite :
  public virtual CoreSite
{

 private :

  /**
   * @brief The position associated with this site.
   */
  int position_;

 public :

  /**
   * @brief Constructor of the AbstractCoreSite object.
   *
   * Construct an 'empty' object, i.e., with no position assocated.
   */
  AbstractCoreSite() :
    position_(0) {}

  /**
   * @brief Constructor of the AbstractCoreSite object.
   *
   * @param position The position of the site.
   */
  AbstractCoreSite(int position) :
    position_(position) {}

  /**
   * @name The copy constructors.
   *
   * @{
   */
  AbstractCoreSite(const CoreSite & site) :
    position_(site.getPosition()) {}

  AbstractCoreSite(const AbstractCoreSite & site) :
    position_(site.position_) {}

  /**
   * @}
   */

  /**
   * @name The assignment operators.
   *
   * @{
   */
  AbstractCoreSite & operator=(const CoreSite & site) {
    position_ = site.getPosition();
    return *this;
  }

  AbstractCoreSite & operator=(const AbstractCoreSite & site) {
    position_ = site.position_;
    return *this;
  }

  /**
   * @}
   */

  /**
   * @name The Clonable interface
   */
  AbstractCoreSite * clone() const = 0;

  /**
   * @}
   */

  // class destructor
  virtual ~AbstractCoreSite() {}

 public :

  virtual int getPosition() const { return position_; }

  virtual void setPosition(int position) { position_ = position; }
  
};

} //end of namespace bpp.

#endif // _CORESITE_H_
