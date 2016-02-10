//
// File: ProbabilisticSymbolList.h
// Created by: Murray Patterson
// Created on: Sun Oct 4 2015
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

#ifndef _PROBABILISTIC_SYMBOLLIST_H_
#define _PROBABILISTIC_SYMBOLLIST_H_

#include "Alphabet/Alphabet.h"
#include <Bpp/Clonable.h>
#include <Bpp/Numeric/DataTable.h>

// From the STL :
#include <string>
#include <vector>

namespace bpp
{

/**
 * @brief The ProbabilisticSymobolList interface.
 *
 * @see Alphabet
 */
class ProbabilisticSymbolList :
  public virtual Clonable
{

 public :

  /**
   * @name The Clonable interface
   *
   * @{
   */
  ProbabilisticSymbolList * clone () const = 0;

  /**
   * @}
   */

  // class destructor
  virtual ~ProbabilisticSymbolList() {}

  /**
   * @brief Get the alphabet associated to with the list.
   *
   * @return A const pointer to the alphabet.
   * @see Alphabet class.
   */
  virtual const Alphabet * getAlphabet() const = 0;

  /**
   * @brief Get the number of elements in the list.
   *
   * @return The number of sites in the list.
   */
  virtual size_t size() const = 0;

  /**
   * @name Modifying the content of the list.
   *
   * @{
   */

  /**
   * @brief Set the entire content of the list.
   *
   * @param list The new content of the list.
   * @throw Exception If the content is internally inconsistent, or is inconsistent with the specified alphabet.
   * @see The ProbabilisticSymbolList constructor for information about the way lists are internally stored.
   */
  virtual void setContent(const DataTable & list) throw (Exception) = 0;

  /**
   * @brief Add an element to the end of the list.
   *
   * @param e The elment to add, given as a vector of string.
   * @throw Exception If the element is internally inconsistent, or is inconsistent with the specified alphabet.
   */
  virtual void addElement(const std::vector<std::string> & element) throw (Exception) = 0;

  /**
   * @brief Delete the element at position 'pos'.
   *
   * @param pos the Position of the element to delete.
   * @throw IndexOutOfBoundsException if position is not in the list.
   */
  virtual void deleteElement(size_t pos) throw (IndexOutOfBoundsException) = 0;

  /**
   * @}
   */
  
  /**
   * @name Retrieving the content of a list.
   *
   * @{
   */

  /**
   * @brief Get the entire content of the list.
   *
   * @return list The content of the list.
   */
  virtual const DataTable & getContent() const = 0;

};

/**
 * @brief A basic ProbabilisticSymbolList object.
 *
 * This is a general purpose container, containing an ordered list of
 * elements.  The states represented by the elements is defined by an
 * alphabet object, which is passed to the list constructor by a
 * pointer.
 *
 * @see Alphabet
 */
class BasicProbabilisticSymbolList :
  public virtual ProbabilisticSymbolList
{

 private :

  /**
   * @brief The Alphabet attribute must be initialized in the constructor and then can never be changed.
   *
   * To apply another alphabet to the list requires creating another
   * list.
   */
  const Alphabet * alphabet_;

 protected :

  /**
   * @brief The list content.
   */
  DataTable content_;

 public :

  /**
   * @brief Build a new void BasicProbabilisticSymbolList object with the specified alphabet.
   *
   * @param alpha the alphabet to use.
   */
  BasicProbabilisticSymbolList(const Alphabet * alpha);

  /**
   * @brief Build a new BasicProbabilisticSymbolList object with the specified alphabet.
   *
   * @param list The content of the site.
   * @param alpha The alphabet to use.
   * @throw If the content is internally inconsistent, or is inconsistent with the specified alphabet.
   */
  BasicProbabilisticSymbolList(const DataTable & list, const Alphabet * alpha) throw (Exception);

  /**
   * @brief The generic copy constructor.
   */
  BasicProbabilisticSymbolList(const ProbabilisticSymbolList & list);

  /**
   * @brief The copy constructor.
   */
  BasicProbabilisticSymbolList(const BasicProbabilisticSymbolList & list);

  /**
   * @brief The generic assignment operator.
   */
  BasicProbabilisticSymbolList & operator=(const ProbabilisticSymbolList & list);

  /**
   * @brief The assignement operator.
   */
  BasicProbabilisticSymbolList & operator=(const BasicProbabilisticSymbolList & list);

  /**
   * @name The Clonable interface
   *
   * @{
   */
  BasicProbabilisticSymbolList * clone() const { return new BasicProbabilisticSymbolList(* this); }

  /**
   * @}
   */

  // class destructor
  virtual ~BasicProbabilisticSymbolList() {}

 public :

  virtual const Alphabet * getAlphabet() const { return alphabet_; }

  virtual size_t size() const { return static_cast<size_t>(content_.getNumberOfRows()); }

  virtual void setContent(const DataTable & list) throw (Exception);

  virtual void addElement(const std::vector<std::string> & element) throw (Exception);

  virtual void deleteElement(size_t pos) throw (IndexOutOfBoundsException) { content_.deleteRow(pos); }

  virtual const DataTable & getContent() const { return content_; }

};

} // end of namespace bpp

#endif // _PROBABILISTIC_SYMBOLLIST_H
