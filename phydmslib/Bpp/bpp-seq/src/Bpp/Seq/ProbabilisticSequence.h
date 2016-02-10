//
// File: ProbabilisticSequence.h
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

#ifndef _PROBABILISTIC_SEQUENCE_H_
#define _PROBABILISTIC_SEQUENCE_H_

#include "ProbabilisticSymbolList.h"
#include "CoreSequence.h"

// From the STL :
#include <string>

namespace bpp
{

/**
 * @brief The probabilisitc sequence interface.
 *
 * This is a general purpose container, containing an ordered list of
 * elements.  The states represented by the elements is defined by an
 * alphabet object.
 *
 * ProbabilisticSequence objects also contain a name attribute and
 * potentially several comment lines.
 *
 * @see Alphabet
 */
class ProbabilisticSequence :
  public virtual CoreSequence,
  public virtual ProbabilisticSymbolList
{

 public :

  /**
   * @name The Clonable interface
   *
   * @{
   */
  ProbabilisticSequence * clone() const = 0;

  /**
   * @}
   */

  // class destructor
  virtual ~ProbabilisticSequence() {}

  /**
   * @name Modifying the content of the probabilistic sequence
   *
   * @{
   */

  /**
   * @brief Set the entire content of the probabilistic sequence.
   *
   * @param content The new content of the sequence.
   * @throw Exception If the content is internally inconsistent, or is inconsistent with the specified alphabet.
   * @see The ProbabilisticSequence constructor for information about the way probabilistic sequences are internally stored.
   */
  virtual void setContent(const DataTable & content) throw (Exception) = 0;

  /**
   * @}
   */

};

/**
 * @brief A basic implementation of the ProbabilisticSequence interface.
 *
 * This is a general purpose container, containing an ordered list of
 * elmements.  The states represented by the elements is defined by an
 * alphabet object, which is passed to the contructor by a pointer.
 *
 * ProbabilisticSequence objects also contain a name attribute and
 * potentially several comment lines.
 *
 * @see Alphabet
 */

class BasicProbabilisticSequence :
  public virtual ProbabilisticSequence,
  public AbstractCoreSequence,
  public BasicProbabilisticSymbolList
{

 public :

  /**
   * @brief Empty constructor : build a void ProbabilisticSequence with just an Alphabet
   *
   * One can use it safely for all types of Alphabet in order to build
   * an empty ProbabilisticSequence, i.e., without name nor sequence
   * data.
   *
   * @param alpha A pointer to the Alphabet to be used with this ProbabilisticSequence.
   */
  BasicProbabilisticSequence(const Alphabet * alpha);

  /**
   * @brief Direct constructor : build a ProbabilisticSequence object from DataTable.
   *
   * One can use it safely for RNA, DNA and protein sequences.
   *
   * @param name The sequence name.
   * @param sequence The entire sequence to parsed as a DataTable.
   * @param alpha A pointer to the alphabet associated with this sequence.
   * @throws Exception if the content is internally inconsistent, or is inconsistent with the specified alphabet.
   */
  BasicProbabilisticSequence(const std::string & name, const DataTable & sequence, const Alphabet * alpha) throw (Exception);

  /**
   * @brief Direct constructor : build a ProbabilisticSequence object from DataTable.
   *
   * One can use it safely for RNA, DNA and protein sequences.
   *
   * @param name The sequence name.
   * @param sequence The entire sequence to parsed as a DataTable.
   * @param comments Comments to add to the sequence.
   * @param alpha A pointer to the alphabet associated with this sequence.
   * @throws Exception if the content is internally inconsistent, or is inconsistent with the specified alphabet.
   */
  BasicProbabilisticSequence(const std::string & name, const DataTable & sequence, const Comments & comments, const Alphabet * alpha) throw (Exception);

  /**
   * @brief The ProbabilisticSequence generic copy constructor.  This does not perform a hard copy of the alphabet object.
   */
  BasicProbabilisticSequence(const ProbabilisticSequence & s);

  /**
   * @brief The copy constructor.  This does not perform a hard copy of the alphabet object.
   */
  BasicProbabilisticSequence(const BasicProbabilisticSequence & s);

  /**
   * @brief The ProbabilisticSequence generic assignment operator.  This does not perform a hard copy of the alphabet object.
   *
   * @return a reference to the assigned BasicProbabilisticSequence.
   */
  BasicProbabilisticSequence & operator=(const ProbabilisticSequence & s);

  /**
   * @brief The assignment operator.  This does not perform a hard cop of the alphabet object.
   *
   * @return A reference to the assigned BasicProbabilisticSequence.
   */
  BasicProbabilisticSequence & operator=(const BasicProbabilisticSequence & s);

  /**
   * @name The Clonable interface
   *
   * @{
   */
  BasicProbabilisticSequence * clone() const { return new BasicProbabilisticSequence(*this); }

  /**
   * @}
   */

  // class destructor
  virtual ~BasicProbabilisticSequence() {}

 public :

  void setContent(const DataTable & content) throw (Exception) { BasicProbabilisticSymbolList::setContent(content); }

  void setToSizeR(size_t newSize) {}; // leave empty for now

  void setToSizeL(size_t newSize) {};

};

} // end of namespace bpp

#endif // _PROBABILISTIC_SEQUENCE_H_
