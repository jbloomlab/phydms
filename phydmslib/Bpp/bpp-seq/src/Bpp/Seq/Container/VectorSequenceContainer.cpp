//
// File VectorSequenceContainer.cpp
// Author : Guillaume Deuchst
//          Julien Dutheil
// Last modification : Wednesday July 30 2003
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

#include "VectorSequenceContainer.h"
#include <Bpp/Text/TextTools.h>

using namespace bpp;
using namespace std;

/** Class constructors: *******************************************************/

VectorSequenceContainer::VectorSequenceContainer(
  const std::vector<const Sequence*>& vs,
  const Alphabet* alpha)
throw (AlphabetMismatchException) :
  AbstractSequenceContainer(alpha),
  sequences_()
{
  for (std::vector<const Sequence*>::const_iterator i = vs.begin(); i < vs.end(); i++)
  {
    addSequence(**i);
  }
}

/** Copy constructors: ********************************************************/

VectorSequenceContainer::VectorSequenceContainer(
  const VectorSequenceContainer& vsc) :
  AbstractSequenceContainer(vsc),
  sequences_()
{
  size_t max = vsc.getNumberOfSequences();
  for (size_t i = 0; i < max; i++)
  {
    addSequence(vsc.getSequence(i), false);
  }
}

VectorSequenceContainer::VectorSequenceContainer(
  const OrderedSequenceContainer& osc) :
  AbstractSequenceContainer(osc.getAlphabet()),
  sequences_()
{
  // Sequences insertion
  for (unsigned int i = 0; i < osc.getNumberOfSequences(); i++)
  {
    addSequence(osc.getSequence(i), false);
  }
}

VectorSequenceContainer::VectorSequenceContainer(
  const SequenceContainer& sc) :
  AbstractSequenceContainer(sc.getAlphabet()),
  sequences_()
{
  // Sequences insertion
  std::vector<std::string> names = sc.getSequencesNames();
  for (unsigned int i = 0; i < names.size(); i++)
  {
    addSequence(sc.getSequence(names[i]), false);
  }

  setGeneralComments(sc.getGeneralComments());
}

/** Assignation operator: *****************************************************/

VectorSequenceContainer& VectorSequenceContainer::operator=(
  const VectorSequenceContainer& vsc)
{
  clear();
  AbstractSequenceContainer::operator=(vsc);

  // Sequences insertion
  size_t max = vsc.getNumberOfSequences();
  for (size_t i = 0; i < max; i++)
  {
    addSequence(vsc.getSequence(i), false);
  }

  return *this;
}

VectorSequenceContainer& VectorSequenceContainer::operator=(
  const OrderedSequenceContainer& osc)
{
  clear();
  AbstractSequenceContainer::operator=(osc);

  // Sequences insertion
  size_t max = osc.getNumberOfSequences();
  for (unsigned int i = 0; i < max; i++)
  {
    addSequence(osc.getSequence(i), false);
  }

  return *this;
}

/******************************************************************************/

VectorSequenceContainer& VectorSequenceContainer::operator=(
  const SequenceContainer& sc)
{
  clear();
  AbstractSequenceContainer::operator=(sc);

  // Seq names:
  std::vector<std::string> names = sc.getSequencesNames();

  for (unsigned int i = 0; i < names.size(); i++)
  {
    addSequence(sc.getSequence(names[i]), false);
  }

  return *this;
}

/******************************************************************************/

const Sequence& VectorSequenceContainer::getSequence(size_t sequenceIndex) const throw (IndexOutOfBoundsException)
{
  // Specified sequence existence verification
  if (sequenceIndex < sequences_.size())
    return *sequences_[sequenceIndex];
  throw IndexOutOfBoundsException("VectorSequenceContainer::getSequence.", sequenceIndex, 0, sequences_.size() - 1);
}

/******************************************************************************/

bool VectorSequenceContainer::hasSequence(const string& name) const
{
  // Specified sequence name research into all sequences
  for (size_t i = 0; i < sequences_.size(); i++)
  {
    if (sequences_[i]->getName() == name)
      return true;
  }
  return false;
}

/******************************************************************************/

const Sequence& VectorSequenceContainer::getSequence(const string& name) const throw (SequenceNotFoundException)
{
  // Specified sequence name research into all sequences
  for (size_t i = 0; i < sequences_.size(); i++)
  {
    if (sequences_[i]->getName() == name)
      return *sequences_[i];
  }
  throw SequenceNotFoundException("VectorSequenceContainer::getSequence : Specified sequence doesn't exist", name);
}

/******************************************************************************/

Sequence& VectorSequenceContainer::getSequence_(size_t sequenceIndex) throw (IndexOutOfBoundsException)
{
  // Specified sequence existence verification
  if (sequenceIndex < sequences_.size())
    return *sequences_[sequenceIndex];
  throw IndexOutOfBoundsException("VectorSequenceContainer::getSequence.", sequenceIndex, 0, sequences_.size() - 1);
}

/******************************************************************************/

Sequence& VectorSequenceContainer::getSequence_(const string& name) throw (SequenceNotFoundException)
{
  // Specified sequence name research into all sequences
  for (size_t i = 0; i < sequences_.size(); i++)
  {
    if (sequences_[i]->getName() == name)
      return *sequences_[i];
  }
  throw SequenceNotFoundException("VectorSequenceContainer::getSequence : Specified sequence doesn't exist", name);
}

/******************************************************************************/

size_t VectorSequenceContainer::getSequencePosition(const string& name) const throw (SequenceNotFoundException)
{
  // Specified sequence name research into all sequences
  for (size_t i = 0; i < sequences_.size(); i++)
  {
    if (sequences_[i]->getName() == name)
      return i;
  }
  throw SequenceNotFoundException("VectorSequenceContainer::getSequencePosition : Specified sequence doesn't exist", name);
}

/******************************************************************************/

void VectorSequenceContainer::setSequence(size_t sequenceIndex, const Sequence& sequence, bool checkName) throw (Exception)
{
  // Sequence's name existence checking
  if (checkName)
  {
    // For all names in vector : throw exception if name already exists
    for (size_t j = 0; j < sequences_.size(); j++)
    {
      if (sequences_[j]->getName() == sequence.getName())
        if (j != sequenceIndex)
          throw Exception("VectorSequenceContainer::setSequence : Sequence's name already exists in container");
    }
  }

  // New sequence's alphabet and sequence container's alphabet matching verification
  if (sequence.getAlphabet()->getAlphabetType() == getAlphabet()->getAlphabetType())
  {
    // Delete old sequence
    delete sequences_[sequenceIndex];
    // New sequence insertion in sequence container
    sequences_[sequenceIndex] = dynamic_cast<Sequence*>(sequence.clone());
  }
  else
    throw AlphabetMismatchException("VectorSequenceContainer::setSequence : Alphabets don't match", getAlphabet(), sequence.getAlphabet());
}

/******************************************************************************/

Sequence* VectorSequenceContainer::removeSequence(size_t sequenceIndex) throw (IndexOutOfBoundsException)
{
  // Copy sequence:
  if (sequenceIndex >= sequences_.size())
    throw IndexOutOfBoundsException("VectorSequenceContainer::removeSequence.", sequenceIndex, 0, sequences_.size() - 1);
  Sequence* old = sequences_[sequenceIndex];
  // Remove pointer toward old sequence:
  sequences_.erase(sequences_.begin() + static_cast<ptrdiff_t>(sequenceIndex));
  // Send copy:
  return old;
}

/******************************************************************************/

void VectorSequenceContainer::deleteSequence(size_t sequenceIndex) throw (IndexOutOfBoundsException)
{
  // Delete sequence
  if (sequenceIndex >= sequences_.size())
    throw IndexOutOfBoundsException("VectorSequenceContainer::deleteSequence.", sequenceIndex, 0, sequences_.size() - 1);
  delete sequences_[sequenceIndex];
  // Remove pointer toward old sequence:
  sequences_.erase(sequences_.begin() + static_cast<ptrdiff_t>(sequenceIndex));
}

/******************************************************************************/

void VectorSequenceContainer::addSequence(const Sequence& sequence, bool checkName) throw (Exception)
{
  // Sequence's name existence checking
  if (checkName)
  {
    // For all names in vector : throw exception if name already exists
    for (size_t i = 0; i < sequences_.size(); i++)
    {
      if (sequences_[i]->getName() == sequence.getName())
        throw Exception("VectorSequenceContainer::addSequence : Sequence '" + sequence.getName() + "' already exists in container");
    }
  }

  // New sequence's alphabet and sequence container's alphabet matching verification
  if (sequence.getAlphabet()->getAlphabetType() == getAlphabet()->getAlphabetType())
  {
    // push_back(new Sequence(sequence.getName(), sequence.getContent(), alphabet));
    sequences_.push_back(dynamic_cast<Sequence*>(sequence.clone()));
  }
  else
    throw AlphabetMismatchException("VectorSequenceContainer::addSequence : Alphabets don't match", getAlphabet(), sequence.getAlphabet());
}

void VectorSequenceContainer::addSequence(const Sequence& sequence, size_t sequenceIndex, bool checkName) throw (Exception)
{
  // Sequence's name existence checking
  if (checkName)
  {
    // For all names in vector : throw exception if name already exists
    for (size_t i = 0; i < sequences_.size(); i++)
    {
      if (sequences_[i]->getName() == sequence.getName())
        throw Exception("VectorSequenceContainer::addSequence : Sequence '" + sequence.getName() + "' already exists in container");
    }
  }

  // New sequence's alphabet and sequence container's alphabet matching verification
  if (sequence.getAlphabet()->getAlphabetType() == getAlphabet()->getAlphabetType())
  {
    // insert(begin() + pos, new Sequence(sequence.getName(), sequence.getContent(), alphabet));
    sequences_.insert(sequences_.begin() + static_cast<ptrdiff_t>(sequenceIndex), dynamic_cast<Sequence*>(sequence.clone()));
  }
  else
    throw AlphabetMismatchException("VectorSequenceContainer::addSequence : Alphabets don't match", getAlphabet(), sequence.getAlphabet());
}

/******************************************************************************/

std::vector<std::string> VectorSequenceContainer::getSequencesNames() const
{
  std::vector<std::string> names;
  for (size_t i = 0; i < sequences_.size(); i++)
  {
    names.push_back(sequences_[i]->getName());
  }
  return names;
}

/******************************************************************************/

void VectorSequenceContainer::setSequencesNames(
  const std::vector<std::string>& names,
  bool checkNames)
throw (Exception)
{
  if (names.size() != getNumberOfSequences())
    throw IndexOutOfBoundsException("VectorSequenceContainer::setSequenceNames : bad number of names", names.size(), getNumberOfSequences(), getNumberOfSequences());
  if (checkNames)
  {
    for (size_t i = 0; i < names.size(); i++)
    {
      // For all names in vector : throw exception if name already exists
      for (size_t j = 0; j < i; j++)
      {
        if (names[j] == names[i])
          throw Exception("VectorSiteContainer::setSequencesNames : Sequence's name already exists in container");
      }
    }
  }
  for (size_t i = 0; i < names.size(); i++)
  {
    sequences_[i]->setName(names[i]);
  }
}

/******************************************************************************/

void VectorSequenceContainer::clear()
{
  // Delete sequences
  for (size_t i = 0; i < sequences_.size(); i++)
  {
    delete sequences_[i];
  }
  // Delete all sequence pointers
  sequences_.clear();
}

/******************************************************************************/

void VectorSequenceContainer::setComments(size_t sequenceIndex, const Comments& comments) throw (IndexOutOfBoundsException)
{
  sequences_[sequenceIndex]->setComments(comments);
}

/******************************************************************************/

VectorSequenceContainer* VectorSequenceContainer::createEmptyContainer() const
{
  VectorSequenceContainer* vsc = new VectorSequenceContainer(getAlphabet());
  vsc->setGeneralComments(getGeneralComments());
  return vsc;
}

/******************************************************************************/

