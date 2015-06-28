//
// File MapSequenceContainer.cpp
// Authors : Guillaume Deuchst
//           Julien Dutheil
//           Sylvain Gaillard
// Last modification : Monday July 19 2004
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

#include "MapSequenceContainer.h"
#include <Bpp/Text/TextTools.h>

using namespace bpp;
using namespace std;

/******************************************************************************/

MapSequenceContainer::MapSequenceContainer(const map<string, Sequence*>& ms, const Alphabet* alpha) :
  AbstractSequenceContainer(alpha), sequences_() 
{
  for (map<string, Sequence*>::const_iterator it = ms.begin(); it != ms.end(); it++)
  {    
    addSequence(it->first, *it->second);
  }
}

/******************************************************************************/

MapSequenceContainer::MapSequenceContainer(const MapSequenceContainer& msc) :
  AbstractSequenceContainer(msc.getAlphabet()), sequences_()
{ 
  for (unsigned int i = 0; i < msc.getNumberOfSequences(); i++)
    addSequence(msc.getKey(i), msc.getSequence(i), false);
}

/******************************************************************************/

MapSequenceContainer& MapSequenceContainer::operator=(const MapSequenceContainer& msc) 
{
  clear();
  AbstractSequenceContainer::operator=(msc);
  
  // Sequences insertion
  vector<string> keys = msc.getKeys();
  for (unsigned int i = 0 ; i < getNumberOfSequences(); i++)
  {
    addSequence(keys[i], msc.getSequence(i), false);
  }

  return * this;
}

/******************************************************************************/

MapSequenceContainer::~MapSequenceContainer() 
{
  clear();
}

/******************************************************************************/

const Sequence& MapSequenceContainer::getSequence(size_t i) const throw (IndexOutOfBoundsException)
{
  // Specified sequence existence verification
  if (i < sequences_.size())
  {
    map<string, Sequence*>::const_iterator it = sequences_.begin();
    for (unsigned int j = 0; j < i; j++) it++;
    return *it->second;
  }
  throw IndexOutOfBoundsException("MapSequenceContainer::getSequence", i, 0, sequences_.size() - 1);
}

/******************************************************************************/

const Sequence& MapSequenceContainer::getSequence(const string& name) const throw (SequenceNotFoundException)
{
  // Specified sequence name research into all sequences
  for (map<string, Sequence*>::const_iterator it = sequences_.begin(); it != sequences_.end(); it++)
    if (it->second->getName() == name)
      return *it->second;
  throw SequenceNotFoundException("MapSequenceContainer::getSequence", name);
}

/******************************************************************************/

bool MapSequenceContainer::hasSequence(const string& name) const
{
  // Specified sequence name research into all sequences
  for (map<string, Sequence*>::const_iterator it = sequences_.begin(); it != sequences_.end(); it++)
    if (it->second->getName() == name)
      return true;
  return false;
}

/******************************************************************************/

Sequence& MapSequenceContainer::getSequence_(size_t i) throw (IndexOutOfBoundsException)
{
  if (i >= sequences_.size())
    throw IndexOutOfBoundsException("MapSequenceContainer::getSequence", i, 0, sequences_.size() - 1);
  map<string, Sequence*>::iterator it = sequences_.begin();
  for (size_t j = 0; j < i; j++) it++;
  return *it->second;
}

/******************************************************************************/

Sequence& MapSequenceContainer::getSequence_(const string& name) throw (SequenceNotFoundException)
{
  // Specified sequence name research into all sequences
  for (map<string, Sequence*>::iterator it = sequences_.begin(); it != sequences_.end(); it++)
    if (it->second->getName() == name)
      return *it->second;
  throw SequenceNotFoundException("MapSequenceContainer::getSequence", name);
}

/******************************************************************************/

const Sequence& MapSequenceContainer::getSequenceByKey(const string& key) const
  throw (SequenceNotFoundException)
{
  map<string, Sequence*>::const_iterator it = sequences_.find(key);
  if (it == sequences_.end())
    throw SequenceNotFoundException("MapSequenceContainer::getSequenceByKey", key);
  return *it->second;
}

/******************************************************************************/

size_t MapSequenceContainer::getSequencePosition(const string& name)
  const throw (SequenceNotFoundException)
{
  // Specified sequence name research into all sequences
  size_t pos = 0;
  for (map<string, Sequence*>::const_iterator it = sequences_.begin(); it != sequences_.end(); it++)
  {
    if (it->second->getName() == name) return pos;
    pos++;
  }

  throw SequenceNotFoundException("MapSequenceContainer::getSequencePosition", name);
}

/******************************************************************************/

void MapSequenceContainer::setSequence(size_t i, const Sequence& sequence, bool checkNames)
    throw (IndexOutOfBoundsException)
{
  // Sequence's name existence checking
  if (checkNames)
  {
    size_t j = 0;
    // For all names in map : throw exception if name already exists
    for (map<string, Sequence*>::const_iterator it = sequences_.begin(); it != sequences_.end(); it++)
    {
      if (it->second->getName() == sequence.getName()) 
        if (j != i) throw Exception("MapSequenceContainer::setSequence : Sequence's name already exists in container");
      j++;
    }
  }

  // New sequence's alphabet and sequence container's alphabet matching verification
  if (sequence.getAlphabet()->getAlphabetType() == getAlphabet()->getAlphabetType())
  {
    // Delete old sequence
    delete sequences_[getKey(i)];
    // New sequence insertion in sequence container
    sequences_[getKey(i)] = dynamic_cast<Sequence*>(sequence.clone());
  }
  else
    throw AlphabetMismatchException("MapSequenceContainer::setSequence", getAlphabet(), sequence.getAlphabet());
}

/******************************************************************************/

void MapSequenceContainer::setSequence(const string& name, const Sequence& sequence, bool checkNames) throw (SequenceNotFoundException)
{
  // Sequence's name existence checking
  if (checkNames)
  {
    // For all names in map : throw exception if name already exists
    for (map<string, Sequence*>::const_iterator it = sequences_.begin(); it != sequences_.end(); it++)
    {
      if (it->second->getName() == name) 
        if (it->second->getName() != name) 
          throw Exception("MapSequenceContainer::setSequence : Sequence's name already exists in container");
    }
  }

  // New sequence's alphabet and sequence container's alphabet matching verification
  if (sequence.getAlphabet()->getAlphabetType() == getAlphabet()->getAlphabetType())
  {
    // Delete old sequence
    delete sequences_[name];
    // New sequence insertion in sequence container
    sequences_[name] = dynamic_cast<Sequence*>(sequence.clone());
  }
  else
    throw AlphabetMismatchException("MapSequenceContainer::setSequence", getAlphabet(), sequence.getAlphabet());
}

/******************************************************************************/

void MapSequenceContainer::setSequenceByKey(const string& key, const Sequence& sequence, bool checkNames) throw (SequenceNotFoundException)
{
  // Sequence's name existence checking
  if (checkNames)
  {
    // For all names in map : throw exception if name already exists
    for (map<string, Sequence*>::const_iterator it = sequences_.begin(); it != sequences_.end(); it++)
    {
      if (it->second->getName() == sequence.getName()) 
        if (it->first != key) 
          throw Exception("MapSequenceContainer::setSequenceByKey : Sequence's name already exists in container");
    }
  }

  // New sequence's alphabet and sequence container's alphabet matching verification
  if (sequence.getAlphabet()->getAlphabetType() == getAlphabet()->getAlphabetType())
  {
    // Delete old sequence
    delete sequences_[key];
    // New sequence insertion in sequence container
    sequences_[key] = dynamic_cast<Sequence*>(sequence.clone());
  }
  else
    throw AlphabetMismatchException("MapSequenceContainer::setSequenceByKey", getAlphabet(), sequence.getAlphabet());
}

/******************************************************************************/

Sequence* MapSequenceContainer::removeSequence(size_t i) throw (IndexOutOfBoundsException)
{
  if (i >= sequences_.size())
    throw IndexOutOfBoundsException("MapSequenceContainer::removeSequence", i, 0, sequences_.size() - 1);
  map<string, Sequence*>::iterator it = sequences_.begin();
  for (size_t j = 0; j < i; j++) it++;
  Sequence* old = it->second;
  sequences_.erase(it);  
  return old;
}

/******************************************************************************/

Sequence* MapSequenceContainer::removeSequence(const string& name) throw (SequenceNotFoundException)
{
  for (map<string, Sequence*>::iterator it = sequences_.begin(); it != sequences_.end(); it++) {
    if (it->second->getName() == name)
    {
      Sequence* old = it->second;
      sequences_.erase(it);
      return old;
    }
  }
  throw SequenceNotFoundException("MapSequenceContainer::removeSequence", name);
}

/******************************************************************************/

Sequence* MapSequenceContainer::removeSequenceByKey(const string& key)throw (SequenceNotFoundException)
{
  map<string, Sequence*>::iterator it = sequences_.find(key);
  if (it == sequences_.end())
    throw SequenceNotFoundException("MapSequenceContainer::removeSequenceByKey", key);
  
  Sequence* old = it->second;
  sequences_.erase(key);
  return old;
}

/******************************************************************************/

void MapSequenceContainer::deleteSequence(size_t i) throw (IndexOutOfBoundsException)
{
  if (i >= sequences_.size())
    throw IndexOutOfBoundsException("MapSequenceContainer::deleteSequence", i, 0, sequences_.size() - 1);
  map<string, Sequence*>::iterator it = sequences_.begin();
  for (size_t j = 0; j < i; j++) it++;
  delete it->second;
  sequences_.erase(it);
}

/******************************************************************************/

void MapSequenceContainer::deleteSequence(const string& name) throw (SequenceNotFoundException)
{
  for (map<string, Sequence*>::iterator it = sequences_.begin(); it != sequences_.end(); it++) {
    if (it->second->getName() == name)
    {
      delete it->second;
      sequences_.erase(it);
      return;
    }
  }
  throw SequenceNotFoundException("MapSequenceContainer::deleteSequence", name);
}

/******************************************************************************/

void MapSequenceContainer::deleteSequenceByKey(const string& key) throw (SequenceNotFoundException)
{
  map<string, Sequence*>::iterator it = sequences_.find(key);
  if (it == sequences_.end())
    throw SequenceNotFoundException("MapSequenceContainer::deleteSequenceByKey", key);  
  delete it->second;
  sequences_.erase(key);
}

/******************************************************************************/

void MapSequenceContainer::addSequence(const string& key, const Sequence& sequence, bool checkNames) throw (Exception)
{
  // Sequence's name existence checking
  if (checkNames)
  {
    // For all names in map : throw exception if name already exists
    for (map<string, Sequence*>::const_iterator it = sequences_.begin(); it != sequences_.end(); it++)
    {
      if (it->second->getName() == sequence.getName()) 
        throw Exception("MapSequenceContainer::addSequence: Sequence '" + sequence.getName() + ", already exists in container");
    }
  }
  
  // Check if the key is not used
  for (map<string, Sequence*>::const_iterator it = sequences_.begin(); it != sequences_.end(); it++)
    if (key == it->first)
      throw Exception("MapSequenceContainer::addSequence: key already in use. (" + key + ")");
  
  // New sequence's alphabet and sequence container's alphabet matching verification
  if (sequence.getAlphabet()->getAlphabetType() == getAlphabet()->getAlphabetType())
    sequences_.insert(make_pair(key, dynamic_cast<Sequence*>(sequence.clone())));
  else throw AlphabetMismatchException("MapSequenceContainer::addSequence", getAlphabet(), sequence.getAlphabet());
}

/******************************************************************************/

vector<string> MapSequenceContainer::getKeys() const
{
  vector<string> keys;
  for (map<string, Sequence*>::const_iterator it = sequences_.begin(); it != sequences_.end(); it++)
    keys.push_back(it->first);
  return keys;
}

/******************************************************************************/

string MapSequenceContainer::getKey(size_t pos) const throw (IndexOutOfBoundsException)
{
  if (pos >= getNumberOfSequences())
    throw IndexOutOfBoundsException("MapSequenceContainer::getKey", pos, 0, sequences_.size() - 1);
  map<string, Sequence*>::const_iterator it = sequences_.begin();
  for (size_t i = 0; i < pos; i++) it++;
  return it->first;
}

/******************************************************************************/

string MapSequenceContainer::getKey(const string& name) const throw (SequenceNotFoundException)
{
  try
  {
    return getKey(getSequencePosition(name));
  }
  catch (SequenceNotFoundException & snfe)
  {
    throw SequenceNotFoundException("MapSequenceContainer::getKey", snfe.getSequenceId());
  }
}
  
/******************************************************************************/

void MapSequenceContainer::setComments(size_t pos, const Comments& comments) throw (IndexOutOfBoundsException)
{
  if (pos >= getNumberOfSequences())
    throw IndexOutOfBoundsException("MapSequenceContainer::setComments", pos, 0, sequences_.size() - 1);
  map<string, Sequence*>::iterator it = sequences_.begin();
  for (size_t i = 0 ; i < pos ; i++) it++;
  it->second->setComments(comments);
}

/******************************************************************************/

vector<string> MapSequenceContainer::getSequencesNames() const
{
  vector<string> names;
  for (map<string, Sequence*>::const_iterator it = sequences_.begin(); it != sequences_.end(); it++)
    names.push_back(it->second->getName());
  return names;
}

/******************************************************************************/

void MapSequenceContainer::setSequencesNames(const vector<string>& names, bool checkNames) throw (Exception)
{
  if (names.size() != getNumberOfSequences())
    throw IndexOutOfBoundsException("MapSequenceContainer::setSequenceNames : bad number of names", names.size(), getNumberOfSequences(), getNumberOfSequences());
  if (checkNames) {
    // check if there is no repeat names in teh vector
    for (size_t i = 0 ; i < names.size() ; i++)
      for (unsigned int j = 0 ; j < i ; j++)
        if (names[j] == names[i])
          throw Exception("MapSequenceContainer::setSequencesNames: Sequence's name already exists in container");
  }
  map<string, Sequence*>::iterator it = sequences_.begin();
  for (size_t i = 0 ; i < names.size() ; i++)
  {
    it->second->setName(names[i]);
    it++;
  }
}

/******************************************************************************/

void MapSequenceContainer::clear()
{
  // Delete sequences
  for (map<string, Sequence *>::iterator it = sequences_.begin(); it != sequences_.end(); it++)
    delete it->second;
  // Delete all sequence pointers
  sequences_.clear();
}

/******************************************************************************/

MapSequenceContainer* MapSequenceContainer::createEmptyContainer() const
{ 
  MapSequenceContainer* msc = new MapSequenceContainer(getAlphabet());
  msc->setGeneralComments(getGeneralComments());
  return msc;
}

/******************************************************************************/

