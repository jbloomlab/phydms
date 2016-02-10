//
// File: VectorSiteContainer.cpp
// Created by: Julien Dutheil
// Created on: Mon Oct  6 11:50:40 2003
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

#include "VectorSiteContainer.h"

#include <iostream>

using namespace std;

#include <Bpp/Text/TextTools.h>

using namespace bpp;

/** Class constructors: *******************************************************/

VectorSiteContainer::VectorSiteContainer(
  const std::vector<const Site*>& vs,
  const Alphabet* alpha,
  bool checkPositions)
throw (Exception) :
  AbstractSequenceContainer(alpha),
  sites_(0),
  names_(0),
  comments_(0),
  sequences_(0)
{
  if (vs.size() == 0)
    throw Exception("VectorSiteContainer::VectorSiteContainer. Empty site set.");
  // Seq names and comments:
  size_t nbSeq = vs[0]->size();
  names_.resize(nbSeq);
  comments_.resize(nbSeq);
  for (size_t i = 0; i < nbSeq; i++)
  {
    names_[i]    = "Seq_" + TextTools::toString(i);
    comments_[i] = new Comments();
  }
  // Now try to add each site:
  for (size_t i = 0; i < vs.size(); i++)
  {
    addSite(*vs[i], checkPositions); // This may throw an exception if position argument already exists or is size is not valid.
  }

  sequences_.resize(nbSeq);
}

/******************************************************************************/

VectorSiteContainer::VectorSiteContainer(size_t size, const Alphabet* alpha) :
  AbstractSequenceContainer(alpha),
  sites_(0),
  names_(size),
  comments_(size),
  sequences_(size)
{
  // Seq names and comments:
  for (size_t i = 0; i < size; i++)
  {
    names_[i]    = string("Seq_") + TextTools::toString(i);
    comments_[i] = new Comments();
  }
}

/******************************************************************************/

VectorSiteContainer::VectorSiteContainer(const std::vector<std::string>& names, const Alphabet* alpha) :
  AbstractSequenceContainer(alpha),
  sites_(0),
  names_(names.size()),
  comments_(names.size()),
  sequences_(names.size())
{
  // Seq names and comments:
  for (size_t i = 0; i < names.size(); i++)
  {
    names_[i]    = names[i];
    comments_[i] = new Comments();
  }
}

/******************************************************************************/

VectorSiteContainer::VectorSiteContainer(const Alphabet* alpha) :
  AbstractSequenceContainer(alpha),
  sites_(0),
  names_(0),
  comments_(0),
  sequences_(0)
{}

/******************************************************************************/

VectorSiteContainer::VectorSiteContainer(const VectorSiteContainer& vsc) :
  AbstractSequenceContainer(vsc),
  sites_(0),
  names_(vsc.names_),
  comments_(vsc.getNumberOfSequences()),
  sequences_(vsc.getNumberOfSequences())
{
  // Now try to add each site:
  for (size_t i = 0; i < vsc.getNumberOfSites(); i++)
  {
    addSite(vsc.getSite(i), false); // We assume that positions are correct.
  }
  // Seq comments:
  for (size_t i = 0; i < vsc.getNumberOfSequences(); i++)
  {
    comments_[i] = new Comments(vsc.getComments(i));
  }
}

/******************************************************************************/

VectorSiteContainer::VectorSiteContainer(const SiteContainer& sc) :
  AbstractSequenceContainer(sc),
  sites_(0),
  names_(sc.getSequencesNames()),
  comments_(sc.getNumberOfSequences()),
  sequences_(sc.getNumberOfSequences())
{
  // Now try to add each site:
  for (size_t i = 0; i < sc.getNumberOfSites(); i++)
  {
    addSite(sc.getSite(i), false); // We assume that positions are correct.
  }
  // Seq comments:
  for (size_t i = 0; i < sc.getNumberOfSequences(); i++)
  {
    comments_[i] = new Comments(sc.getComments(i));
  }
}

/******************************************************************************/

VectorSiteContainer::VectorSiteContainer(const OrderedSequenceContainer& osc) :
  AbstractSequenceContainer(osc),
  sites_(0),
  names_(0),
  comments_(0),
  sequences_(0)
{
  for (size_t i = 0; i < osc.getNumberOfSequences(); i++)
  {
    addSequence(osc.getSequence(i), false);
  }
  reindexSites();
}

/******************************************************************************/

VectorSiteContainer::VectorSiteContainer(const SequenceContainer& sc) :
  AbstractSequenceContainer(sc),
  sites_(0),
  names_(0),
  comments_(0),
  sequences_(0)
{
  vector<string> names = sc.getSequencesNames();
  for (size_t i = 0; i < names.size(); i++)
  {
    addSequence(sc.getSequence(names[i]), false);
  }
  reindexSites();
}

/******************************************************************************/

VectorSiteContainer& VectorSiteContainer::operator=(const VectorSiteContainer& vsc)
{
  clear();
  AbstractSequenceContainer::operator=(vsc);
  // Seq names:
  names_.resize(vsc.getNumberOfSequences());
  setSequencesNames(vsc.getSequencesNames(), true);
  // Now try to add each site:
  for (size_t i = 0; i < vsc.getNumberOfSites(); i++)
  {
    addSite(vsc.getSite(i), false); // We assume that positions are correct.
  }
  // Seq comments:
  size_t nbSeq = vsc.getNumberOfSequences();
  comments_.resize(nbSeq);
  for (size_t i = 0; i < nbSeq; i++)
  {
    comments_[i] = new Comments(vsc.getComments(i));
  }
  sequences_.resize(nbSeq);

  return *this;
}

/******************************************************************************/

VectorSiteContainer& VectorSiteContainer::operator=(const SiteContainer& sc)
{
  clear();
  AbstractSequenceContainer::operator=(sc);
  // Seq names:
  names_.resize(sc.getNumberOfSequences());
  setSequencesNames(sc.getSequencesNames(), true);
  // Now try to add each site:
  for (size_t i = 0; i < sc.getNumberOfSites(); i++)
  {
    addSite(sc.getSite(i), false); // We assume that positions are correct.
  }
  // Seq comments:
  size_t nbSeq = sc.getNumberOfSequences();
  comments_.resize(nbSeq);
  for (size_t i = 0; i < nbSeq; i++)
  {
    comments_[i] = new Comments(sc.getComments(i));
  }
  sequences_.resize(nbSeq);

  return *this;
}

/******************************************************************************/

VectorSiteContainer& VectorSiteContainer::operator=(const OrderedSequenceContainer& osc)
{
  clear();
  AbstractSequenceContainer::operator=(osc);

  size_t nbSeq = osc.getNumberOfSequences();
  for (size_t i = 0; i < nbSeq; i++)
  {
    addSequence(osc.getSequence(i), false);
  }
  reindexSites();

  return *this;
}

/******************************************************************************/

VectorSiteContainer& VectorSiteContainer::operator=(const SequenceContainer& sc)
{
  clear();
  AbstractSequenceContainer::operator=(sc);

  vector<string> names = sc.getSequencesNames();
  for (size_t i = 0; i < names.size(); i++)
  {
    addSequence(sc.getSequence(names[i]), false);
  }
  reindexSites();

  return *this;
}

/******************************************************************************/

const Site& VectorSiteContainer::getSite(size_t i) const throw (IndexOutOfBoundsException)
{
  if (i >= getNumberOfSites())
    throw IndexOutOfBoundsException("VectorSiteContainer::getSite.", i, 0, getNumberOfSites() - 1);
  return *sites_[i];
}

/******************************************************************************/

void VectorSiteContainer::setSite(size_t pos, const Site& site, bool checkPositions) throw (Exception)
{
  if (pos >= getNumberOfSites())
    throw IndexOutOfBoundsException("VectorSiteContainer::setSite.", pos, 0, getNumberOfSites() - 1);

  // Check size:
  if (site.size() != getNumberOfSequences())
    throw SiteException("AlignedSequenceContainer::addSite. Site does not have the appropriate length", &site);

  // New site's alphabet and site container's alphabet matching verification
  if (site.getAlphabet()->getAlphabetType() != getAlphabet()->getAlphabetType())
    throw AlphabetMismatchException("VectorSiteContainer::setSite", getAlphabet(), site.getAlphabet());

  // Check position:
  if (checkPositions)
  {
    int position = site.getPosition();
    // For all positions in vector : throw exception if position already exists
    for (size_t i = 0; i < sites_.size(); i++)
    {
      if (sites_[i]->getPosition() == position)
        throw SiteException("VectorSiteContainer::setSite: Site position already exists in container", &site);
    }
  }
  delete sites_[pos];
  sites_[pos] = dynamic_cast<Site*>(site.clone());
}

/******************************************************************************/

Site* VectorSiteContainer::removeSite(size_t i) throw (IndexOutOfBoundsException)
{
  if (i >= getNumberOfSites())
    throw IndexOutOfBoundsException("VectorSiteContainer::removeSite.", i, 0, getNumberOfSites() - 1);
  Site* site = sites_[i];
  sites_.erase(sites_.begin() + static_cast<ptrdiff_t>(i));
  return site;
}

/******************************************************************************/

void VectorSiteContainer::deleteSite(size_t i) throw (IndexOutOfBoundsException)
{
  if (i >= getNumberOfSites())
    throw IndexOutOfBoundsException("VectorSiteContainer::deleteSite.", i, 0, getNumberOfSites() - 1);
  delete sites_[i];
  sites_.erase(sites_.begin() + static_cast<ptrdiff_t>(i));
}

/******************************************************************************/

void VectorSiteContainer::deleteSites(size_t siteIndex, size_t length) throw (IndexOutOfBoundsException)
{
  if (siteIndex + length > getNumberOfSites())
    throw IndexOutOfBoundsException("VectorSiteContainer::deleteSites.", siteIndex + length, 0, getNumberOfSites() - 1);
  for (size_t i = siteIndex; i < siteIndex + length; ++i)
  {
    delete sites_[i];
  }
  sites_.erase(sites_.begin() + static_cast<ptrdiff_t>(siteIndex), sites_.begin() + static_cast<ptrdiff_t>(siteIndex + length));
}

/******************************************************************************/

void VectorSiteContainer::addSite(const Site& site, bool checkPositions) throw (Exception)
{
  // Check size:
  if (site.size() != getNumberOfSequences())
    throw SiteException("VectorSiteContainer::addSite. Site does not have the appropriate length", &site);

  // New site's alphabet and site container's alphabet matching verification
  if (site.getAlphabet()->getAlphabetType() != getAlphabet()->getAlphabetType())
  {
    throw AlphabetMismatchException("VectorSiteContainer::addSite", getAlphabet(), site.getAlphabet());
  }

  // Check position:
  if (checkPositions)
  {
    int position = site.getPosition();
    // For all positions in vector : throw exception if position already exists
    for (size_t i = 0; i < sites_.size(); i++)
    {
      if (sites_[i]->getPosition() == position)
        throw SiteException("VectorSiteContainer::addSite. Site position already exists in container", &site);
    }
  }

  sites_.push_back(dynamic_cast<Site*>(site.clone()));
}

/******************************************************************************/

void VectorSiteContainer::addSite(const Site& site, int position, bool checkPositions) throw (Exception)
{
  // Check size:
  if (site.size() != getNumberOfSequences())
    throw SiteException("VectorSiteContainer::addSite. Site does not have the appropriate length", &site);

  // New site's alphabet and site container's alphabet matching verification
  if (site.getAlphabet()->getAlphabetType() != getAlphabet()->getAlphabetType())
  {
    throw AlphabetMismatchException("VectorSiteContainer::addSite", getAlphabet(), site.getAlphabet());
  }

  // Check position:
  if (checkPositions)
  {
    // For all positions in vector : throw exception if position already exists
    for (size_t i = 0; i < sites_.size(); i++)
    {
      if (sites_[i]->getPosition() == position)
        throw SiteException("VectorSiteContainer::addSite. Site position already exists in container", &site);
    }
  }
  Site* copy = dynamic_cast<Site*>(site.clone());
  copy->setPosition(position);
  sites_.push_back(copy);
}

/******************************************************************************/

void VectorSiteContainer::addSite(const Site& site, size_t siteIndex, bool checkPositions) throw (Exception)
{
  if (siteIndex >= getNumberOfSites())
    throw IndexOutOfBoundsException("VectorSiteContainer::addSite", siteIndex, 0, getNumberOfSites() - 1);

  // Check size:
  if (site.size() != getNumberOfSequences())
    throw SiteException("VectorSiteContainer::addSite. Site does not have the appropriate length", &site);

  // New site's alphabet and site container's alphabet matching verification
  if (site.getAlphabet()->getAlphabetType() != getAlphabet()->getAlphabetType())
  {
    throw AlphabetMismatchException("VectorSiteContainer::addSite", getAlphabet(), site.getAlphabet());
  }

  // Check position:
  if (checkPositions)
  {
    int position = site.getPosition();
    // For all positions in vector : throw exception if position already exists
    for (size_t i = 0; i < sites_.size(); i++)
    {
      if (sites_[i]->getPosition() == position)
        throw SiteException("VectorSiteContainer::addSite. Site position already exists in container", &site);
    }
  }

  // insert(begin() + pos, new Site(site));
  sites_.insert(sites_.begin() + static_cast<ptrdiff_t>(siteIndex), dynamic_cast<Site*>(site.clone()));
}

/******************************************************************************/

void VectorSiteContainer::addSite(const Site& site, size_t siteIndex, int position, bool checkPositions) throw (Exception)
{
  if (siteIndex >= getNumberOfSites())
    throw IndexOutOfBoundsException("VectorSiteContainer::addSite", siteIndex, 0, getNumberOfSites() - 1);

  // Check size:
  if (site.size() != getNumberOfSequences())
    throw SiteException("VectorSiteContainer::addSite. Site does not have the appropriate length", &site);

  // New site's alphabet and site container's alphabet matching verification
  if (site.getAlphabet()->getAlphabetType() != getAlphabet()->getAlphabetType())
  {
    throw AlphabetMismatchException("VectorSiteContainer::addSite", getAlphabet(), site.getAlphabet());
  }

  // Check position:
  if (checkPositions)
  {
    // For all positions in vector : throw exception if position already exists
    for (size_t i = 0; i < sites_.size(); i++)
    {
      if (sites_[i]->getPosition() == position)
        throw SiteException("VectorSiteContainer::addSite. Site position already exists in container", &site);
    }
  }

  Site* copy = dynamic_cast<Site*>(site.clone());
  copy->setPosition(position);
  sites_.insert(sites_.begin() + static_cast<ptrdiff_t>(siteIndex), copy);
}

/******************************************************************************/

size_t VectorSiteContainer::getNumberOfSites() const
{
  return sites_.size();
}

/******************************************************************************/

void VectorSiteContainer::reindexSites()
{
  int pos = 1; // first position is 1.
  for (vector<Site*>::iterator i = sites_.begin(); i < sites_.end(); i++)
  {
    (*i)->setPosition(pos++);
  }
}

/******************************************************************************/

Vint VectorSiteContainer::getSitePositions() const
{
  Vint positions(sites_.size());
  for (size_t i = 0; i < sites_.size(); i++)
  {
    positions[i] = sites_[i]->getPosition();
  }
  return positions;
}

void VectorSiteContainer::setSitePositions(Vint vPositions)
{
  if (vPositions.size() != sites_.size())
    throw BadSizeException("VectorSiteContainer::setSitePositions bad size of positions vector", vPositions.size(), sites_.size());
  
  size_t pos = 0; // first position is 1.
  for (vector<Site*>::iterator i = sites_.begin(); i < sites_.end(); i++)
  {
    (*i)->setPosition(vPositions[pos++]);
  }
}

/******************************************************************************/

const Sequence& VectorSiteContainer::getSequence(size_t i) const throw (IndexOutOfBoundsException)
{
  if (i >= getNumberOfSequences())
    throw IndexOutOfBoundsException("VectorSiteContainer::getSequence.", i, 0, getNumberOfSequences() - 1);

  // Main loop : for all sites
  size_t n = getNumberOfSites();
  vector<int> sequence(n);
  for (size_t j = 0; j < n; j++)
  {
    sequence[j] = (*sites_[j])[i];
  }
  if (sequences_[i])
    delete sequences_[i];
  sequences_[i] = new BasicSequence(names_[i], sequence, *comments_[i], getAlphabet());
  return *sequences_[i];
}

/******************************************************************************/

const Sequence& VectorSiteContainer::getSequence(const string& name) const throw (SequenceNotFoundException)
{
  // Look for sequence name:
  size_t pos = getSequencePosition(name);
  return getSequence(pos);
}

/******************************************************************************/

bool VectorSiteContainer::hasSequence(const string& name) const
{
  // Look for sequence name:
  for (size_t pos = 0; pos < names_.size(); pos++)
  {
    if (names_[pos] == name)
      return true;
  }
  return false;
}

/******************************************************************************/

size_t VectorSiteContainer::getSequencePosition(const string& name) const throw (SequenceNotFoundException)
{
  // Look for sequence name:
  for (size_t pos = 0; pos < names_.size(); pos++)
  {
    if (names_[pos] == name)
      return pos;
  }
  throw SequenceNotFoundException("VectorSiteContainer::getSequencePosition().", name);
}

/******************************************************************************/

void VectorSiteContainer::setSequence(const string& name, const Sequence& sequence, bool checkNames) throw (Exception)
{
  // Look for sequence name:
  size_t pos = getSequencePosition(name);
  setSequence(pos, sequence, checkNames);
}

/******************************************************************************/

void VectorSiteContainer::setSequence(size_t pos, const Sequence& sequence, bool checkNames)
throw (Exception)
{
  if (pos >= getNumberOfSequences())
    throw IndexOutOfBoundsException("VectorSiteContainer::setSequence", pos, 0, getNumberOfSequences() - 1);

  // New sequence's alphabet and site container's alphabet matching verification
  if (sequence.getAlphabet()->getAlphabetType() != getAlphabet()->getAlphabetType())
    throw AlphabetMismatchException("VectorSiteContainer::addSite", getAlphabet(), sequence.getAlphabet());

  // If the container has only one sequence, we set the size to the size of this sequence:
  if (getNumberOfSequences() == 1)
    realloc(sequence.size());

  if (sequence.size() != sites_.size())
    throw SequenceException("VectorSiteContainer::setSequence. Sequence has not the appropriate length.", &sequence);

  if (checkNames)
  {
    for (size_t i = 0; i < names_.size(); i++)
    {
      if (i != pos && sequence.getName() == names_[i])
        throw SequenceException("VectorSiteContainer::settSequence. Name already exists in container.", &sequence);
    }
  }
  // Update name:
  names_[pos] = sequence.getName();
  // Update elements at each site:
  for (size_t i = 0; i < sites_.size(); i++)
  {
    sites_[i]->setElement(pos, sequence.getValue(i));
  }
  // Update comments:
  if (comments_[pos])
    delete comments_[pos];
  comments_[pos] = new Comments(sequence.getComments());
  // Update sequences:
  if (sequences_[pos])
    delete sequences_[pos];
  sequences_[pos] = 0;
}

/******************************************************************************/

Sequence* VectorSiteContainer::removeSequence(size_t i) throw (IndexOutOfBoundsException)
{
  if (i >= getNumberOfSequences())
    throw IndexOutOfBoundsException("VectorSiteContainer::removeSequence.", i, 0, getNumberOfSequences() - 1);

  getSequence(i); // Actuallizes pointer.
  Sequence* sequence = sequences_[i];
  for (size_t j = 0; j < sites_.size(); j++)
  {
    // For each site:
    sites_[j]->deleteElement(i);
  }

  // Now actualize names and comments:
  names_.erase(names_.begin() + static_cast<ptrdiff_t>(i));
  if (comments_[i])
    delete comments_[i];
  comments_.erase(comments_.begin() + static_cast<ptrdiff_t>(i));
  // We remove the sequence, so the destruction of the sequence is up to the user:
  // if (sequences_[i] != 0) delete sequences_[i];
  sequences_.erase(sequences_.begin() + static_cast<ptrdiff_t>(i));
  return sequence;
}

/******************************************************************************/

Sequence* VectorSiteContainer::removeSequence(const string& name) throw (SequenceNotFoundException)
{
  // Look for sequence name:
  size_t pos = getSequencePosition(name);
  return removeSequence(pos);
}

/******************************************************************************/

void VectorSiteContainer::deleteSequence(size_t i) throw (IndexOutOfBoundsException)
{
  if (i >= getNumberOfSequences())
    throw IndexOutOfBoundsException("VectorSiteContainer::demeteSequence.", i, 0, getNumberOfSequences() - 1);
  for (size_t j = 0; j < sites_.size(); j++)
  {
    sites_[j]->deleteElement(i);
  }

  // Now actualize names and comments:
  names_.erase(names_.begin() + static_cast<ptrdiff_t>(i));
  if (comments_[i])
    delete comments_[i];
  comments_.erase(comments_.begin() + static_cast<ptrdiff_t>(i));
  if (sequences_[i])
    delete sequences_[i];
  sequences_.erase(sequences_.begin() + static_cast<ptrdiff_t>(i));
}

/******************************************************************************/

void VectorSiteContainer::deleteSequence(const string& name) throw (SequenceNotFoundException)
{
  // Look for sequence name:
  size_t pos = getSequencePosition(name);
  deleteSequence(pos);
}

/******************************************************************************/

void VectorSiteContainer::addSequence(const Sequence& sequence, bool checkNames) throw (Exception)
{
  // If the container has no sequence, we set the size to the size of this sequence:
  if (getNumberOfSequences() == 0)
    realloc(sequence.size());

  // New sequence's alphabet and site container's alphabet matching verification
  if (sequence.getAlphabet()->getAlphabetType() != getAlphabet()->getAlphabetType())
    throw AlphabetMismatchException("VectorSiteContainer::addSequence", getAlphabet(), sequence.getAlphabet());

  if (sequence.size() != sites_.size())
    throw SequenceException("VectorSiteContainer::addSequence. Sequence has not the appropriate length: " + TextTools::toString(sequence.size()) + ", should be " + TextTools::toString(sites_.size()) + ".", &sequence);

  if (checkNames)
  {
    for (size_t i = 0; i < names_.size(); i++)
    {
      if (sequence.getName() == names_[i])
        throw SequenceException("VectorSiteContainer::addSequence. Name already exists in container.", &sequence);
    }
  }

  // Append name:
  names_.push_back(sequence.getName());

  // Append elements at each site:
  for (size_t i = 0; i < sites_.size(); i++)
  {
    sites_[i]->addElement(sequence.getValue(i));
  }

  // Append comments:
  comments_.push_back(new Comments(sequence.getComments()));

  // Sequences pointers:
  sequences_.push_back(0);
}

/******************************************************************************/

void VectorSiteContainer::addSequence(
  const Sequence& sequence,
  size_t pos,
  bool checkNames)
throw (Exception)
{
  if (pos >= getNumberOfSequences())
    throw IndexOutOfBoundsException("VectorSiteContainer::addSequence.", pos, 0, getNumberOfSequences() - 1);
  if (sequence.size() != sites_.size())
    throw SequenceNotAlignedException("VectorSiteContainer::setSequence", &sequence);

  // New sequence's alphabet and site container's alphabet matching verification
  if (sequence.getAlphabet()->getAlphabetType() != getAlphabet()->getAlphabetType())
  {
    throw AlphabetMismatchException("VectorSiteContainer::addSite", getAlphabet(), sequence.getAlphabet());
  }

  if (checkNames)
  {
    for (size_t i = 0; i < names_.size(); i++)
    {
      if (sequence.getName() == names_[i])
        throw SequenceException("VectorSiteContainer::addSequence. Name already exists in container.", &sequence);
    }
  }

  for (size_t i = 0; i < sites_.size(); i++)
  {
    // For each site:
    sites_[i]->addElement(pos, sequence.getValue(i));
  }
  // Actualize names and comments:
  names_.insert(names_.begin() + static_cast<ptrdiff_t>(pos), sequence.getName());
  comments_.insert(comments_.begin() + static_cast<ptrdiff_t>(pos), new Comments(sequence.getComments()));
  sequences_.insert(sequences_.begin() + static_cast<ptrdiff_t>(pos), 0);
}

/******************************************************************************/

void VectorSiteContainer::clear()
{
  // Must delete all sites in the container:
  for (size_t i = 0; i < sites_.size(); i++)
  {
    delete sites_[i];
  }

  // must delete all comments too:
  for (size_t i = 0; i < comments_.size(); i++)
  {
    if (comments_[i] != 0)
      delete comments_[i];
  }

  // Delete all sequences retrieved:
  for (size_t i = 0; i < sequences_.size(); i++)
  {
    if (sequences_[i] != 0)
      delete (sequences_[i]);
  }

  // Delete all sites pointers
  sites_.clear();
  names_.clear();
  comments_.clear();
  sequences_.clear();
}

/******************************************************************************/

void VectorSiteContainer::realloc(size_t n)
{
  clear();
  sites_.resize(n);
  for (size_t i = 0; i < n; i++)
  {
    sites_[i] = new Site(getAlphabet());
  }
  reindexSites();
}

/******************************************************************************/

vector<string> VectorSiteContainer::getSequencesNames() const
{
  return names_;
}

/******************************************************************************/

void VectorSiteContainer::setSequencesNames(
  const vector<string>& names,
  bool checkNames)
throw (Exception)
{
  if (names.size() != getNumberOfSequences())
    throw IndexOutOfBoundsException("VectorSiteContainer::setSequenceNames: bad number of names.", names.size(), getNumberOfSequences(), getNumberOfSequences());
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
    names_[i] = names[i];
  }
}

/******************************************************************************/

void VectorSiteContainer::setComments(size_t sequenceIndex, const Comments& comments) throw (IndexOutOfBoundsException)
{
  comments_[sequenceIndex] = new Comments(comments);
}

/******************************************************************************/

VectorSiteContainer* VectorSiteContainer::createEmptyContainer() const
{
  VectorSiteContainer* vsc = new VectorSiteContainer(getAlphabet());
  vsc->setGeneralComments(getGeneralComments());
  return vsc;
}

/******************************************************************************/

