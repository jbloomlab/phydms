//
// File: AlignedSequenceContainer.cpp
// Created by: Guillaume Deuchst
//             Julien Dutheil
// Created on: Friday August 22 2003
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

#include "AlignedSequenceContainer.h"

#include <Bpp/Text/TextTools.h>

using namespace bpp;

// From the STL:
#include <iostream>

using namespace std;

/***************************************************************************/

AlignedSequenceContainer::AlignedSequenceContainer(const OrderedSequenceContainer& osc) throw (SequenceNotAlignedException) :
  VectorSequenceContainer(osc.getAlphabet()),
  // We can't call the copy constructor because we want to use the overloaded addSequence method !!!
  positions_(),
  length_(),
  sites_()
{
  // Initializing
  for (unsigned int i = 0; i < osc.getNumberOfSequences(); i++)
  {
    addSequence(osc.getSequence(i), true);
  }

  if (osc.getNumberOfSequences() > 0)
    length_ = getSequence(0).size();  // the overloaded
  else
    length_ = 0;

  reindexSites();
  sites_.resize(length_);
  setGeneralComments(osc.getGeneralComments());
}

/***************************************************************************/

AlignedSequenceContainer& AlignedSequenceContainer::operator=(const AlignedSequenceContainer& asc)
{
  VectorSequenceContainer::operator=(asc);

  // Initializing
  length_    = asc.getNumberOfSites();
  positions_ = asc.getSitePositions();
  sites_.resize(length_);

  return *this;
}

/***************************************************************************/

AlignedSequenceContainer& AlignedSequenceContainer::operator=(const SiteContainer& sc)
{
  VectorSequenceContainer::operator=(sc);

  // Initializing
  length_    = sc.getNumberOfSites();
  positions_ = sc.getSitePositions();
  sites_.resize(length_);

  return *this;
}

/***************************************************************************/

AlignedSequenceContainer& AlignedSequenceContainer::operator=(const OrderedSequenceContainer& osc) throw (SequenceNotAlignedException)
{
  VectorSequenceContainer::operator=(osc);

  // Initializing
  length_ = 0;
  reindexSites();
  sites_.resize(length_);

  return *this;
}

/** Class destructor: *********************************************************/

AlignedSequenceContainer::~AlignedSequenceContainer()
{
  // delete all sites:

  for (unsigned int i = 0; i < sites_.size(); i++)
  {
    if (sites_[i])
      delete sites_[i];
  }
}

/***************************************************************************/

const Site& AlignedSequenceContainer::getSite(size_t i) const throw (IndexOutOfBoundsException)
{
  if (i >= length_)
    throw IndexOutOfBoundsException("AlignedSequenceContainer::getSite", i, 0, getNumberOfSites() - 1);

  // Main loop : for all sequences
  size_t n = getNumberOfSequences();
  std::vector<int> site(n);
  for (size_t j = 0; j < n; j++)
  {
    site[j] = getSequence(j)[i];
  }

  if (sites_[i])
    delete sites_[i];
  sites_[i] = new Site(site, getAlphabet(), positions_[i]);
  return *sites_[i];
}

/******************************************************************************/

void AlignedSequenceContainer::setSite(size_t pos, const Site& site, bool checkPositions) throw (Exception)
{
  // New site's alphabet and site container's alphabet matching verification
  if (pos >= getNumberOfSites())
    throw IndexOutOfBoundsException("AlignedSequenceContainer::setSite", pos, 0, getNumberOfSites() - 1);
  if (site.getAlphabet()->getAlphabetType() != getAlphabet()->getAlphabetType())
    throw AlphabetMismatchException("AlignedSequenceContainer::setSite", getAlphabet(), site.getAlphabet());

  // Check size:
  if (site.size() != getNumberOfSequences())
    throw SiteException("AlignedSequenceContainer::setSite, site does not have the appropriate length", &site);

  // Check position:
  int position = site.getPosition();
  if (checkPositions)
  {
    // For all positions in vector : throw exception if position already exists
    for (size_t i = 0; i < positions_.size(); i++)
    {
      if (positions_[i] == position)
        throw SiteException("AlignedSequenceContainer::setSite: Site position already exists in container", &site);
    }
  }

  // For all sequences
  for (size_t j = 0; j < getNumberOfSequences(); j++)
  {
    getSequence_(j).setElement(pos, site[j]);
  }
  positions_[pos] = site.getPosition();
}

/******************************************************************************/

Site* AlignedSequenceContainer::removeSite(size_t pos) throw (IndexOutOfBoundsException)
{
  if (pos >= getNumberOfSites())
    throw IndexOutOfBoundsException("AlignedSequenceContainer::removeSite", pos, 0, getNumberOfSites() - 1);

  // Get old site
  getSite(pos); // Creates the site!
  Site* old = sites_[pos];

  // For all sequences
  for (size_t j = 0; j < getNumberOfSequences(); j++)
  {
    getSequence_(j).deleteElement(pos);
  }

  // Delete site's position
  positions_.erase(positions_.begin() + static_cast<ptrdiff_t>(pos));
  length_--;

  // Actualizes the 'sites' vector:
  if (sites_[pos])
    delete sites_[pos];
  sites_.erase(sites_.begin() + static_cast<ptrdiff_t>(pos));

  // Send result
  return old;
}

/******************************************************************************/

void AlignedSequenceContainer::deleteSite(size_t pos) throw (IndexOutOfBoundsException)
{
  if (pos >= getNumberOfSites())
    throw IndexOutOfBoundsException("AlignedSequenceContainer::deleteSite", pos, 0, getNumberOfSites() - 1);

  // For all sequences
  for (size_t j = 0; j < getNumberOfSequences(); j++)
  {
    getSequence_(j).deleteElement(pos);
  }

  // Delete site's position
  positions_.erase(positions_.begin() + static_cast<ptrdiff_t>(pos));
  length_--;

  // Actualizes the 'sites' vector:
  if (sites_[pos])
    delete sites_[pos];
  sites_.erase(sites_.begin() + static_cast<ptrdiff_t>(pos));
}

/******************************************************************************/

void AlignedSequenceContainer::deleteSites(size_t siteIndex, size_t length) throw (IndexOutOfBoundsException, Exception)
{
  if (siteIndex + length > getNumberOfSites())
    throw IndexOutOfBoundsException("AlignedSequenceContainer::deleteSites", siteIndex + length, 0, getNumberOfSites() - 1);

  // For all sequences
  for (size_t j = 0; j < getNumberOfSequences(); j++)
  {
    getSequence_(j).deleteElements(siteIndex, length);
  }

  // Delete site's siteIndexition
  positions_.erase(positions_.begin() + static_cast<ptrdiff_t>(siteIndex),
      positions_.begin() + static_cast<ptrdiff_t>(siteIndex + length));
  length_ -= length;

  // Actualizes the 'sites' vector:
  for (size_t i = siteIndex; i < siteIndex + length; ++i)
  {
    if (sites_[i])
      delete sites_[i];
  }
  sites_.erase(sites_.begin() + static_cast<ptrdiff_t>(siteIndex),
      sites_.begin() + static_cast<ptrdiff_t>(siteIndex + length));
}

/******************************************************************************/

void AlignedSequenceContainer::addSite(const Site& site, bool checkPositions) throw (Exception)
{
  // New site's alphabet and site container's alphabet matching verification
  if (site.getAlphabet()->getAlphabetType() != getAlphabet()->getAlphabetType())
    throw AlphabetMismatchException("AlignedSequenceContainer::addSite");

  // Check size:
  if (site.size() != getNumberOfSequences())
    throw SiteException("AlignedSequenceContainer::addSite, site does not have the appropriate length", &site);

  // Check position:
  int position = site.getPosition();
  if (checkPositions)
  {
    // For all positions in vector : throw exception if position already exists
    for (size_t i = 0; i < positions_.size(); i++)
    {
      if (positions_[i] == position)
        throw SiteException("AlignedSequenceContainer::addSite: Site position already exists in container", &site);
    }
  }

  // For all sequences
  for (size_t j = 0; j < getNumberOfSequences(); j++)
  {
    getSequence_(j).addElement(site[j]);
  }

  length_++;
  positions_.push_back(position);

  // Actualizes the 'sites' vector:
  sites_.push_back(0);
}

/******************************************************************************/

void AlignedSequenceContainer::addSite(const Site& site, int position, bool checkPositions) throw (Exception)
{
  // New site's alphabet and site container's alphabet matching verification
  if (site.getAlphabet()->getAlphabetType() != getAlphabet()->getAlphabetType())
    throw AlphabetMismatchException("AlignedSequenceContainer::addSite");

  // Check size:
  if (site.size() != getNumberOfSequences())
    throw SiteException("AlignedSequenceContainer::addSite, site does not have the appropriate length", &site);

  // Check position:
  if (checkPositions)
  {
    // For all positions in vector : throw exception if position already exists
    for (size_t i = 0; i < positions_.size(); i++)
    {
      if (positions_[i] == position)
        throw SiteException("AlignedSequenceContainer::addSite: Site position already exists in container", &site);
    }
  }

  // For all sequences
  for (size_t j = 0; j < getNumberOfSequences(); j++)
  {
    getSequence_(j).addElement(site[j]);
  }

  length_++;
  positions_.push_back(position);

  // Actualizes the 'sites' vector:
  sites_.push_back(0);
}

/******************************************************************************/

void AlignedSequenceContainer::addSite(const Site& site, size_t siteIndex, bool checkPositions) throw (Exception)
{
  if (siteIndex >= getNumberOfSites())
    throw IndexOutOfBoundsException("AlignedSequenceContainer::addSite", siteIndex, 0, getNumberOfSites() - 1);

  // New site's alphabet and site container's alphabet matching verification
  if (site.getAlphabet()->getAlphabetType() != getAlphabet()->getAlphabetType())
    throw AlphabetMismatchException("AlignedSequenceContainer::addSite", getAlphabet(), site.getAlphabet());

  // Check size:
  if (site.size() != getNumberOfSequences())
    throw SiteException("AlignedSequenceContainer::addSite, site does not have the appropriate length", &site);

  // Check position:
  int position = site.getPosition();
  if (checkPositions)
  {
    // For all positions in vector : throw exception if position already exists
    for (size_t i = 0; i < positions_.size(); i++)
    {
      if (positions_[i] == position)
        throw SiteException("AlignedSequenceContainer::addSite: Site position already exists in container", &site);
    }
  }

  // For all sequences
  for (size_t j = 0; j < getNumberOfSequences(); j++)
  {
    getSequence_(j).addElement(siteIndex, site[j]);
  }

  length_++;
  positions_.insert(positions_.begin() + static_cast<ptrdiff_t>(siteIndex), position);

  // Actualizes the 'sites' vector:
  sites_.insert(sites_.begin() + static_cast<ptrdiff_t>(siteIndex), 0);
}

/******************************************************************************/

void AlignedSequenceContainer::addSite(const Site& site, size_t siteIndex, int position, bool checkPositions) throw (Exception)
{
  if (siteIndex >= getNumberOfSites())
    throw IndexOutOfBoundsException("AlignedSequenceContainer::addSite", siteIndex, 0, getNumberOfSites() - 1);

  // New site's alphabet and site container's alphabet matching verification
  if (site.getAlphabet()->getAlphabetType() != getAlphabet()->getAlphabetType())
    throw AlphabetMismatchException("AlignedSequenceContainer::addSite", getAlphabet(), site.getAlphabet());

  // Check size:
  if (site.size() != getNumberOfSequences())
    throw SiteException("AlignedSequenceContainer::addSite, site does not have the appropriate length", &site);

  // Check position:
  if (checkPositions)
  {
    // For all positions in vector : throw exception if position already exists
    for (size_t i = 0; i < positions_.size(); i++)
    {
      if (positions_[i] == position)
        throw SiteException("AlignedSequenceContainer::addSite: Site position already exists in container", &site);
    }
  }

  // For all sequences
  for (size_t j = 0; j < getNumberOfSequences(); j++)
  {
    getSequence_(j).addElement(siteIndex, site[j]);
  }

  length_++;
  positions_.insert(positions_.begin() + static_cast<ptrdiff_t>(siteIndex), position);

  // Actualizes the 'sites' vector:
  sites_.insert(sites_.begin() + static_cast<ptrdiff_t>(siteIndex), 0);
}

/******************************************************************************/

void AlignedSequenceContainer::reindexSites()
{
  positions_.resize(length_);
  for (size_t i = 0; i < length_; i++)
  {
    positions_[i] = static_cast<int>(i + 1); // start with 1.
  }
}

/******************************************************************************/

void AlignedSequenceContainer::setSequence(size_t i, const Sequence& sequence, bool checkName) throw (Exception)
{
  if (i >= getNumberOfSequences())
    throw IndexOutOfBoundsException("AlignedSequenceContainer::setSequence", i, 0, getNumberOfSequences() - 1);
  // if container has only one sequence
  if (getNumberOfSequences() == 1)
    length_ = sequence.size();
  if (checkSize_(sequence))
    VectorSequenceContainer::setSequence(i, sequence, checkName);
  else
    throw SequenceNotAlignedException("AlignedSequenceContainer::setSequence", &sequence);
}

/******************************************************************************/

void AlignedSequenceContainer::setSequence(const string& name, const Sequence& sequence, bool checkName) throw (Exception)
{
  // if container has only one sequence
  if (getNumberOfSequences() == 1)
    length_ = sequence.size();
  if (checkSize_(sequence))
    VectorSequenceContainer::setSequence(name, sequence, checkName);
  else
    throw SequenceNotAlignedException("AlignedSequenceContainer::setSequence", &sequence);
}

/******************************************************************************/

void AlignedSequenceContainer::addSequence(const Sequence& sequence, bool checkName) throw (Exception)
{
  // if container has only one sequence
  if (length_ == 0)
  {
    length_ = sequence.size();
    sites_.resize(length_);
    reindexSites();
  }
  if (checkSize_(sequence))
    VectorSequenceContainer::addSequence(sequence, checkName);
  else
    throw SequenceNotAlignedException("AlignedSequenceContainer::addSequence", &sequence);
}

/******************************************************************************/

void AlignedSequenceContainer::addSequence(const Sequence& sequence, size_t i, bool checkName) throw (Exception)
{
  if (i >= getNumberOfSequences())
    throw IndexOutOfBoundsException("AlignedSequenceContainer::addSequence", i, 0, getNumberOfSequences() - 1);
  // if container has only one sequence
  if (length_ == 0)
    length_ = sequence.size();
  if (checkSize_(sequence))
    VectorSequenceContainer::addSequence(sequence, i, checkName);
  else
    throw SequenceNotAlignedException("AlignedSequenceContainer::addSequence", &sequence);
}

/******************************************************************************/

void AlignedSequenceContainer::clear()
{
  length_ = 0;
  VectorSequenceContainer::clear();
}

/******************************************************************************/

AlignedSequenceContainer* AlignedSequenceContainer::createEmptyContainer() const
{
  AlignedSequenceContainer* asc = new AlignedSequenceContainer(getAlphabet());
  asc->setGeneralComments(getGeneralComments());
  return asc;
}

/******************************************************************************/


