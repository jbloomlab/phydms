//
// File: CompressedCompressedVectorSiteContainer.cpp
// Created by: Julien Dutheil
// Created on: Wed Dec  16 12:08 2009
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

#include "CompressedVectorSiteContainer.h"
#include <Bpp/Text/TextTools.h>

#include <iostream>

using namespace std;

using namespace bpp;

/** Class constructors: *******************************************************/

CompressedVectorSiteContainer::CompressedVectorSiteContainer(
  const std::vector<const Site*>& vs,
  const Alphabet* alpha)
throw (Exception) :
  AbstractSequenceContainer(alpha),
  sites_(0),
  index_(0),
  names_(0),
  comments_(0),
  sequences_(0)
{
  if (vs.size() == 0) throw Exception("CompressedVectorSiteContainer::CompressedVectorSiteContainer. Empty site set.");
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
    addSite(*vs[i]); // This may throw an exception if position argument already exists or is size is not valid.
  }

  sequences_.resize(nbSeq);
}

/******************************************************************************/

CompressedVectorSiteContainer::CompressedVectorSiteContainer(size_t size, const Alphabet* alpha) :
  AbstractSequenceContainer(alpha),
  sites_(0),
  index_(0),
  names_(size),
  comments_(size),
  sequences_(size)
{
  // Seq names and comments:
  for (size_t i = 0; i < size; i++)
  {
    names_[i]    = "Seq_" + TextTools::toString(i);
    comments_[i] = new Comments();
  }
}

/******************************************************************************/

CompressedVectorSiteContainer::CompressedVectorSiteContainer(const std::vector<std::string>& names, const Alphabet* alpha) :
  AbstractSequenceContainer(alpha),
  sites_(0),
  index_(0),
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

CompressedVectorSiteContainer::CompressedVectorSiteContainer(const Alphabet* alpha) :
  AbstractSequenceContainer(alpha),
  sites_(0),
  index_(0),
  names_(0),
  comments_(0),
  sequences_(0)
{}

/******************************************************************************/

CompressedVectorSiteContainer::CompressedVectorSiteContainer(const CompressedVectorSiteContainer& vsc) :
  AbstractSequenceContainer(vsc),
  sites_(vsc.sites_.size()),
  index_(vsc.index_),
  names_(vsc.names_),
  comments_(vsc.getNumberOfSequences()),
  sequences_(vsc.getNumberOfSequences())
{
  // Now try to add each site:
  sites_.resize(vsc.sites_.size());
  for (size_t i = 0; i < vsc.sites_.size(); i++)
  {
    sites_[i] = dynamic_cast<Site*>(vsc.sites_[i]->clone());
  }
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

CompressedVectorSiteContainer::CompressedVectorSiteContainer(const SiteContainer& sc) :
  AbstractSequenceContainer(sc.getAlphabet()),
  sites_(0),
  index_(0),
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

CompressedVectorSiteContainer& CompressedVectorSiteContainer::operator=(const CompressedVectorSiteContainer& vsc)
{
  AbstractSequenceContainer::operator=(vsc);
  // Seq names:
  names_ = vsc.names_;
  // Now try to add each site:
  sites_.resize(vsc.sites_.size());
  for (size_t i = 0; i < vsc.sites_.size(); i++)
  {
    sites_[i] = dynamic_cast<Site*>(vsc.sites_[i]->clone());
  }
  index_ = vsc.index_;
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

CompressedVectorSiteContainer& CompressedVectorSiteContainer::operator=(const SiteContainer& sc)
{
  clear();
  AbstractSequenceContainer::operator=(sc);
  // Seq names:
  names_ = sc.getSequencesNames();
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

const Site& CompressedVectorSiteContainer::getSite(size_t i) const throw (IndexOutOfBoundsException)
{
  if (i >= getNumberOfSites())
    throw IndexOutOfBoundsException("CompressedVectorSiteContainer::getSite.", i, 0, getNumberOfSites() - 1);
  return *sites_[index_[i]];
}

/******************************************************************************/

void CompressedVectorSiteContainer::setSite(size_t pos, const Site& site, bool checkPositions) throw (Exception)
{
  if (pos >= getNumberOfSites()) throw IndexOutOfBoundsException("CompressedVectorSiteContainer::setSite.", pos, 0, getNumberOfSites() - 1);

  // Check size:
  if (site.size() != getNumberOfSequences()) throw SiteException("AlignedSequenceContainer::addSite. Site does not have the appropriate length", &site);

  // New site's alphabet and site container's alphabet matching verification
  if (site.getAlphabet()->getAlphabetType() != getAlphabet()->getAlphabetType())
    throw AlphabetMismatchException("CompressedVectorSiteContainer::setSite", getAlphabet(), site.getAlphabet());
  
  size_t current = index_[pos];
  size_t siteIndex = getSiteIndex_(site);
  if (siteIndex == current)
  {
    //Nothing to do here, this is the same site.
  }
  else if (siteIndex < sites_.size())
  {
    //The new site is already in the list, si we just update the index:
    index_[pos] = siteIndex;

    //We have to check if the previous pattern was unique, and if so, remove it and update indices:
    bool test = true;
    for (size_t i = 0; test && i < index_.size(); ++i)
    {
      if (index_[i] == current)
      {
        //There is another site, so nothing to do...
        test = false;
      }
    }
    if (test)
    {
      //There was no other site pointing toward this pattern, so we remove it.
      delete sites_[current];
      sites_.erase(sites_.begin() + static_cast<ptrdiff_t>(current));
      //Now we have to correct all indices:
      for (size_t i = 0; i < index_.size(); ++i)
      {
        if (index_[i] > current) index_[i]--;
      }
    }
  }
  else
  {
    //This is a new pattern, and we have to add it to the list...
    Site* copy = dynamic_cast<Site*>(site.clone());

    //Now we have to check if the previous pattern was unique, and if so,
    //replace it with the new one. Otherwise, add the new site at the end of the list.
    bool test = true;
    for (size_t i = 0; test && i < index_.size(); ++i)
    {
      if (i != pos && index_[i] == current)
      {
        //There is another site, so nothing to do...
        test = false;
      }
    }
    if (test)
    {
      //There was no other site pointing toward this pattern, so we remove it.
      delete sites_[current];
      sites_[current] = copy;
    }
    else
    {
      //We add the site at the end:
      sites_.push_back(copy);
      index_[pos] = siteIndex;
    }
  }
}

/******************************************************************************/

Site* CompressedVectorSiteContainer::removeSite(size_t i) throw (IndexOutOfBoundsException)
{
  if (i >= getNumberOfSites()) throw IndexOutOfBoundsException("CompressedVectorSiteContainer::removeSite.", i, 0, getNumberOfSites() - 1);
  //Here we return a copy of the site, as it will not necessarily be removed from the set, so we don't want to delete it.
  Site* site = dynamic_cast<Site *>(sites_[index_[i]]->clone());
  deleteSite(i);
  return site;
}

/******************************************************************************/

void CompressedVectorSiteContainer::deleteSite(size_t siteIndex) throw (IndexOutOfBoundsException)
{
  if (siteIndex >= getNumberOfSites())
    throw IndexOutOfBoundsException("CompressedVectorSiteContainer::deleteSite.", siteIndex, 0, getNumberOfSites() - 1);
  //Here we need to check whether the pattern corresponding to this site is unique:
  size_t current = index_[siteIndex];
  bool test = true;
  for (size_t j = 0; test && j < index_.size(); ++j)
  {
    if (j != siteIndex && index_[j] == current)
    {
      //There is a nother site, so nothing to...
      test = false;
    }
  }
  if (test)
  {
    //There was no other site pointing toward this pattern, so we remove it.
    delete sites_[current];
    sites_.erase(sites_.begin() + static_cast<ptrdiff_t>(current));
    //Now we have to correct all indices:
    for (size_t j = 0; j < index_.size(); ++j)
    {
      if (index_[j] > current) index_[j]--;
    }
  }
  index_.erase(index_.begin() + static_cast<ptrdiff_t>(siteIndex));
}

/******************************************************************************/

void CompressedVectorSiteContainer::deleteSites(size_t siteIndex, size_t length) throw (IndexOutOfBoundsException)
{
  //This may be optimized later:
  for (size_t i = 0; i < length; ++i) {
    deleteSite(siteIndex + i);
  }
}

/******************************************************************************/

void CompressedVectorSiteContainer::addSite(const Site& site, bool checkPositions) throw (Exception)
{
  // Check size:
  if (site.size() != getNumberOfSequences()) throw SiteException("CompressedVectorSiteContainer::addSite. Site does not have the appropriate length", &site);

  // New site's alphabet and site container's alphabet matching verification
  if (site.getAlphabet()->getAlphabetType() != getAlphabet()->getAlphabetType())
  {
    throw AlphabetMismatchException("CompressedVectorSiteContainer::addSite", getAlphabet(), site.getAlphabet());
  }

  size_t siteIndex = getSiteIndex_(site);
  if (siteIndex == sites_.size())
  {
    //This is a new pattern:
    Site* copy = dynamic_cast<Site*>(site.clone());
    sites_.push_back(copy);
  }
  index_.push_back(siteIndex);
}

/******************************************************************************/

void CompressedVectorSiteContainer::addSite(const Site& site, size_t siteIndex, bool checkPositions) throw (Exception)
{
  if (siteIndex >= getNumberOfSites()) throw IndexOutOfBoundsException("CompressedVectorSiteContainer::addSite", siteIndex, 0, getNumberOfSites() - 1);

  // Check size:
  if (site.size() != getNumberOfSequences()) throw SiteException("CompressedVectorSiteContainer::addSite. Site does not have the appropriate length", &site);

  // New site's alphabet and site container's alphabet matching verification
  if (site.getAlphabet()->getAlphabetType() != getAlphabet()->getAlphabetType())
  {
    throw AlphabetMismatchException("CompressedVectorSiteContainer::addSite", getAlphabet(), site.getAlphabet());
  }

  size_t index = getSiteIndex_(site);
  if (index == sites_.size())
  {
    //This is a new pattern:
    Site* copy = dynamic_cast<Site*>(site.clone());
    sites_.push_back(copy);
  }
  index_.insert(index_.begin() + static_cast<ptrdiff_t>(siteIndex), index);
}

/******************************************************************************/

void CompressedVectorSiteContainer::reindexSites()
{
  int pos = 1; // first position is 1.
  for (vector<Site*>::iterator i = sites_.begin(); i < sites_.end(); i++)
  {
    (*i)->setPosition(pos++);
  }
}

/******************************************************************************/

Vint CompressedVectorSiteContainer::getSitePositions() const
{
  size_t n = getNumberOfSites();
  Vint positions(n);
  for (size_t i = 0; i < n; i++)
  {
    positions[i] = sites_[index_[i]]->getPosition();
  }
  return positions;
}

/******************************************************************************/

const Sequence& CompressedVectorSiteContainer::getSequence(size_t i) const throw (IndexOutOfBoundsException)
{
  if (i >= getNumberOfSequences()) throw IndexOutOfBoundsException("CompressedVectorSiteContainer::getSequence.", i, 0, getNumberOfSequences() - 1);

  // Main loop : for all sites
  size_t n = getNumberOfSites();
  vector<int> sequence(n);
  for (size_t j = 0; j < n; j++)
  {
    sequence[j] = (*sites_[index_[j]])[i];
  }
  if (sequences_[i]) delete sequences_[i];
  sequences_[i] = new BasicSequence(names_[i], sequence, *comments_[i], getAlphabet());
  return *sequences_[i];
}

/******************************************************************************/

const Sequence& CompressedVectorSiteContainer::getSequence(const std::string& name) const throw (SequenceNotFoundException)
{
  // Look for sequence name:
  size_t pos = getSequencePosition(name);
  return getSequence(pos);
}

/******************************************************************************/

bool CompressedVectorSiteContainer::hasSequence(const string& name) const
{
  //Look for sequence name:
  for (size_t pos = 0; pos < names_.size(); pos++) {
    if (names_[pos] == name) return true;
  }
  return false;
}

/******************************************************************************/

size_t CompressedVectorSiteContainer::getSequencePosition(const std::string& name) const throw (SequenceNotFoundException)
{
  // Look for sequence name:
  for (size_t pos = 0; pos < names_.size(); pos++)
  {
    if (names_[pos] == name) return pos;
  }
  throw SequenceNotFoundException("CompressedVectorSiteContainer::getSequencePosition().", name);
}

/******************************************************************************/

void CompressedVectorSiteContainer::clear()
{
  // Must delete all sites in the container:
  for (size_t i = 0; i < sites_.size(); i++)
  {
    delete sites_[i];
  }

  // must delete all comments too:
  for (size_t i = 0; i < comments_.size(); i++)
  {
    if (comments_[i]) delete comments_[i];
  }

  // Delete all sequences retrieved:
  for (size_t i = 0; i < sequences_.size(); i++)
  {
    if (sequences_[i]) delete (sequences_[i]);
  }

  // Delete all sites pointers
  sites_.clear();
  index_.clear();
  names_.clear();
  comments_.clear();
  sequences_.clear();
}

/******************************************************************************/

vector<string> CompressedVectorSiteContainer::getSequencesNames() const
{
   vector<string> seqnames(names_.size());
  for (size_t i = 0; i < names_.size(); i++)
  {
    seqnames[i] = names_[i];
  }
  return seqnames;
}

/******************************************************************************/

void CompressedVectorSiteContainer::setSequencesNames(
  const vector<string>& names,
  bool checkNames)
throw (Exception)
{
  if (names.size() != getNumberOfSequences())
    throw IndexOutOfBoundsException("CompressedVectorSiteContainer::setSequenceNames: bad number of names.", names.size(), getNumberOfSequences(), getNumberOfSequences());
  if (checkNames)
  {
    for (size_t i = 0; i < names.size(); i++)
    {
      // For all names in vector : throw exception if name already exists
      for (size_t j = 0; j < i; j++)
      {
        if (names[j] == names[i])
          throw Exception("CompressedVectorSiteContainer::setSequencesNames : Sequence's name already exists in container");
      }
    }
  }
  for (size_t i = 0; i < names.size(); i++)
  {
    names_[i] = names[i];
  }
}

/******************************************************************************/

void CompressedVectorSiteContainer::setComments(size_t sequenceIndex, const Comments& comments) throw (IndexOutOfBoundsException)
{
  comments_[sequenceIndex] = new Comments(comments);
}

/******************************************************************************/

CompressedVectorSiteContainer* CompressedVectorSiteContainer::createEmptyContainer() const
{
   CompressedVectorSiteContainer* vsc = new CompressedVectorSiteContainer(getAlphabet());
   vsc->setGeneralComments(getGeneralComments());
  return vsc;
}

/******************************************************************************/

size_t CompressedVectorSiteContainer::getSiteIndex_(const Site& site)
{
  size_t pos = sites_.size();
  bool test;
  for (size_t i = 0; i < sites_.size(); ++i)
  {
    test = true;
    for (size_t j = 0; test && j < site.size(); ++j) //site is supposed to have the correct size, that is the same as all the ones in the container.
    {
      if (site[j] != (*sites_[i])[j])
        test = false;
    }
    if (test)
    {
      pos = i;
      break;
    }
  }
  return pos;
}

/******************************************************************************/

