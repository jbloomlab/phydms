//
// File: VectorProbabilisticSiteContainer.cpp
// Created by: Murray Patterson
// Created on: Mon Oct 19 2015
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

#include "VectorProbabilisticSiteContainer.h"

#include <Bpp/Text/TextTools.h>
#include <Bpp/Numeric/DataTable.h>

using namespace bpp;

/********************************************************************************/

VectorProbabilisticSiteContainer::VectorProbabilisticSiteContainer(const Alphabet * alpha) :
  VectorSiteContainer(alpha),
  p_sites_(0),
  p_sequences_(0)
{}

/********************************************************************************/

const ProbabilisticSite & VectorProbabilisticSiteContainer::getProbabilisticSite(std::size_t i) const throw (IndexOutOfBoundsException)
{
  if(i >= getNumberOfProbabilisticSites())
    throw IndexOutOfBoundsException("VectorProbabilisticSiteContainer::getProbabilisticSite.", i, 0, getNumberOfProbabilisticSites() - 1);

  return *p_sites_[i];
}

/********************************************************************************/

void VectorProbabilisticSiteContainer::addSite(const ProbabilisticSite & site, bool checkPosition) throw (Exception)
{
  // check size :
  if(site.size() != getNumberOfProbabilisticSequences())
    throw Exception("VectorProbabilisticSiteContainer::addSite. Site does not have the appropriate length: " + TextTools::toString(site.size()) + ", should be " + TextTools::toString(getNumberOfProbabilisticSequences()) + ".");

  // new site's alphabet and site container's alphabet must match :
  if(site.getAlphabet()->getAlphabetType() != getAlphabet()->getAlphabetType())
    throw AlphabetMismatchException("VectorProbabilisticSiteContainer::addSite.", getAlphabet(), site.getAlphabet());

  // check position :
  if(checkPosition) {

    int position = site.getPosition();
    // for all positions in vector : throw exception if position already exists
    for(std::size_t i = 0; i < p_sites_.size(); ++i)
      if(p_sites_[i]->getPosition() == position)
	throw Exception("VectorSiteContainer::addSite. Site position: " + TextTools::toString(position) + ", already exists in container.");
  }

  p_sites_.push_back(dynamic_cast<ProbabilisticSite *>(site.clone()));
}

/********************************************************************************/

const ProbabilisticSequence & VectorProbabilisticSiteContainer::getProbabilisticSequence(std::size_t i) const throw (IndexOutOfBoundsException)
{

  if(i >= getNumberOfProbabilisticSequences())
    throw IndexOutOfBoundsException("VectorProbabilisticSiteContainer::getProbabilisticSequence.", i, 0, getNumberOfProbabilisticSequences() - 1);

  // main loop : for all sites
  std::size_t n = getNumberOfProbabilisticSites();
  DataTable sequence(getAlphabet()->getResolvedChars());
  for(std::size_t j = 0; j < n; ++j)
    sequence.addRow(p_sites_[j]->getContent().getRow(i));

  if(p_sequences_[i])
    delete p_sequences_[i];

  p_sequences_[i] = new BasicProbabilisticSequence(names_[i], sequence, *comments_[i], getAlphabet());

  return *p_sequences_[i];
}

/********************************************************************************/

void VectorProbabilisticSiteContainer::addSequence(const ProbabilisticSequence & sequence, bool checkName) throw (Exception)
{

  // if the container has no sequence, we set the size to the size of this sequence :
  if(getNumberOfProbabilisticSequences() == 0)
    pRealloc(sequence.size());

  // new sequence's alphabet and site container's alphabet must match :
  if(sequence.getAlphabet()->getAlphabetType() != getAlphabet()->getAlphabetType())
    throw AlphabetMismatchException("VectorProbabilisticSiteContainer::addSequence.", getAlphabet(), sequence.getAlphabet());

  if(sequence.size() != p_sites_.size())
    throw Exception("VectorProbabilisticSiteContainer::addSequence. Sequence does not have the appropriate length: " + TextTools::toString(sequence.size()) + ", should be " + TextTools::toString(p_sites_.size()) + ".");

  // check name :
  if(checkName)
    for(std::size_t i = 0; i < names_.size(); ++i)
      if(sequence.getName() == names_[i])
	throw Exception("VectorProbabilisticSiteContainer::addSequence. Name: " + sequence.getName() + ", already exists in the container.");

  // append name :
  names_.push_back(sequence.getName());

  // append elements at each site :
  for(size_t i = 0; i < p_sites_.size(); ++i)
    p_sites_[i]->addElement(sequence.getContent().getRow(i));

  // append comments :
  comments_.push_back(new Comments(sequence.getComments()));

  // sequence pointers :
  p_sequences_.push_back(0);
}

/********************************************************************************/

void VectorProbabilisticSiteContainer::pClear()
{
  clear(); // call VectorSiteContainer clear

  // now clear all probabilistic sites / sequences
  for(std::size_t i = 0; i < p_sites_.size(); ++i)
    delete p_sites_[i];

  for(std::size_t i = 0; i < p_sequences_.size(); ++i)
    delete p_sequences_[i];

  // and delete the corresponding pointers
  p_sites_.clear();
  p_sequences_.clear();
}  

/********************************************************************************/

void VectorProbabilisticSiteContainer::reindexpSites()
{
  int pos = 1; // start at position 1
  std::vector<ProbabilisticSite *>::iterator i = p_sites_.begin();
  for(; i != p_sites_.end(); ++i)
    (*i)->setPosition(++pos);
}

/********************************************************************************/

void VectorProbabilisticSiteContainer::pRealloc(std::size_t n)
{
  pClear();
  p_sites_.resize(n);

  for(std::size_t i = 0; i < n; ++i)
    p_sites_[i] = new BasicProbabilisticSite(getAlphabet());

  reindexpSites();
}
