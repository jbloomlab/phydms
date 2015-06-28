//
// File: SiteContainerIterator.cpp
// Created by: Julien Dutheil
// Created on: Sun Oct 19 12:47:16 2003
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

#include "SiteContainerIterator.h"
#include "../SiteTools.h"

using namespace bpp;

// From the STL:
#include <iostream>

using namespace std;

/******************************************************************************/
	
AbstractSiteContainerIterator::AbstractSiteContainerIterator(const SiteContainer& sites) :
  sites_(&sites),
  currentPosition_(0)
{}

/******************************************************************************/
	
SimpleSiteContainerIterator::SimpleSiteContainerIterator(const SiteContainer& sites): AbstractSiteContainerIterator(sites) {}

const Site* SimpleSiteContainerIterator::nextSite()
{
	const Site* s = &sites_->getSite(static_cast<size_t>(currentPosition_));
	currentPosition_++;
	return s;
}

bool SimpleSiteContainerIterator::hasMoreSites() const
{
	return currentPosition_ < static_cast<int>(sites_->getNumberOfSites());
}

/******************************************************************************/
	
NoGapSiteContainerIterator::NoGapSiteContainerIterator(const SiteContainer& sites): AbstractSiteContainerIterator(sites)
{
	currentPosition_ = nextSiteWithoutGapPosition(-1);
}

const Site* NoGapSiteContainerIterator::nextSite()
{
	const Site* s = &sites_->getSite(static_cast<size_t>(currentPosition_));
	currentPosition_ = nextSiteWithoutGapPosition(currentPosition_);
	return s;
}

bool NoGapSiteContainerIterator::hasMoreSites() const
{
	return currentPosition_ < static_cast<int>(sites_->getNumberOfSites());
}

int NoGapSiteContainerIterator::nextSiteWithoutGapPosition(int current) const
{
	size_t position = static_cast<size_t>(current + 1);
	while (position < sites_->getNumberOfSites() && SiteTools::hasGap(sites_->getSite(position)))
    position++;
	return static_cast<int>(position);
}

int NoGapSiteContainerIterator::previousSiteWithoutGapPosition(int current) const
{
	int position = current - 1;
	while (position >= 0 && SiteTools::hasGap(sites_->getSite(static_cast<size_t>(position))))
    position--;
	return position;
}

/******************************************************************************/
	
CompleteSiteContainerIterator::CompleteSiteContainerIterator(const SiteContainer & sites): AbstractSiteContainerIterator(sites)
{
	currentPosition_ = nextCompleteSitePosition(-1);
}

const Site* CompleteSiteContainerIterator::nextSite()
{
	const Site* s = &sites_->getSite(static_cast<size_t>(currentPosition_));
	currentPosition_ = nextCompleteSitePosition(currentPosition_);
	return s;
}

bool CompleteSiteContainerIterator::hasMoreSites() const
{
	return currentPosition_ < static_cast<int>(sites_->getNumberOfSites());
}

int CompleteSiteContainerIterator::nextCompleteSitePosition(int current) const
{
  size_t position = static_cast<size_t>(current + 1);
	while (position < sites_->getNumberOfSites() && !SiteTools::isComplete(sites_->getSite(position)))
    position++;
	return static_cast<int>(position);
}

int CompleteSiteContainerIterator::previousCompleteSitePosition(int current) const
{
  int position = current - 1;
	while (position >= 0 && !SiteTools::isComplete(sites_->getSite(static_cast<size_t>(position))))
    position --;
	return position;
}

/******************************************************************************/

