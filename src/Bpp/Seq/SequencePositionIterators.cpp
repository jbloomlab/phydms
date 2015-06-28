// 
// File:    SequencePositionIterators.cpp
// Author:  Sylvain Gaillard
// Created: 23/06/2009 11:38:27
// 

/*
Copyright or © or Copr. Bio++ Development Team, (June 23, 2009)

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

#include "SequencePositionIterators.h"

using namespace bpp;
using namespace std; // for the STL

/******************************************************************************/

bool AbstractSequencePositionIterator::operator==(const SequencePositionIterator& it) const {
  return this->getPosition() == it.getPosition();
}

/******************************************************************************/

bool AbstractSequencePositionIterator::operator!=(const SequencePositionIterator& it) const {
  return this->getPosition() != it.getPosition();
}

/******************************************************************************/

void AbstractSequencePositionIterator::setPosition(unsigned int pos) {
  this->currentPosition_ = pos;
}

/******************************************************************************/

const Sequence & AbstractSequencePositionIterator::getSequence() const {
  return * (this->sequence_);
}

/******************************************************************************/

unsigned int AbstractSequencePositionIterator::getPosition() const {
  return this->currentPosition_;
}

/******************************************************************************/

int AbstractSequencePositionIterator::getValue() const {
  return this->sequence_->getValue(this->currentPosition_);
}

/******************************************************************************/

string AbstractSequencePositionIterator::getChar() const {
  return this->sequence_->getChar(this->currentPosition_);
}


//===============================
// SimpleSequencePositionIterator
//===============================
/******************************************************************************/

SimpleSequencePositionIterator::SimpleSequencePositionIterator(const SequencePositionIterator& it):
  AbstractSequencePositionIterator(it.getSequence(), it.getPosition()) {};

/******************************************************************************/

SimpleSequencePositionIterator& SimpleSequencePositionIterator::operator++() {
  this->setPosition(this->getPosition() + 1);
  return *this;
}

/******************************************************************************/

SimpleSequencePositionIterator SimpleSequencePositionIterator::operator++(int i) {
  SimpleSequencePositionIterator ans = *this;
  ++(*this);
  return ans;
}

/******************************************************************************/

SimpleSequencePositionIterator& SimpleSequencePositionIterator::operator+=(int i) {
  if (i > 0)
    this->setPosition(this->getPosition() + static_cast<unsigned int>(i));
  else if (i < 0) {
    unsigned int d = static_cast<unsigned int>(-i);
    if (d > this->getPosition())
      throw Exception("SimpleSequencePositionIterator::operator+=. Negative increment too large.");
    else
      this->setPosition(this->getPosition() - d);
  }
  return *this;
}

/******************************************************************************/

SimpleSequencePositionIterator& SimpleSequencePositionIterator::operator-=(int i) {
  return (*this) += -i;
}

/******************************************************************************/

SimpleSequencePositionIterator SimpleSequencePositionIterator::operator+(int i) const {
  SimpleSequencePositionIterator res(*this);
  res += i;
  return res;
}

/******************************************************************************/

SimpleSequencePositionIterator SimpleSequencePositionIterator::operator-(int i) const {
  return (*this) + (- i);
}

/******************************************************************************/

bool SimpleSequencePositionIterator::hasMorePositions() const {
  return (this->getPosition() < this->getSequence().size());
}
/******************************************************************************/
