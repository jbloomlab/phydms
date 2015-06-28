//
// File: SymbolList.cpp
// Created by: Julien Dutheil
// Created on: Fri Apr 9 2005
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

#include "SymbolList.h"
#include "StringSequenceTools.h"

using namespace bpp;

using namespace std;

/****************************************************************************************/

BasicSymbolList::BasicSymbolList(const std::vector<string>& list, const Alphabet* alpha) throw (BadCharException) :
	alphabet_(alpha), content_()
{
	setContent(list);
}

BasicSymbolList::BasicSymbolList(const std::vector<int>& list, const Alphabet* alpha) throw (BadIntException) :
	alphabet_(alpha), content_()
{
	setContent(list);
}

/****************************************************************************************/

BasicSymbolList::BasicSymbolList(const SymbolList& list):
  alphabet_(list.getAlphabet()), content_(list.size())
{
  for (size_t i = 0; i < list.size(); ++i)
    content_[i] = list[i];
}

BasicSymbolList::BasicSymbolList(const BasicSymbolList& list):
  alphabet_(list.alphabet_), content_(list.content_) {}

BasicSymbolList& BasicSymbolList::operator=(const SymbolList& list)
{
	content_.resize(list.size());
  for (size_t i = 0; i < list.size(); ++i)
    content_[i] = list[i];
	alphabet_ = list.getAlphabet();
	return *this;
}

BasicSymbolList& BasicSymbolList::operator=(const BasicSymbolList& list)
{
	content_  = list.content_;
	alphabet_ = list.alphabet_;
	return *this;
}

/****************************************************************************************/

void BasicSymbolList::setContent(const vector<string>& list) throw (BadCharException)
{
	// Check list for incorrect characters
	vector<int> coded(list.size());
	for (size_t i = 0; i < list.size(); i++)
		if(!alphabet_->isCharInAlphabet(list[i])) throw BadCharException(list[i], "BasicSymbolList::setContent", alphabet_);

  for (size_t i = 0; i < list.size(); i++) 
		coded[i] = alphabet_->charToInt(list[i]);
	
  //BasicSymbolList is valid:
	content_ = coded;
};

/****************************************************************************************/

void BasicSymbolList::setContent(const vector<int>& list) throw (BadIntException)
{
	// Check list for incorrect characters
	for (size_t i = 0; i < list.size(); i++)
		if(!alphabet_->isIntInAlphabet(list[i]))
      throw BadIntException(list[i], "BasicSymbolList::setContent", alphabet_);
	
  //Sequence is valid:
	content_ = list;
};

/****************************************************************************************/

string BasicSymbolList::toString() const
{
	return StringSequenceTools::decodeSequence(content_, alphabet_);
};

/****************************************************************************************/

void BasicSymbolList::addElement(const string& c) throw (BadCharException)
{
	content_.push_back(alphabet_->charToInt(c));
}

/****************************************************************************************/

void BasicSymbolList::addElement(size_t pos, const string& c) throw (BadCharException, IndexOutOfBoundsException)
{
  if(pos >= content_.size()) throw IndexOutOfBoundsException("BasicSymbolList::addElement. Invalid position.", pos, 0, size() - 1);
	content_.insert(content_.begin() + static_cast<ptrdiff_t>(pos), alphabet_->charToInt(c));
}

/****************************************************************************************/

void BasicSymbolList::setElement(size_t pos, const string& c) throw (BadCharException, IndexOutOfBoundsException)
{
	if(pos >= content_.size())
    throw IndexOutOfBoundsException("BasicSymbolList::setElement. Invalid position.", pos, 0, size() - 1);
	content_[pos] = alphabet_->charToInt(c);
}

/****************************************************************************************/

string BasicSymbolList::getChar(size_t pos) const throw (IndexOutOfBoundsException)
{
	if(pos >= content_.size())
    throw IndexOutOfBoundsException("BasicSymbolList::getChar. Invalid position.", pos, 0, size() - 1);
	string c = "";
	try
  {
		c = alphabet_->intToChar(content_[pos]);
	}
  catch(BadIntException bie)
  {
		//This should never happen!
	}
	return c;
}

/****************************************************************************************/

void BasicSymbolList::deleteElement(size_t pos) throw (IndexOutOfBoundsException)
{
	if(pos >= content_.size())
    throw IndexOutOfBoundsException("BasicSymbolList::deleteElement. Invalid position.", pos, 0, size() - 1);
	content_.erase(content_.begin() + static_cast<ptrdiff_t>(pos));
}

/****************************************************************************************/

void BasicSymbolList::deleteElements(size_t pos, size_t len) throw (IndexOutOfBoundsException)
{
	if (pos + len > content_.size())
    throw IndexOutOfBoundsException("BasicSymbolList::deleteElements. Invalid position.", pos + len, 0, size() - 1);
	 content_.erase(content_.begin() + static_cast<ptrdiff_t>(pos), content_.begin() + static_cast<ptrdiff_t>(pos + len));
}

/****************************************************************************************/

void BasicSymbolList::addElement(int v) throw (BadIntException)
{
	//test:
	alphabet_->intToChar(v);
	content_.push_back(v);
}

/****************************************************************************************/

void BasicSymbolList::addElement(size_t pos, int v) throw (BadIntException, IndexOutOfBoundsException)
{
	//test:
	if(pos >= content_.size())
    throw IndexOutOfBoundsException("BasicSymbolList::addElement. Invalid position.", pos, 0, size() - 1);
	alphabet_->intToChar(v);
	content_.insert(content_.begin() + static_cast<ptrdiff_t>(pos), v);
}

/****************************************************************************************/

void BasicSymbolList::setElement(size_t pos, int v) throw (BadIntException, IndexOutOfBoundsException)
{
	//test:
  if(pos >= content_.size())
    throw IndexOutOfBoundsException("BasicSymbolList::setElement. Invalid position.", pos, 0, size() - 1);
	alphabet_->intToChar(v);
	content_[pos] = v;
}

/****************************************************************************************/

int BasicSymbolList::getValue(size_t pos) const throw (IndexOutOfBoundsException)
{
  if(pos >= content_.size())
    throw IndexOutOfBoundsException("BasicSymbolList::getValue. Invalid position.", pos, 0, size() - 1);
	return content_[pos];
}

/****************************************************************************************/


/****************************************************************************************/

EdSymbolList::EdSymbolList(const std::vector<string>& list, const Alphabet* alpha) throw (BadCharException) :
	alphabet_(alpha), propagateEvents_(true), content_(), listeners_()
{
	setContent(list);
}

EdSymbolList::EdSymbolList(const std::vector<int>& list, const Alphabet* alpha) throw (BadIntException) :
	alphabet_(alpha), propagateEvents_(true), content_(), listeners_()
{
	setContent(list);
}

/****************************************************************************************/

EdSymbolList::EdSymbolList(const SymbolList& list):
  alphabet_(list.getAlphabet()), propagateEvents_(true), content_(list.size()), listeners_()
{
  for (size_t i = 0; i < list.size(); ++i) {
    content_[i] = list[i];
  }
}

EdSymbolList::EdSymbolList(const EdSymbolList& list):
  alphabet_(list.getAlphabet()), propagateEvents_(list.propagateEvents_), content_(list.content_), listeners_(list.listeners_)
{
  for (size_t i = 0; i < listeners_.size(); ++i)
    if (!list.listeners_[i]->isShared())
      listeners_[i] = dynamic_cast<SymbolListListener*>(list.listeners_[i]->clone());
}

EdSymbolList& EdSymbolList::operator=(const SymbolList& list)
{
	content_.resize(list.size());
  for (size_t i = 0; i < list.size(); ++i) {
    content_[i] = list[i];
  }
	alphabet_        = list.getAlphabet();
  propagateEvents_ = true;
  for (size_t i = 0; i < listeners_.size(); ++i)
    if (!listeners_[i]->isShared())
     delete listeners_[i];
  listeners_.clear();
	return *this;
}

EdSymbolList& EdSymbolList::operator=(const EdSymbolList& list)
{
	content_         = list.getContent();
	alphabet_        = list.getAlphabet();
  propagateEvents_ = list.propagateEvents_;
  for (size_t i = 0; i < listeners_.size(); ++i)
    delete listeners_[i];
  listeners_ = list.listeners_;
  for (size_t i = 0; i < listeners_.size(); ++i)
    if (!list.listeners_[i]->isShared())
      listeners_[i] = dynamic_cast<SymbolListListener*>(list.listeners_[i]->clone());
	return *this;
}

/****************************************************************************************/

void EdSymbolList::setContent(const vector<string>& list) throw (BadCharException)
{
  SymbolListEditionEvent event(this);
  fireBeforeSequenceChanged(event);

  // Check list for incorrect characters
	vector<int> coded(list.size());
	for (size_t i = 0; i < list.size(); i++)
		if (!alphabet_->isCharInAlphabet(list[i])) throw BadCharException(list[i], "EdSymbolList::setContent", alphabet_);
	
  for (size_t i = 0; i < list.size(); i++) 
		coded[i] = alphabet_->charToInt(list[i]);
	
  //SymbolList is valid:
	content_ = coded;
  fireAfterSequenceChanged(event);
};

/****************************************************************************************/

void EdSymbolList::setContent(const vector<int>& list) throw (BadIntException)
{
  SymbolListEditionEvent event(this);
  fireBeforeSequenceChanged(event);

	// Check list for incorrect characters
	for (size_t i = 0; i < list.size(); i++)
		if(!alphabet_->isIntInAlphabet(list[i]))
      throw BadIntException(list[i], "EdSymbolList::setContent", alphabet_);
	
  //Sequence is valid:
	content_ = list;
  fireAfterSequenceChanged(event);
};

/****************************************************************************************/

string EdSymbolList::toString() const
{
	return StringSequenceTools::decodeSequence(content_, alphabet_);
};

/****************************************************************************************/

void EdSymbolList::addElement(const string& c) throw (BadCharException)
{
  SymbolListInsertionEvent event(this, size(), 1);
  fireBeforeSequenceInserted(event);
	content_.push_back(alphabet_->charToInt(c));
  fireAfterSequenceInserted(event);
}

/****************************************************************************************/

void EdSymbolList::addElement(size_t pos, const string& c) throw (BadCharException, IndexOutOfBoundsException)
{
  if (pos >= content_.size()) throw IndexOutOfBoundsException("EdSymbolList::addElement. Invalid position.", pos, 0, size() - 1);
  SymbolListInsertionEvent event(this, pos, 1);
  fireBeforeSequenceInserted(event);
	content_.insert(content_.begin() + static_cast<ptrdiff_t>(pos), alphabet_->charToInt(c));
  fireAfterSequenceInserted(event);
}

/****************************************************************************************/

void EdSymbolList::setElement(size_t pos, const string& c) throw (BadCharException, IndexOutOfBoundsException)
{
	if (pos >= content_.size())
    throw IndexOutOfBoundsException("EdSymbolList::setElement. Invalid position.", pos, 0, size() - 1);
  SymbolListSubstitutionEvent event(this, pos, pos);
  fireBeforeSequenceSubstituted(event);
	content_[pos] = alphabet_->charToInt(c);
  fireAfterSequenceSubstituted(event);
}

/****************************************************************************************/

string EdSymbolList::getChar(size_t pos) const throw (IndexOutOfBoundsException)
{
	if (pos >= content_.size())
    throw IndexOutOfBoundsException("EdSymbolList::getChar. Invalid position.", pos, 0, size() - 1);
	string c = "";
	try {
		c = alphabet_->intToChar(content_[pos]);
	} catch(BadIntException bie) {
		//This should never happen!
	}
	return c;
}

/****************************************************************************************/

void EdSymbolList::deleteElement(size_t pos) throw (IndexOutOfBoundsException)
{
	if (pos >= content_.size())
    throw IndexOutOfBoundsException("EdSymbolList::deleteElement. Invalid position.", pos, 0, size() - 1);
  SymbolListDeletionEvent event(this, pos, 1);
  fireBeforeSequenceDeleted(event);
	content_.erase(content_.begin() + static_cast<ptrdiff_t>(pos));
  fireAfterSequenceDeleted(event);
}

/****************************************************************************************/

void EdSymbolList::deleteElements(size_t pos, size_t len) throw (IndexOutOfBoundsException)
{
	if(pos + len > content_.size())
    throw IndexOutOfBoundsException("EdSymbolList::deleteElements. Invalid position.", pos + len, 0, size() - 1);
  SymbolListDeletionEvent event(this, pos, len);
  fireBeforeSequenceDeleted(event);
	content_.erase(content_.begin() + static_cast<ptrdiff_t>(pos), content_.begin() + static_cast<ptrdiff_t>(pos + len));
  fireAfterSequenceDeleted(event);
}


/****************************************************************************************/

void EdSymbolList::addElement(int v) throw (BadIntException)
{
  SymbolListInsertionEvent event(this, size(), 1);
  fireBeforeSequenceInserted(event);
	//test:
	alphabet_->intToChar(v);
	content_.push_back(v);
  fireAfterSequenceInserted(event);
}

/****************************************************************************************/

void EdSymbolList::addElement(size_t pos, int v) throw (BadIntException, IndexOutOfBoundsException)
{
	//test:
	if (pos >= content_.size())
    throw IndexOutOfBoundsException("EdSymbolList::addElement. Invalid position.", pos, 0, size() - 1);
  SymbolListInsertionEvent event(this, pos, 1);
  fireBeforeSequenceInserted(event);
	alphabet_->intToChar(v);
	content_.insert(content_.begin() + static_cast<ptrdiff_t>(pos), v);
  fireAfterSequenceInserted(event);
}

/****************************************************************************************/

void EdSymbolList::setElement(size_t pos, int v) throw (BadIntException, IndexOutOfBoundsException)
{
	//test:
  if (pos >= content_.size())
    throw IndexOutOfBoundsException("EdSymbolList::setElement. Invalid position.", pos, 0, size() - 1);
  SymbolListSubstitutionEvent event(this, pos, pos);
  fireBeforeSequenceSubstituted(event);
	alphabet_->intToChar(v);
	content_[pos] = v;
  fireAfterSequenceSubstituted(event);
}

/****************************************************************************************/

int EdSymbolList::getValue(size_t pos) const throw (IndexOutOfBoundsException)
{
  if (pos >= content_.size())
    throw IndexOutOfBoundsException("EdSymbolList::getValue. Invalid position.", pos, 0, size() - 1);
	return content_[pos];
}

/****************************************************************************************/

