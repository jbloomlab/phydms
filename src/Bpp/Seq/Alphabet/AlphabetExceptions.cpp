//
// File: AlphabetExceptions.cpp
// Created by: Julien Dutheil
// Created on: Mon Nov  3 16:41:53 2003
//

/*
Copyright or © or Copr. CNRS, (November 17, 2004)

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

#include "AlphabetExceptions.h"
#include "Alphabet.h"
#include <Bpp/Text/TextTools.h>

using namespace bpp;

using namespace std;

/******************************************************************************
 *                         Alphabet exceptions:                               *
 ******************************************************************************/

AlphabetException::AlphabetException(const std::string& text, const Alphabet* alpha) :
	Exception("AlphabetException: " + text + (alpha ? "(" + (alpha->getAlphabetType()) + ")" : string(""))),
	alphabet_(alpha) {}
		
/******************************************************************************/

BadCharException::BadCharException(const std::string& badChar, const std::string& text, const Alphabet* alpha) :
	AlphabetException("BadCharException: " + badChar + ". " + text, alpha),
	c_(badChar) {}
		
string BadCharException::getBadChar() const { return c_; }

/******************************************************************************/

BadIntException::BadIntException(int badInt, const std::string& text, const Alphabet* alpha) :
	AlphabetException("BadIntException: " + TextTools::toString(badInt) + ". " + text, alpha),
	i_(badInt) {}
		
int BadIntException::getBadInt() const { return i_; }

/******************************************************************************/
		
AlphabetMismatchException::AlphabetMismatchException(const std::string& text, const Alphabet* alpha1, const Alphabet* alpha2) :
	Exception("AlphabetMismatchException: " + text + (alpha1 != 0 && alpha2 != 0 ? "(" + alpha1->getAlphabetType() + ", " + alpha2->getAlphabetType() + ")" : string(""))),
	alphabet1_(alpha1),
	alphabet2_(alpha2) {}
		
vector<const Alphabet*> AlphabetMismatchException::getAlphabets() const
{
	vector<const Alphabet*> v(2);
	v[0] = alphabet1_;
	v[1] = alphabet2_;
	return v;
}

/******************************************************************************/

CharStateNotSupportedException::CharStateNotSupportedException(const string & text, const Alphabet * alpha) :
  AlphabetException("CharStateNotSupportedException: " + text, alpha) {};

/******************************************************************************/
