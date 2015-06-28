//
// File: StringSequenceTools.cpp
// Created by: Julien Dutheil
// Created on: Sun Nov 30 11:29:07 2003
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

#include "StringSequenceTools.h"

#include "Alphabet/AlphabetTools.h"
#include "Alphabet/DNA.h"
#include "Alphabet/RNA.h"
#include "Alphabet/ProteicAlphabet.h"
#include <Bpp/Text/TextTools.h>
#include <Bpp/Numeric/Random/RandomTools.h>

using namespace bpp;

// From the STL:
#include <map>
#include <ctype.h>
#include <algorithm>
#include <iostream>

using namespace std;

/****************************************************************************************/

string StringSequenceTools::subseq(const string& sequence, size_t begin, size_t end) throw (Exception)
{
  // Checking interval
  if (end < begin)
    throw Exception ("StringSequenceTools::subseq: Invalid interval");

  // Copy sequence
  string temp(sequence);

  // Truncate sequence
  temp.erase(temp.begin() + static_cast<ptrdiff_t>(end + 1), temp.end());
  temp.erase(temp.begin(), temp.begin() + static_cast<ptrdiff_t>(begin));

  // Send result
  return temp;
}

/****************************************************************************************/

string StringSequenceTools::setToSizeR(const string& sequence, size_t size)
{
  return TextTools::resizeRight(sequence, size, '-');
}

string StringSequenceTools::setToSizeL(const string& sequence, size_t size)
{
  return TextTools::resizeLeft(sequence, size, '-');
}

/****************************************************************************************/

string StringSequenceTools::deleteChar(const string& sequence, char chars)
{
  // Copy sequence
  string result(sequence);

  // Search and delete specified char
  for (unsigned int i = 0; i < result.size(); i++)
  {
    if (result[i] == chars)
      result.erase(result.begin() + i);
  }

  return result;
}

/****************************************************************************************/

string StringSequenceTools::deleteChar(const string& sequence, string chars)
{
  // Copy sequence
  string result(sequence);

  // For all characters to delete
  for (unsigned int i = 0; i < chars.size(); i++)
  {
    // Search and delete char
    for (unsigned int j = 0; j < result.size(); j++)
    {
      if (result[j] == chars[i])
        result.erase(result.begin() + j);
    }
  }

  return result;
}

/****************************************************************************************/

string* StringSequenceTools::reverse(const string& sequence)
{
  // Initializing
  string* result = new string;

  // Main loop : reverse all characters of sequence
  size_t size = sequence.size();
  for (size_t i = 0; i < size; i++)
  {
    *result += sequence[size - i - 1];
  }

  // Send result
  return result;
}

/****************************************************************************************/

string* StringSequenceTools::complement(const string& sequence)
{
  // Initializing
  string* result = new string;

  // Main loop : completement all characters
  size_t size = sequence.size();
  for (unsigned int i = 0; i < size; i++)
  {
    switch (sequence[i])
    {
    case 'A': *result += 'T';
      break;
    case 'C': *result += 'G';
      break;
    case 'G': *result += 'C';
      break;
    case 'T': *result += 'A';
      break;
    case 'M': *result += 'K';
      break;
    case 'R': *result += 'Y';
      break;
    case 'Y': *result += 'R';
      break;
    case 'K': *result += 'M';
      break;
    case 'V': *result += 'B';
      break;
    case 'H': *result += 'D';
      break;
    case 'D': *result += 'H';
      break;
    case 'B': *result += 'V';
      break;
    default: *result += sequence[i];
      break;
    }
  }

  // Send new sequence
  return result;
}

/****************************************************************************************/

double StringSequenceTools::getGCcontent(const string& sequence, size_t pos, size_t window) throw (BadIntegerException, Exception)
{
  // Frequency counts for nucleotids A, C, G, T
  map<char, double> counts;

  // Window size checking
  if (window < sequence.size())
    throw BadIntegerException("StringSequenceTools::getGCContent : specified window too high", static_cast<int>(window));

  // For last nucleotides
  if (pos + window > sequence.size())
  {
    pos = sequence.size() - window;
  }

  // Main loop
  for (size_t i = pos; i < pos + window; i++)
  {
    switch (toupper(sequence[i]))
    {
    case 'A': counts['A'] += 1;
      break;
    case 'C': counts['C'] += 1;
      break;
    case 'G': counts['G'] += 1;
      break;
    case 'T': counts['T'] += 1;
      break;
    case 'M': counts['A'] += 0.5;
      counts['C'] += 0.5;
      break;
    case 'R': counts['A'] += 0.5;
      counts['G'] += 0.5;
      break;
    case 'W': counts['A'] += 0.5;
      counts['T'] += 0.5;
      break;
    case 'S': counts['C'] += 0.5;
      counts['G'] += 0.5;
      break;
    case 'Y': counts['C'] += 0.5;
      counts['T'] += 0.5;
      break;
    case 'K': counts['G'] += 0.5;
      counts['T'] += 0.5;
      break;
    case 'V': counts['A'] += 0.34;
      counts['C'] += 0.34;
      counts['G'] += 0.34;
      break;
    case 'H': counts['A'] += 0.34;
      counts['C'] += 0.34;
      counts['T'] += 0.34;
      break;
    case 'D': counts['A'] += 0.34;
      counts['G'] += 0.34;
      counts['T'] += 0.34;
      break;
    case 'B': counts['C'] += 0.34;
      counts['G'] += 0.34;
      counts['T'] += 0.34;
      break;
    case '-': throw Exception("StringSequenceTools::getGCContent : Gap found in sequence");
      break;
    // Unresolved bases
    default: counts['A'] += 0.25;
      counts['C'] += 0.25;
      counts['G'] += 0.25;
      counts['T'] += 0.25;
    }
  }

  // Calculate and send GC rate
  return (counts['G'] + counts['C']) / static_cast<double>(window);
}

/****************************************************************************************/

vector<int> StringSequenceTools::codeSequence(const string& sequence, const Alphabet* alphabet)
throw (BadCharException)
{
  unsigned int size = AlphabetTools::getAlphabetCodingSize(alphabet); // Warning, an exception may be casted here!
  vector<int> code(static_cast<size_t>(floor(static_cast<double>(sequence.size()) / static_cast<double>(size))));
  size_t pos = 0;
  size_t count = 0;
  while (pos + size <= sequence.size())
  {
    code[count] = alphabet->charToInt(sequence.substr(pos, size));
    count++;
    pos += size;
  }
  return code;
}

/****************************************************************************************/

string StringSequenceTools::decodeSequence(const vector<int>& sequence, const Alphabet* alphabet) throw (BadIntException)
{
  string result = "";
  for (unsigned int i = 0; i < sequence.size(); i++)
  {
    result += alphabet->intToChar(sequence[i]);
  }
  return result;
}

/****************************************************************************************/

Alphabet* StringSequenceTools::getAlphabetFromSequence(const std::string& sequence)
throw (EmptySequenceException, SequenceException, AlphabetException)
{
  // empty sequence test
  if (sequence.size() == 0)
  {
    throw EmptySequenceException("Sequence::getAlphabetFromSequence : Empty sequence string");
  }

  // initialisation
  bool p = false; // indicates that a protein specific character is found
  bool r = false; // indicates that a RNA specific character is found
  bool u = false; // indicates that an unknown character is found
  bool pd = false; // Protein or DNA (T)

  // Main loop : for all character in sequence
  for (unsigned int i = 0; i < sequence.size(); i++)
  {
    // Character analyse
    switch (AlphabetTools::getType(sequence[i]))
    {
    case 0: u = true; break;
    case 3: p = true; break;
    case 2: r = true; break;
    case 5: pd = true; break;
    }
  }

  if (u)
    throw AlphabetException ("Sequence::getAlphabetFromSequence : Unknow character detected in specified sequence");
  if (r && pd)
    throw SequenceException ("Sequence::getAlphabetFromSequence : Both 'T' and 'U' in the same sequence!");
  if (r && p)
    throw SequenceException ("Sequence::getAlphabetFromSequence : Protein character and 'U' in the same sequence!");
  if (p)
    return new ProteicAlphabet();
  if (r)
    return new RNA();
  return new DNA();
}

/****************************************************************************************/

