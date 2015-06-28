//
// File: MaseTools.cpp
// Created by: Julien Dutheil
// Created on: Tue Apr  1 09:16:59 2003
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

#include "MaseTools.h"
#include "../Container/VectorSequenceContainer.h"
#include "../Container/AlignedSequenceContainer.h"
#include "../Container/SequenceContainerTools.h"
#include <Bpp/Text/StringTokenizer.h>
#include <Bpp/Text/TextTools.h>

#include <iostream>

using namespace std;
using namespace bpp;

SiteSelection MaseTools::getSiteSet(const Comments& maseFileHeader, const string& setName) throw (IOException)
{
  SiteSelection selection;
  for (size_t i = 0; i < maseFileHeader.size(); i++)
  {
    string current = maseFileHeader[i];
    string::size_type index = current.find("# of");
    if (index < current.npos)
    {
      StringTokenizer st(string(current.begin() + static_cast<ptrdiff_t>(index + 4), current.end()), " \t=;");
      st.nextToken(); // skip next word: may be 'regions' or 'segments' or else ;-)
      size_t numberOfSegments = TextTools::to<size_t>(st.nextToken());
      string name = st.unparseRemainingTokens();
      if (name == setName)
      {
        // cout << numberOfSegments << " segments found." << endl;
        // Then look for the set definition:
        i++; // next line.
        size_t counter = 0;
        while (i < maseFileHeader.size())
        {
          current = maseFileHeader[i++];
          StringTokenizer st2(current);
          // st.nextToken(); //Skip ';;'
          while (st2.hasMoreToken())
          {
            StringTokenizer st3(st2.nextToken(), ",");
            size_t begin = TextTools::to<size_t>(st3.nextToken());
            size_t end   = TextTools::to<size_t>(st3.nextToken());
            // WARNING!!! In the mase+ format, sites are numbered from 1 to nbSites,
            // Whereas in SiteContainer the index begins at 0.
            for (size_t j = begin; j <= end; j++)
            {
              selection.push_back(j - 1); // bounds included.
            }
            counter++;
            if (counter == numberOfSegments)
              return selection;
          }
        }
      }
    }
  }
  if (selection.size() == 0)
  {
    throw IOException("Site set " + setName + " has not been found in the sequence file.");
  }
  return selection;
}

/******************************************************************************/

SequenceSelection MaseTools::getSequenceSet(const Comments& maseFileHeader, const string& setName) throw (IOException)
{
  SequenceSelection selection;
  for (size_t i = 0; i < maseFileHeader.size(); i++)
  {
    string current = maseFileHeader[i];

    string::size_type index = current.find("@ of");
    if (index < current.npos)
    {
      StringTokenizer st(string(current.begin() + static_cast<ptrdiff_t>(index + 4), current.end()), " \t=;");
      st.nextToken(); // skip next word: may be 'sequences' or else ;-)
      size_t numberOfSequences = TextTools::to<size_t>(st.nextToken());
      string name = st.unparseRemainingTokens();
      size_t counter = 0;
      if (name == setName)
      {
        // cout << numberOfSequences << " segments found." << endl;
        // Then look for the set definition:
        i++; // next line.
        while (i < maseFileHeader.size())
        {
          current = maseFileHeader[i++];
          StringTokenizer st2(current, ",");
          while (st2.hasMoreToken())
          {
            int seqIndex = TextTools::toInt(st2.nextToken());
            // WARNING!!! In the mase+ format, sequences are numbered from 1 to nbSequences,
            // Whereas in SequenceContainer the index begins at 0.
            selection.push_back(static_cast<size_t>(seqIndex - 1)); // bounds included.
            counter++;
            if (counter == numberOfSequences)
              return selection;
          }
        }
      }
    }
  }
  if (selection.size() == 0)
  {
    throw IOException("Sequence set " + setName + " has not been found in the sequence file.");
  }
  return selection;
}

/******************************************************************************/

SiteContainer* MaseTools::getSelectedSites(
  const SiteContainer& sequences,
  const string& setName) throw (IOException)
{
  SiteSelection ss = getSiteSet(sequences.getGeneralComments(), setName);
  return SiteContainerTools::getSelectedPositions(sequences, ss);
}

/******************************************************************************/

SequenceContainer* MaseTools::getSelectedSequences(
  const OrderedSequenceContainer& sequences,
  const std::string& setName) throw (IOException)
{
  SequenceSelection ss = getSequenceSet(sequences.getGeneralComments(), setName);
  VectorSequenceContainer* cont = new VectorSequenceContainer(sequences.getAlphabet());
  SequenceContainerTools::getSelectedSequences(sequences, ss, *cont);
  return cont;
}

/******************************************************************************/

map<string, size_t> MaseTools::getAvailableSiteSelections(const Comments& maseHeader)
{
  map<string, size_t> selections;
  for (size_t i = 0; i < maseHeader.size(); i++)
  {
    string current = maseHeader[i];

    string::size_type index = current.find("# of");
    if (index < current.npos)
    {
      StringTokenizer st(string(current.begin() + static_cast<ptrdiff_t>(index + 4), current.end()), " \t\n\f\r=;");
      st.nextToken(); // skip next word: may be 'sequences' or else ;-)
      size_t numberOfSegments = TextTools::to<size_t>(st.nextToken());
      string name = st.nextToken();
      while (st.hasMoreToken())
      {
        name += " " + st.nextToken();
      }
      size_t counter = 0;
      size_t nbSites = 0;
      while (i < maseHeader.size())
      {
        i++;
        current = maseHeader[i];
        StringTokenizer st2(current);
        // st.nextToken(); //Skip ';;'
        while (st2.hasMoreToken())
        {
          StringTokenizer st3(st2.nextToken(), ",");
          size_t begin = TextTools::to<size_t>(st3.nextToken());
          size_t end   = TextTools::to<size_t>(st3.nextToken());
          counter++;
          nbSites += end - begin + 1;
        }
        if (counter == numberOfSegments)
        {
          selections[name] = nbSites;
          break;
        }
      }
    }
  }
  return selections;
}

/******************************************************************************/

map<string, size_t> MaseTools::getAvailableSequenceSelections(const Comments& maseHeader)
{
  map<string, size_t> selections;
  for (size_t i = 0; i < maseHeader.size(); i++)
  {
    string current = maseHeader[i];

    string::size_type index = current.find("@ of");
    if (index < current.npos)
    {
      StringTokenizer st(string(current.begin() + static_cast<ptrdiff_t>(index + 4), current.end()), " \t\n\f\r=;");
      st.nextToken(); // skip next word: may be 'sequences' or else ;-)
      size_t numberOfSequences = TextTools::fromString<size_t>(st.nextToken());
      string name = st.nextToken();
      while (st.hasMoreToken())
      {
        name += st.nextToken();
      }
      selections[name] = numberOfSequences;
    }
  }
  return selections;
}

/******************************************************************************/

size_t MaseTools::getPhase(const Comments& maseFileHeader, const string& setName) throw (Exception)
{
  size_t phase = 0;
  string::size_type index = 0;
  for (size_t i = 0; i < maseFileHeader.size(); i++)
  {
    string current = maseFileHeader[i];

    index = current.find("# of");
    if (index < current.npos)
    {
      StringTokenizer st(string(current.begin() + static_cast<ptrdiff_t>(index + 12), current.end()), " \t\n\f\r=;");
      // size_t numberOfSegments = TextTools::toInt(st.nextToken());
      // cout << "Number of regions: " << st.nextToken() << endl;
      string name;
      while (st.hasMoreToken())
      {
        name = st.nextToken();
        // cout << "Name of regions: " << name << endl;
      }
      if (name == setName)
      {
        return phase;
      }
    }

    index = current.find("/codon_start");
    if (index < current.npos)
    {
      StringTokenizer st(string(current.begin() + static_cast<ptrdiff_t>(index + 12), current.end()), " \t\n\f\r=;");
      phase = TextTools::to<size_t>(st.nextToken());
    }
  }
  throw Exception("PolymorphismSequenceContainer::getPhase: no /codon_start found, or site selection missing.");
}

/******************************************************************************/

