//
// File: Phylip.cpp
// Created by: Julien Dutheil
// Created on: Mon Oct 27 12:22:56 2003
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

#include "Phylip.h"
#include "../Container/SequenceContainerTools.h"
#include <Bpp/Text/TextTools.h>
#include <Bpp/Text/StringTokenizer.h>
#include <Bpp/Io/FileTools.h>

using namespace bpp;

// From the STL:
#include <sstream>

using namespace std;

/******************************************************************************/

const std::vector<std::string> Phylip::splitNameAndSequence(const std::string& s) const throw (Exception)
{
  vector<string> v(2);
  if (extended_)
  {
    string::size_type index = s.find(namesSplit_);
    if(index == string::npos) throw Exception("No sequence name found.");
    v[0] = TextTools::removeSurroundingWhiteSpaces(s.substr(0, index));
    v[1] = TextTools::removeFirstWhiteSpaces      (s.substr(index + namesSplit_.size())); //There may be more than 2 white spaces.
  }
  else
  {
    v[0] = TextTools::removeSurroundingWhiteSpaces(s.substr(0, 10));
    v[1] = s.substr(10);
  }
  return v;
}  

/******************************************************************************/

void Phylip::readSequential(std::istream& in, SiteContainer& asc) const throw (Exception)
{
  string temp;
  
  //Ignore first line:
  getline(in, temp, '\n');  // Copy current line in temporary string
  temp = TextTools::removeSurroundingWhiteSpaces(FileTools::getNextLine(in));
  string name = "";
  string seq  = "";
  
  while (!in.eof())
  {
    // Read each sequence:
    vector<string> v;
    bool hasName = true;
    try
    { 
      v = splitNameAndSequence(temp);
    }
    catch (Exception & e)
    {
      hasName = false;
    }
    if (hasName)
    {
      // a new sequence is found:
      if (!TextTools::isEmpty(name)) //If this is not the first sequence!
      {
        // Add the previous sequence to the container:
        asc.addSequence(BasicSequence(name, seq, asc.getAlphabet()), checkNames_);
      }
      name = v[0];
      seq  = v[1];
    }
    else
    {
      //No sequence name found.
      if (TextTools::isEmpty(name))
        throw Exception("First sequence in file has no name!");
      seq += TextTools::removeWhiteSpaces(temp);
    }
    //while(!TextTools::isEmpty(temp))
    //{
    //  //Sequences are separated by at least one blank line:
    //  getline(in, temp, '\n');  // read next line in file.
    //  seq += TextTools::removeWhiteSpaces(temp);      
    //}
    //end of this sequence:
    temp = TextTools::removeSurroundingWhiteSpaces(FileTools::getNextLine(in));

  }
  // Add last sequence:
  asc.addSequence(BasicSequence(name, seq, asc.getAlphabet()), checkNames_);
}

/******************************************************************************/

void Phylip::readInterleaved(std::istream& in, SiteContainer& asc) const throw (Exception)
{
  string temp;
  
  //Read first line:
  getline(in, temp, '\n'); // Copy current line in temporary string
  StringTokenizer st(temp);
  unsigned int nbSequences = TextTools::to<unsigned int>(st.nextToken());
  //int nbSites     = TextTools::toInt(st.nextToken());
  temp = FileTools::getNextLine(in);
  
  vector<string> names, seqs;
  // Read first block:
  for (unsigned int i = 0; i < nbSequences && !in.eof() && !TextTools::isEmpty(temp); i++)
  {
    vector<string> v = splitNameAndSequence(temp);
    names.push_back(v[0]);
    seqs.push_back(v[1]);
    getline(in, temp, '\n');  // read next line in file.
  }
  
  //Then read all other blocks:
  temp = FileTools::getNextLine(in);
  while (!in.eof())
  {
    for (unsigned int i = 0; i < names.size(); i++)
    {
      if (TextTools::isEmpty(temp))
        throw IOException("Phylip::readInterleaved. Bad file,there are not the same number of sequence in each block.");
      seqs[i] += TextTools::removeWhiteSpaces(temp);      
      getline(in, temp, '\n');  // read next line in file.
    }
    temp = FileTools::getNextLine(in);
  }
  for (unsigned int i = 0; i < names.size(); i++)
  {
    asc.addSequence(BasicSequence(names[i], seqs[i], asc.getAlphabet()), checkNames_);
  }
}
  
/******************************************************************************/

void Phylip::appendAlignmentFromStream(std::istream& input, SiteContainer& vsc) const throw (Exception)
{
  // Checking the existence of specified file
  if (!input) { throw IOException ("Phylip::read: fail to open file"); }
  
  if(sequential_) readSequential (input, vsc);
  else            readInterleaved(input, vsc);
}

/******************************************************************************/

unsigned int Phylip::getNumberOfSequences(const std::string& path) const throw (IOException)
{
  // Checking the existence of specified file
  ifstream file (path.c_str(), ios::in);
  if (! file) { throw IOException ("Phylip::getNumberOfSequences: failed to open file"); }
  string firstLine = FileTools::getNextLine(file);
  StringTokenizer st(firstLine, " \t");
  istringstream iss(st.nextToken());
  unsigned int nb;
  iss >> nb;
  file.close();
  return nb;
}
 
/******************************************************************************/

std::vector<std::string> Phylip::getSizedNames(const std::vector<std::string>& names) const
{
  vector<string> sizedNames(names.size());
  if (extended_)
  {
    //Add 6 white spaces to the larger name and align other names.
    //First, determine the size of the wider name:
    size_t sizeMax = 0;
    for (size_t i = 0; i < names.size(); i++)
      if (names[i].size() > sizeMax) sizeMax = names[i].size();
    //Quite easy ;-) Now update all lengths:
    for (size_t i = 0; i < names.size(); i++)
      sizedNames[i] = TextTools::resizeRight(names[i], sizeMax) + namesSplit_;  
  }
  else
  {
    //We trunc all names to ten characters:
    for(unsigned int i = 0; i < names.size(); i++) sizedNames[i] = TextTools::resizeRight(names[i], 10);
    cout << "Warning: names have been truncated to 10 characters. They may be ambiguous sequence names then." << endl;
  }
  return sizedNames;
}

/******************************************************************************/

void Phylip::writeSequential(std::ostream& out, const SequenceContainer& sc) const
{
  //cout << "Write sequential" << endl;
  size_t numberOfSites = sc.getSequence(sc.getSequencesNames()[0]).size() * sc.getAlphabet()->getStateCodingSize();
  out << sc.getNumberOfSequences() << " " << numberOfSites << endl;
  
  vector<string> seqNames = sc.getSequencesNames();
  vector<string> names = getSizedNames(seqNames);
  for (size_t i = 0; i < seqNames.size(); ++i)
  {
    vector<string> seq = TextTools::split(sc.toString(seqNames[i]), charsByLine_);
    out << names[i] << seq[0] << endl;
    for (size_t j = 1; j < seq.size(); ++j)
    {
      out << string(names[i].size(), ' ') << seq[j] << endl;
    }
    out << endl;
  }
}

void Phylip::writeInterleaved(std::ostream& out, const SequenceContainer& sc) const
{
  //cout << "Write interleaved;" << endl;
  size_t numberOfSites = sc.getSequence(sc.getSequencesNames()[0]).size() * sc.getAlphabet()->getStateCodingSize();
  out << sc.getNumberOfSequences() << " " << numberOfSites << endl;
  
  vector<string> seqNames = sc.getSequencesNames();
  vector<string> names = getSizedNames(seqNames);
  //Split sequences:
  vector< vector<string> > seqs(sc.getNumberOfSequences());
  for (size_t i = 0; i < seqNames.size(); ++i)
  {
    seqs[i] = TextTools::split(sc.toString(seqNames[i]), charsByLine_);
  }
  //Write first block:
  for (size_t i = 0; i < names.size(); ++i)
  {
    out << names[i] << seqs[i][0] << endl;
  }
  out << endl;
  //Write other blocks:
  for (size_t j = 1; j < seqs[0].size(); ++j)
  {
    for (size_t i = 0; i < sc.getNumberOfSequences(); ++i)
    {
      out << seqs[i][j] << endl;
    }
    out << endl;
  }
}

/******************************************************************************/

void Phylip::writeAlignment(std::ostream& output, const SiteContainer& sc) const throw (Exception)
{
  //First must check if all sequences are aligned:
  if (sc.getNumberOfSequences() == 0)
    throw Exception("Phylip::write. SequenceContainer appear to contain no sequence.");
  
  // Checking the existence of specified file, and possibility to open it in write mode
  if (!output) { throw IOException ("Phylip::write : failed to open file"); }

  if (sequential_) writeSequential (output, sc);
  else             writeInterleaved(output, sc);
}

/******************************************************************************/

const std::string Phylip::getFormatName() const { return "Phylip file, " + string(extended_ ? "extended," : "") + string(sequential_ ? "sequential" : "interleaved"); }

/******************************************************************************/

const std::string Phylip::getFormatDescription() const
{
  return "Phylip file format, sequential and interleaved. PAML extension also supported.";
}

/******************************************************************************/

