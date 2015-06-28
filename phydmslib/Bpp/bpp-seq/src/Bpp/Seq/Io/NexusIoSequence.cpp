//
// File: NexusIOSequence.cpp
// Created by: Julien Dutheil
// Created on: Wed May 27 16:15 2009
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

#include "NexusIoSequence.h"
#include "NexusTools.h"
#include "../Container/SiteContainerTools.h"
#include "../Alphabet/AlphabetTools.h"
#include <Bpp/Text/TextTools.h>
#include <Bpp/Text/KeyvalTools.h>
#include <Bpp/Io/FileTools.h>

using namespace bpp;

// From the STL:
#include <sstream>

using namespace std;

/******************************************************************************/

const std::vector<std::string> NexusIOSequence::splitNameAndSequence_(const std::string& s) const throw (Exception)
{
  vector<string> v(2);
  string::size_type index = s.find(" ");
  if(index == string::npos) throw Exception("NexusIOSequence::splitNameAndSequence_(). No sequence name found.");
  v[0] = TextTools::removeSurroundingWhiteSpaces(s.substr(0, index));
  v[1] = TextTools::removeFirstWhiteSpaces(s.substr(index + 1));
  return v;
}  

  
/******************************************************************************/

void NexusIOSequence::appendAlignmentFromStream(std::istream& input, SiteContainer& vsc) const throw (Exception)
{
  // Checking the existence of specified file
  if (!input) { throw IOException ("NexusIOSequence::read(). Fail to open file"); }

  //Look for the DATA block:
  string line = "";
  while (TextTools::toUpper(line) != "BEGIN DATA;")
  {
    if (input.eof())
      throw Exception("NexusIOSequence::appendFromStream(). No data block was found.");
    line = TextTools::removeSurroundingWhiteSpaces(FileTools::getNextLine(input));
  }

  //Look for the DIMENSIONS command:
  string cmdName = "", cmdArgs = "";
  while (cmdName != "DIMENSIONS")
  {
    if (input.eof())
      throw Exception("NexusIOSequence::appendFromStream(). No DIMENSIONS command was found.");
    NexusTools::getNextCommand(input, cmdName, cmdArgs);
    cmdName = TextTools::toUpper(cmdName);
  }
  map<string, string> args;
  KeyvalTools::multipleKeyvals(cmdArgs, args, " ");
  map<string, string> argsUp;
  for (map<string, string>::iterator it = args.begin(); it != args.end(); it++)
    argsUp[TextTools::toUpper(it->first)] = it->second;
  if (argsUp["NTAX"] == "")
    throw Exception("NexusIOSequence::appendFromStream(). DIMENSIONS command does not have a NTAX argument.");
  unsigned int ntax = TextTools::to<unsigned int>(argsUp["NTAX"]);

  //Look for the FORMAT command:
  while (cmdName != "FORMAT")
  {
    if (input.eof())
      throw Exception("NexusIOSequence::appendFromStream(). No FORMAT command was found.");
    NexusTools::getNextCommand(input, cmdName, cmdArgs);
    cmdName = TextTools::toUpper(cmdName);
  }
  if (TextTools::hasSubstring(cmdArgs, "TRANSPOSE"))
    throw Exception("NexusIOSequence::appendFromStream(). TRANSPOSE option is not supported.");

  //Check if the alignment is dotted or not:
  bool matchChar = TextTools::hasSubstring(TextTools::toUpper(cmdArgs), "MATCHCHAR");

  SiteContainer* alignment = 0;
  if (matchChar)
    alignment = new AlignedSequenceContainer(&AlphabetTools::DEFAULT_ALPHABET);
  else
    alignment = &vsc;

  //Look for the MATRIX command:
  line = "";
  while (!TextTools::startsWith(TextTools::toUpper(line), "MATRIX"))
  {
    if (input.eof())
      throw Exception("NexusIOSequence::appendFromStream(). No MATRIX command was found.");
    line = TextTools::removeSurroundingWhiteSpaces(FileTools::getNextLine(input));
  }
  line = FileTools::getNextLine(input);

  vector<string> names, seqs;
  // Read first block:
  bool commandFinished = false;
  for (unsigned int i = 0; i < ntax && !input.eof(); i++)
  {
    if (TextTools::endsWith(line, ";"))
    {
      if (i < ntax - 1)
        throw IOException("NexusIOSequence::appendFromStream. Early end of MATRIX command, some sequences are missing.");
      else 
      {
        commandFinished = true;
        line = line.substr(0, line.size() - 1); //Remove trailing semi-colon.
      }
    }
    vector<string> v = splitNameAndSequence_(line);
    names.push_back(v[0]);
    seqs.push_back(v[1]);
    line = FileTools::getNextLine(input);
  }
  
  //Then read all other blocks:
  commandFinished = TextTools::removeSurroundingWhiteSpaces(line) == ";"; //In case the end of command is on a separate line.
  while (!commandFinished)
  {
    for (unsigned int i = 0; i < ntax && !input.eof(); i++)
    {
      if (TextTools::endsWith(line, ";"))
      {
        if (i < ntax - 1)
          throw IOException("NexusIOSequence::appendFromStream. Early end of MATRIX command, some sequences are missing.");
        else 
        {
          commandFinished = true;
          line = line.substr(0, line.size() - 1); //Remove trailing semi-colon.
        }
      }

      vector<string> v = splitNameAndSequence_(line);
      if (v[0] != names[i])
        throw IOException("NexusIOSequence::appendFromStream. Bad file, the sequences are not in the same order in interleaved blocks, or one taxon is missing.");
      seqs[i] += v[1];      
      line = FileTools::getNextLine(input);
      commandFinished = TextTools::removeSurroundingWhiteSpaces(line) == ";"; //In case the end of command is on a separate line.
    }
  }
  for (unsigned int i = 0; i < names.size(); i++)
  {
    alignment->addSequence(BasicSequence(names[i], seqs[i], vsc.getAlphabet()), checkNames_);
  }

  if (matchChar)
  {
    //Now we resolve the alignment:
    SiteContainer* resolvedAlignment =
      SiteContainerTools::resolveDottedAlignment(*alignment, vsc.getAlphabet());
    delete alignment;
    for (unsigned int i = 0; i < resolvedAlignment->getNumberOfSequences(); i++)
    {
      vsc.addSequence(resolvedAlignment->getSequence(i), false);
    }
    delete resolvedAlignment;
  }
}

/******************************************************************************/

const std::string NexusIOSequence::getFormatName() const { return "Nexus"; }

/******************************************************************************/

const std::string NexusIOSequence::getFormatDescription() const
{
  return "Nexus file format.";
}

/******************************************************************************/

