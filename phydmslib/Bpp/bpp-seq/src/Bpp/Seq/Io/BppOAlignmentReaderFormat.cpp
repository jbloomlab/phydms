//
// File: BppOAlignmentReaderFormat.cpp
// Created by: Julien Dutheil
// Created on: Friday September 15th, 22:06
//

/*
  Copyright or © or Copr. Bio++ Development Team, (November 16, 2004)

  This software is a computer program whose purpose is to provide classes
  for phylogenetic data analysis.

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

#include "BppOAlignmentReaderFormat.h"
#include "Mase.h"
#include "Phylip.h"
#include "Fasta.h"
#include "Clustal.h"
#include "Dcse.h"
#include "NexusIoSequence.h"

#include <Bpp/Text/KeyvalTools.h>

#include <string>
#include <memory>

using namespace bpp;
using namespace std;

IAlignment* BppOAlignmentReaderFormat::read(const std::string& description) throw (Exception)
{
  unparsedArguments_.clear();
  string format = "";
  KeyvalTools::parseProcedure(description, format, unparsedArguments_);
  auto_ptr<IAlignment> iAln;
  if (format == "Mase")
  {
    iAln.reset(new Mase());
  }
  else if (format == "Phylip")
  {
    bool sequential = true, extended = true;
    string split = "  ";
    string order = ApplicationTools::getStringParameter("order", unparsedArguments_, "sequential", "", true, warningLevel_);
    if (order == "sequential")
      sequential = true;
    else if (order == "interleaved")
      sequential = false;
    else
      throw Exception("BppOAlignmentReaderFormat::read. Invalid argument 'order' for phylip format: " + order);
    
    string type = ApplicationTools::getStringParameter("type", unparsedArguments_, "extended", "", true, warningLevel_);
    if (type == "extended")
    {
      extended = true;
      split = ApplicationTools::getStringParameter("split", unparsedArguments_, "spaces", "", true, warningLevel_);
      if (split == "spaces")
        split = "  ";
      else if (split == "tab")
        split = "\t";
      else
        throw Exception("BppOAlignmentReaderFormat::read. Invalid argument 'split' for phylip format: " + split);
    }
    else if (type == "classic")
      extended = false;
    else
      throw Exception("BppOAlignmentReaderFormat::read. Invalid argument 'type' for phylip format: " + type);
    
    iAln.reset(new Phylip(extended, sequential, 100, true, split));
  }
  else if (format == "Fasta")
  {
    bool strictNames = ApplicationTools::getBooleanParameter("strict_names", unparsedArguments_, false, "", true, warningLevel_);
    bool extended    = ApplicationTools::getBooleanParameter("extended", unparsedArguments_, false, "", true, warningLevel_);
    iAln.reset(new Fasta(100, true, extended, strictNames));
  }
  else if (format == "Clustal")
  {
    unsigned int extraSpaces = ApplicationTools::getParameter<unsigned int>("extraSpaces", unparsedArguments_, 0, "", true, warningLevel_);
    iAln.reset(new Clustal(true, extraSpaces));
  }
  else if (format == "Dcse")
  {
    iAln.reset(new DCSE());
  }
  else if (format == "Nexus")
  {
    iAln.reset(new NexusIOSequence());
  }
  else
  {
    throw Exception("Sequence format '" + format + "' unknown.");
  }

  return iAln.release();
}

