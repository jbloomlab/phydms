//
// File: BppOAlphabetIndex2Format.cpp
// Created by: Julien Dutheil
// Created on: Thursday Februar 07th, 19:26
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

#include "BppOAlphabetIndex2Format.h"
#include "BppOAlphabetIndex1Format.h"
#include "../Alphabet/AlphabetTools.h"
#include "../AlphabetIndex/BLOSUM50.h"
#include "../AlphabetIndex/GranthamAAChemicalDistance.h"
#include "../AlphabetIndex/MiyataAAChemicalDistance.h"
#include "../AlphabetIndex/SimpleIndexDistance.h"
#include "../AlphabetIndex/AAIndex2Entry.h"
#include "../AlphabetIndex/AlphabetIndex1.h"

#include <Bpp/Text/KeyvalTools.h>
#include <Bpp/App/ApplicationTools.h>

#include <string>
#include <memory>

using namespace bpp;
using namespace std;

AlphabetIndex2* BppOAlphabetIndex2Format::read(const std::string& description) throw (Exception)
{
  if (description != "None")
  {
    string name;
    map<string, string> args;
    KeyvalTools::parseProcedure(description, name, args);
    if (verbose_)
      ApplicationTools::displayResult(message_, description);

    //Currently, only protein indices are supported:
    if (!AlphabetTools::isProteicAlphabet(alphabet_))
        throw Exception("BppOAlphabetIndex2Format::read. This index is only supported with a protein alphabet.");
    if (name == "Blosum50")
    {
      return new BLOSUM50();
    }
    else if (name == "Grantham")
    {
      bool sym = ApplicationTools::getBooleanParameter("symmetrical", args, true, "", true, 1);
      GranthamAAChemicalDistance* M = new GranthamAAChemicalDistance();
      M->setSymmetric(sym);
      if (!sym) M->setPC1Sign(true);
      return M;
    }
    else if (name == "Miyata")
    {
      bool sym = ApplicationTools::getBooleanParameter("symmetrical", args, true, "", true, 1);
      MiyataAAChemicalDistance* M = new MiyataAAChemicalDistance();
      M->setSymmetric(sym);
      return M;
    }
    else if (name == "Diff")
    {
      string index1Desc = ApplicationTools::getStringParameter("index1", args, "None", "", true, 1);
      bool sym = ApplicationTools::getBooleanParameter("symmetrical", args, true, "", true);
      BppOAlphabetIndex1Format index1Reader(alphabet_, "" , false);
      AlphabetIndex1* index1 = index1Reader.read(index1Desc);
      if (index1) {
        SimpleIndexDistance* M = new SimpleIndexDistance(index1);
        M->setSymmetric(sym);
        return M;
      } else {
        throw Exception("BppOAlphabetIndex2Format::read. Diff: index1 should be provided.");
      }
    }
    else if (name == "User")
    {
      bool sym = ApplicationTools::getBooleanParameter("symmetrical", args, true, "", true, 1);
      string aax2FilePath = ApplicationTools::getAFilePath("file", args, true, true, "", false);
      ifstream aax2File(aax2FilePath.c_str(), ios::in);
      AAIndex2Entry* M = new AAIndex2Entry(aax2File, sym);
      aax2File.close();
      return M;
    }
    else
    {
      throw Exception("Invalid index2 '" + name + "'.");
    }
  }
  else
  {
    return 0;
  }
}

