//
// File: BppOAlphabetIndex1Format.cpp
// Created by: Julien Dutheil
// Created on: Thursday Februar 07th, 16:30
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

#include "BppOAlphabetIndex1Format.h"
#include "../Alphabet/AlphabetTools.h"
#include "../AlphabetIndex/GranthamAAPolarityIndex.h"
#include "../AlphabetIndex/GranthamAAVolumeIndex.h"
#include "../AlphabetIndex/KleinAANetChargeIndex.h"
#include "../AlphabetIndex/AAChouFasmanAHelixIndex.h"
#include "../AlphabetIndex/AAChouFasmanBSheetIndex.h"
#include "../AlphabetIndex/AAChouFasmanTurnIndex.h"
#include "../AlphabetIndex/AAChenGuHuangHydrophobicityIndex.h"
#include "../AlphabetIndex/AASurfaceIndex.h"
#include "../AlphabetIndex/AAMassIndex.h"
#include "../AlphabetIndex/AAVolumeIndex.h"
#include "../AlphabetIndex/AAChargeIndex.h"
#include "../AlphabetIndex/AASEAInf10Index.h"
#include "../AlphabetIndex/AASEA1030Index.h"
#include "../AlphabetIndex/AASEASup30Index.h"
#include "../AlphabetIndex/AAIndex1Entry.h"

#include <Bpp/Text/KeyvalTools.h>
#include <Bpp/App/ApplicationTools.h>

#include <string>
#include <memory>

using namespace bpp;
using namespace std;

AlphabetIndex1* BppOAlphabetIndex1Format::read(const std::string& description) throw (Exception)
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
        throw Exception("BppOAlphabetIndex1Format::read. This index is only supported with a protein alphabet.");
    if (name == "GranthamPolarity")
    {
      return new GranthamAAPolarityIndex();
    }
    else if (name == "GranthamVolume")
    {
      return new GranthamAAVolumeIndex();
    }
    else if (name == "KleinCharge")
    {
      return new KleinAANetChargeIndex();
    }
    else if (name == "ChouFasmanAHelix")
    {
      return new AAChouFasmanAHelixIndex();
    }
    else if (name == "ChouFasmanBSheet")
    {
      return new AAChouFasmanBSheetIndex();
    }
    else if (name == "ChouFasmanTurn")
    {
      return new AAChouFasmanTurnIndex();
    }
    else if (name == "ChenGuHuangHydrophobicity")
    {
      return new AAChenGuHuangHydrophobicityIndex();
    }
    else if (name == "Surface")
    {
      return new AASurfaceIndex();
    }
    else if (name == "Mass")
    {
      return new AAMassIndex();
    }
    else if (name == "Volume")
    {
      return new AAVolumeIndex();
    }
    else if (name == "Charge")
    {
      return new AAChargeIndex();
    }
    else if (name == "SEAMedium")
    {
      return new AASEA1030Index();
    }
    else if (name == "SEAHigh")
    {
      return new AASEASup30Index();
    }
    else if (name == "SEALow")
    {
      return new AASEAInf10Index();
    }
    else if (name == "User")
    {
      string aax1FilePath = ApplicationTools::getAFilePath("file", args, true, true, "", false);
      ifstream aax1File(aax1FilePath.c_str(), ios::in);
      AAIndex1Entry* I = new AAIndex1Entry (aax1File);
      aax1File.close();
      return I;
    }
    else
    {
      throw Exception("Invalid index1 '" + name + "'.");
    }
  }
  else
  {
    return 0;
  }
}

