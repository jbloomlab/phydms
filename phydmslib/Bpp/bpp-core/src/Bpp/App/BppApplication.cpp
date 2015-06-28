//
// File: BppApplication.cpp
// Created by: Julien Dutheil
// Created on: Sat Aug 08 08:21 2009
//

/*
Copyright or © or Copr. Bio++ Development Team, (November 17, 2004)

This software is a computer program whose purpose is to provide basal and 
utilitary classes. This file belongs to the Bio++ Project.

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

#include "BppApplication.h"
#include "../Utils/AttributesTools.h"
#include "../Numeric/Random/RandomTools.h"
#include "ApplicationTools.h"

// From the STL:
#include <iostream>

using namespace bpp;
using namespace std;

BppApplication::BppApplication(int argc, char* argv[], const std::string& name): appName_(name), params_(), timerStarted_(false)
{
  cout << "Parsing options:" << endl;  
  params_ = AttributesTools::parseOptions(argc, argv);
  ApplicationTools::warningLevel = ApplicationTools::getIntParameter("--warning", params_, 0, "", true, 3);
  bool noint = ApplicationTools::getBooleanParameter("--noninteractive", params_, false, "", true, 3);
  ApplicationTools::interactive = !noint;
  long seed = ApplicationTools::getParameter<long>("--seed", params_, -1, "", true, 3);
  if (seed >= 0) {
    RandomTools::setSeed(seed);
    ApplicationTools::displayResult("Random seed set to", seed);
  }
}

void BppApplication::startTimer()
{
  ApplicationTools::startTimer();
  timerStarted_ = true;
}

void BppApplication::done()
{
  cout << appName_ << "'s done. Bye." << endl;
  if (timerStarted_)
    ApplicationTools::displayTime("Total execution time:");
}

