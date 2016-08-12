//
// File: test_bfgs.cpp
// Created by: Julien Dutheil
// Created on: Mon Jan 31 19:45 2011
//

/*
Copyright or © or Copr. Bio++ Development Team, (November 17, 2004)

This software is a computer program whose purpose is to provide classes
for numerical calculus. This file is part of the Bio++ project.

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

#include "../src/Bpp/Graph/TreeGraph.h"
#include "../src/Bpp/Graph/AssociationTreeGraphObserver.h"


#include <vector>
#include <iostream>
using namespace bpp;
using namespace std;

int main() {
  SimpleAssociationTreeGraphObserver<string,unsigned int,bpp::SimpleTreeGraph<SimpleGraph> > grObs(true);
  string zero = "zero";
  string one = "one";
  string two = "two";
  string three = "three";
  unsigned int r1 = 1;
  unsigned int r3 = 3;
  cout << "Creating node zero." << endl;
  grObs.createNode(&zero);
  grObs.getGraph()->outputToDot(std::cout,"myTestDirGrObs");
  cout << "Creating node one from the number zero." << endl;
  grObs.createNode(&one);
  grObs.link(&zero,&one,&r1);
  grObs.getGraph()->outputToDot(std::cout,"myTestDirGrObs");
  cout << "The father branch of zero is " << *(grObs.getBranchToFather(&one)) << endl;
  cout << "The father node of one is " << *(grObs.getFather(&one)) << endl;
  
//   cout << "unlink 0->1 and relink with setFather()";
//   grObs.unlink(&zero,&one);
//   grObs.setFather(&one,&zero);
//FIXME: setFather only works with valid tree (wich cannot be true…)
  grObs.getGraph()->outputToDot(std::cout,"myTestDirGrObs");
  
  
  cout << "Creating node two from the number one." << endl;
  grObs.createNode(&one,&two);
  grObs.getGraph()->outputToDot(std::cout,"myTestDirGrObs");
  cout << "Linking two to zero." << endl;
  grObs.link(&two,&zero,&r3);
  grObs.getGraph()->outputToDot(std::cout,"myTestDirGrObs");  
  cout << "Linking one to three." << endl;
  grObs.createNode(&one,&three);
  grObs.getGraph()->outputToDot(std::cout,"myTestDirGrObs");
  

  cout << "Linking three to zero." << endl;
  grObs.link(&three,&zero);
  grObs.getGraph()->outputToDot(std::cout,"myTestDirGrObs");
  // so now we have zero -> one -> two -> zero ...
  vector<string*> fromOne = grObs.getOutgoingNeighbors(&zero);
  vector<string*> fromThree = grObs.getOutgoingNeighbors(&two);
  bool test = (*(fromOne.begin()) == &one) && (*(fromThree.begin()) == &zero);
  grObs.getGraph()->outputToDot(std::cout,"myTestDirGrObs");
  
  
  
  cout << "Is this a tree?\n    " << (grObs.isValid()? "TRUE":"FALSE") << endl;
  
  // the tree must be considered as unvalid at this point
  test &= !grObs.isValid();
  
  
  cout << "Removing 3->0"<< endl;
  grObs.unlink(&two,&zero);
  grObs.unlink(&three,&zero);
  grObs.getGraph()->outputToDot(std::cout,"myTestDirGrObs");
  
  cout << "Is this a tree?\n    " << (grObs.isValid()? "TRUE":"FALSE") << endl;
  // the tree is valid at this point
  test &= grObs.isValid();
  
  cout << "Re-linking 3->0"<< endl;
  grObs.link(&two,&zero);
  grObs.getGraph()->outputToDot(std::cout,"myTestDirGrObs");
  cout << "Is this a tree?\n    " << (grObs.isValid()? "TRUE":"FALSE") << endl;
  // the tree must be considered as unvalid at this point
  test &= !grObs.isValid();
  
  cout << "Test " << (test? "passed":"failed") << endl;
  
  return (test ? 0 : 1);
}
