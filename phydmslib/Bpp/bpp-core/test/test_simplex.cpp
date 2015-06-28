//
// File: test_simplex.cpp
// Created by: Laurent Guéguen
// Created on: vendredi 5 juillet 2013, à 11h 08
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

#include <Bpp/Numeric/Prob/Simplex.h>
#include <Bpp/Text/TextTools.h>

#include <iostream>
#include <limits>

using namespace bpp;
using namespace std;

int main()
{
  try
  {
    vector<double> prob;
    prob.push_back(0.1);
    prob.push_back(0.2);
    prob.push_back(0.3);
    prob.push_back(0.15);
    prob.push_back(0.1);
    prob.push_back(0.05);
    prob.push_back(0.1);

    vector<Simplex*> vpsi;
    vpsi.push_back(new Simplex(prob, 1));
    vpsi.push_back(new Simplex(prob, 2));
    vpsi.push_back(new Simplex(prob, 3));

    for (size_t i = 0; i < 3; i++)
    {
      cout << "Method " << i + 1 << endl;
      for (size_t j = 0; j < prob.size() - 1; j++)
      {
        cout << vpsi[i]->getParameterValue("theta" + TextTools::toString(j + 1)) << "\t";
      }
      cout << endl;
    }


    cout << "Prob:";
    for (size_t j = 0; j < prob.size(); j++)
    {
      cout << prob[j] << "\t";
    }
    cout << endl;
    for (size_t i = 0; i < 3; i++)
    {
      for (size_t j = 0; j < prob.size() - 1; j++)
      {
        vpsi[i]->setParameterValue("theta" + TextTools::toString(j + 1),
                                   vpsi[i]->getParameterValue("theta" + TextTools::toString(j + 1)) + 0.1);
      }
      for (size_t j = 0; j < prob.size() - 1; j++)
      {
        vpsi[i]->setParameterValue("theta" + TextTools::toString(j + 1),
                                   vpsi[i]->getParameterValue("theta" + TextTools::toString(j + 1)) - 0.1);
      }
      cout << "Method " << i + 1 << endl;
      cout << "prob\t";
      for (size_t j = 0; j < prob.size(); j++)
      {
        cout << vpsi[i]->prob(j) << "\t";
      }
      cout << endl;
    }

    for (size_t i = 0; i < 3; i++)
    {
      delete vpsi[i];
    }
    return 0;
  }
  catch (Exception& ex)
  {
    cout << "failed :(" << endl;
    return 1;
  }
}

