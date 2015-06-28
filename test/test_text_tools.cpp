//
// File: test_text_tools.cpp
// Created by: Julien Dutheil
// Created on: Mon Nov 5 16:12 2011
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

#include <Bpp/Text/TextTools.h>
#include <Bpp/Text/StringTokenizer.h>
#include <iostream>

using namespace bpp;
using namespace std;

int main()
{
  cout << "Testing string conversion..." << endl;
  if ( TextTools::isDecimalNumber("aazz")) { cout << "aazz is not a decimal number!" << endl; return 1; }
  if ( TextTools::isDecimalNumber("-aazz")) { cout << "-aazz is not a decimal number!" << endl; return 1; }
  if ( TextTools::isDecimalNumber("-3.45z")) { cout << "-3.45z is not a decimal number!" << endl; return 1; }
  if (!TextTools::isDecimalNumber("0")) { cout << "0 is a decimal number!" << endl; return 1; }
  if (!TextTools::isDecimalNumber("123")) { cout << "123 is a decimal number!" << endl; return 1; }
  if (!TextTools::isDecimalNumber("-123")) { cout << "-123 is a decimal number!" << endl; return 1; }
  if (!TextTools::isDecimalNumber("-123.456")) { cout << "-123.456 is a decimal number!" << endl; return 1; }
  if (!TextTools::isDecimalInteger("123456")) { cout << "123456 is a decimal integer!" << endl; return 1; }
  if (!TextTools::isDecimalInteger("-7890")) { cout << "-7890 is a decimal integer!" << endl; return 1; }
  if (!TextTools::isDecimalNumber("-123.456e-5")) { cout << "-123.456e-5 is a decimal number!" << endl; return 1; }
  if ( TextTools::isDecimalNumber("-123.456e-5.8")) { cout << "-123.456e-5.8 is not a decimal number!" << endl; return 1; }
  if (!TextTools::isDecimalInteger("-123e6")) { cout << "-123e6 is a decimal integer!" << endl; return 1; }
  if ( TextTools::isDecimalInteger("-123.456e5")) { cout << "-123.456e5 is not a decimal integer!" << endl; return 1; }
  if ( TextTools::isDecimalInteger("-123e-6")) { cout << "-123e-6 is not a decimal integer!" << endl; return 1; }

  cout << "Testing string tokenizer..." << endl;
  string t;
  StringTokenizer st1(" aaazzer  aeerd a    eer", " \t", false, false);
  if (st1.numberOfRemainingTokens() != 4) return 1;
  cout << (t = st1.nextToken()) << endl;
  if (t != "aaazzer") return 1;
  cout << (t = st1.nextToken()) << endl;
  if (t != "aeerd") return 1;
  cout << (t = st1.nextToken()) << endl;
  if (t != "a") return 1;
  cout << (t = st1.nextToken()) << endl;
  if (t != "eer") return 1;

  StringTokenizer st2(" aaazzer  aeerd a    eer", " \t", false, true);
  if (st2.numberOfRemainingTokens() != 8) return 1;
  cout << (t = st2.nextToken()) << endl;
  if (t != "aaazzer") return 1;
  cout << (t = st2.nextToken()) << endl;
  if (t != "") return 1;
  cout << (t = st2.nextToken()) << endl;
  if (t != "aeerd") return 1;
  cout << (t = st2.nextToken()) << endl;
  if (t != "a") return 1;
  cout << (t = st2.nextToken()) << endl;
  if (t != "") return 1;
  cout << (t = st2.nextToken()) << endl;
  if (t != "") return 1;
  cout << (t = st2.nextToken()) << endl;
  if (t != "") return 1;
  cout << (t = st2.nextToken()) << endl;
  if (t != "eer") return 1;
 
  StringTokenizer st3(" aaazzer  aeerd a    eer", " \t", true, false);
  if (st3.numberOfRemainingTokens() != 1) return 1;
  cout << (t = st3.nextToken()) << endl;
  if (t != " aaazzer  aeerd a    eer") return 1;

  StringTokenizer st4(" aaazzer  aeerd a    eer", " \t", false, false);
  cout << st4.nextToken() << "+";
  cout << st4.unparseRemainingTokens() << endl;
  cout << st4.nextToken() << "+";
  cout << st4.unparseRemainingTokens() << endl;
  cout << st4.nextToken() << "+";
  cout << st4.unparseRemainingTokens() << endl;
  cout << st4.nextToken() << "+";
  cout << st4.unparseRemainingTokens() << endl;
  return 0;
}
