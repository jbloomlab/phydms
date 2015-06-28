//
// File: test_derivative1.cpp
// Created by: Julien Dutheil
// Created on: Thu Oct 28 12:49 2010
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

#include <Bpp/Numeric/Random/RandomTools.h>
#include <Bpp/Numeric/Function/PowellMultiDimensions.h>
#include <Bpp/Numeric/Function/ReparametrizationFunctionWrapper.h>
#include <vector>
#include <iostream>

using namespace bpp;
using namespace std;

class MyFunction:
  public virtual Function,
  public AbstractParametrizable
{
  private:
    double fval_;
 
  public:
    MyFunction() : AbstractParametrizable(""), fval_(0) {
      //We declare parameters here:
      addParameter_(new Parameter("x", 0, new IntervalConstraint(-1, 7, true, true), true));
      addParameter_(new Parameter("y", 0, new IntervalConstraint(-4, 4, true, true), true));
      fireParameterChanged(getParameters());
    }
 
    MyFunction* clone() const { return new MyFunction(*this); }
 
  public:
    void setParameters(const ParameterList& pl) 
        throw (ParameterNotFoundException, ConstraintException, Exception)
    {
      matchParametersValues(pl);
    }
    double getValue() const throw (Exception) { return fval_; }
 
    void fireParameterChanged(const ParameterList& pl) {
      double x = getParameterValue("x");
      double y = getParameterValue("y");
      fval_ = cos(x) + sin(y);
    }
};

int main() {
  MyFunction f;
  ReparametrizationFunctionWrapper fw(&f);
  ParameterList pl = fw.getParameters();
  PowellMultiDimensions optimizer(&fw);
  optimizer.init(pl);
  optimizer.optimize();
  double minf = f.getValue();
  double x = f.getParameterValue("x");
  double y = f.getParameterValue("y");
  cout << "x=" << x << endl;
  cout << "y=" << y << endl;
  cout << "f=" << minf << endl;

  cout << setprecision(20) << (abs(x - 3.141593) + abs(y + 1.570796)) << endl;
  bool test = (abs(x - 3.141593) + abs(y + 1.570796)) < optimizer.getStopCondition()->getTolerance();
  return (test ? 0 : 1);
}
