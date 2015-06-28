//
// File: PolynomialFunction.cpp
// Created by: Julien Dutheil
// Created on: Wed Oct 27 18:46 2010
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

#include <Bpp/Numeric/Function/Functions.h>
#include <Bpp/Numeric/AbstractParametrizable.h>

#include <map>
#include <string>

using namespace bpp;
using namespace std;

class PolynomialFunction1:
  public virtual Function,
  public AbstractParametrizable
{
  private:
    double fval_;
 
  public:
    PolynomialFunction1() : AbstractParametrizable(""), fval_(0) {
      //We declare parameters here:
      addParameter_(new Parameter("x", 0));
      addParameter_(new Parameter("y", 0));
      addParameter_(new Parameter("z", 0));
      fireParameterChanged(getParameters());
    }
 
    PolynomialFunction1* clone() const { return new PolynomialFunction1(*this); }
 
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
      double z = getParameterValue("z");
      fval_ = (x-5)*(x-5) + (y+2)*(y+2) + (z-3)*(z-3);
    }
};

class PolynomialFunction1Der1:
  public PolynomialFunction1,
  public virtual DerivableFirstOrder
{
  protected:
    bool compFirstDer_;
    mutable map<string, double> firstDer_;

  public:
    PolynomialFunction1Der1(): compFirstDer_(true), firstDer_() {
      //Need to compute derivatives:
      fireParameterChanged(getParameters());
    }
 
    PolynomialFunction1Der1* clone() const { return new PolynomialFunction1Der1(*this); }
 
  public:
    void setParameters(const ParameterList& pl) 
        throw (ParameterNotFoundException, ConstraintException, Exception)
    {
      matchParametersValues(pl);
    }
 
    void fireParameterChanged(const ParameterList& pl) {
      PolynomialFunction1::fireParameterChanged(pl);
      if (compFirstDer_) {
        double x = getParameterValue("x");
        double y = getParameterValue("y");
        double z = getParameterValue("z");
        firstDer_["x"] = 2 * (x - 5);
        firstDer_["y"] = 2 * (y + 2);
        firstDer_["z"] = 2 * (z - 3);
      }
    }

    void enableFirstOrderDerivatives(bool yn) { compFirstDer_ = yn; }
    bool enableFirstOrderDerivatives() const { return compFirstDer_; }

    double getFirstOrderDerivative(const std::string& variable) const throw (Exception) {
      if (!compFirstDer_)
        throw Exception("PolynomialFunction1Der1::getFirstOrderDerivative. First order derivatives are not computed.");
      return firstDer_[variable];
    }
};

