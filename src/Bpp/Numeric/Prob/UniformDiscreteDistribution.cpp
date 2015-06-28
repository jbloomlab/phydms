//
// File: UniformDiscreteDistribution.cpp
// Created by: Laurent Guéguen
// Created on: April 2010
//

/*
Copyright or © or Copr. CNRS, (2010)

This software is a computer program whose purpose is to provide classes
for numerical calculus.

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

#include "UniformDiscreteDistribution.h"
#include "../Random/RandomTools.h"
#include "../../Utils/MapTools.h"

using namespace bpp;

// From the STL:
#include <cmath>

using namespace std;

/** Constructor: **************************************************************/

UniformDiscreteDistribution::UniformDiscreteDistribution(unsigned int n, double min, double max) :
  AbstractParameterAliasable("Uniform."),
  AbstractDiscreteDistribution(n,"Uniform."),
  min_((min<max)?min:max),
  max_((min<max)?max:min) 
{
  intMinMax_.setLowerBound(min_,false);
  intMinMax_.setUpperBound(max_,false);

  discretize();
}

UniformDiscreteDistribution::UniformDiscreteDistribution(const UniformDiscreteDistribution& udd) : 
  AbstractParameterAliasable(udd),
  AbstractDiscreteDistribution(udd),
  min_(udd.min_),
  max_(udd.max_)
{
}

UniformDiscreteDistribution& UniformDiscreteDistribution::operator=(const UniformDiscreteDistribution& udd) 
{
  AbstractParameterAliasable::operator=(udd);
  AbstractDiscreteDistribution::operator=(udd);
  min_=udd.min_;
  max_=udd.max_;

  return *this;
}

UniformDiscreteDistribution::~UniformDiscreteDistribution() {}

/******************************************************************************/

void UniformDiscreteDistribution::fireParameterChanged(const ParameterList& parameters)
{
}





