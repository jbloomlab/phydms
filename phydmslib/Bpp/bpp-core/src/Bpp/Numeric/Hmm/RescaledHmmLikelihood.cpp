//
// File: RescaledHmmLikelihood.cpp
// Created by: Julien Dutheil
// Created on: Fri Oct 26 11:57 2007
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

#include "RescaledHmmLikelihood.h"

#include "../../App/ApplicationTools.h"

// from the STL:
#include <iostream>
#include <algorithm>
using namespace bpp;
using namespace std;

RescaledHmmLikelihood::RescaledHmmLikelihood(
    HmmStateAlphabet* hiddenAlphabet,
    HmmTransitionMatrix* transitionMatrix,
    HmmEmissionProbabilities* emissionProbabilities,
    const std::string& prefix) throw (Exception):
  AbstractHmmLikelihood(),
  AbstractParametrizable(prefix),
  hiddenAlphabet_(hiddenAlphabet),
  transitionMatrix_(transitionMatrix),
  emissionProbabilities_(emissionProbabilities),
  likelihood_(),
  dLikelihood_(),
  d2Likelihood_(),
  backLikelihood_(),
  backLikelihoodUpToDate_(false),
  scales_(),
  dScales_(),
  d2Scales_(),
  logLik_(),
  breakPoints_(),
  nbStates_(),
  nbSites_()
{
  if (!hiddenAlphabet)        throw Exception("RescaledHmmLikelihood: null pointer passed for HmmStateAlphabet.");
  if (!transitionMatrix)      throw Exception("RescaledHmmLikelihood: null pointer passed for HmmTransitionMatrix.");
  if (!emissionProbabilities) throw Exception("RescaledHmmLikelihood: null pointer passed for HmmEmissionProbabilities.");
  if (!hiddenAlphabet_->worksWith(transitionMatrix->getHmmStateAlphabet()))
    throw Exception("RescaledHmmLikelihood: HmmTransitionMatrix and HmmEmissionProbabilities should point toward the same HmmStateAlphabet object.");
  if (!hiddenAlphabet_->worksWith(emissionProbabilities->getHmmStateAlphabet()))
    throw Exception("RescaledHmmLikelihood: HmmTransitionMatrix and HmmEmissionProbabilities should point toward the same HmmStateAlphabet object.");
  nbStates_ = hiddenAlphabet_->getNumberOfStates();
  nbSites_ = emissionProbabilities_->getNumberOfPositions();

  //Manage parameters:
  addParameters_(hiddenAlphabet_->getParameters());
  addParameters_(transitionMatrix_->getParameters());
  addParameters_(emissionProbabilities_->getParameters());

  //Init arrays:
  likelihood_.resize(nbSites_ * nbStates_);
  
  scales_.resize(nbSites_);
  
  //Compute:
  computeForward_();
}

void RescaledHmmLikelihood::setNamespace(const std::string& nameSpace)
{
  AbstractParametrizable::setNamespace(nameSpace);

  hiddenAlphabet_->setNamespace(nameSpace);
  transitionMatrix_->setNamespace(nameSpace);
  emissionProbabilities_->setNamespace(nameSpace);
}

void RescaledHmmLikelihood::fireParameterChanged(const ParameterList& pl)
{
  bool alphabetChanged    = hiddenAlphabet_->matchParametersValues(pl);
  bool transitionsChanged = transitionMatrix_->matchParametersValues(pl);
  bool emissionChanged    = emissionProbabilities_->matchParametersValues(pl);
  // these lines are necessary because the transitions and emissions can depend on the alphabet.
  // we could use a StateChangeEvent, but this would result in computing some calculations twice in some cases
  // (when both the alphabet and other parameter changed).
  if (alphabetChanged && !transitionsChanged) transitionMatrix_->setParametersValues(transitionMatrix_->getParameters());
  if (alphabetChanged && !emissionChanged) emissionProbabilities_->setParametersValues(emissionProbabilities_->getParameters());
  
  computeForward_();
  backLikelihoodUpToDate_=false;
}

/***************************************************************************************************************************/

void RescaledHmmLikelihood::computeForward_()
{
  double x;
  vector<double> tmp(nbStates_);
  vector<double> lScales(nbSites_);
  vector<double> trans(nbStates_ * nbStates_);

  //Transition probabilities:
  for (size_t i = 0; i < nbStates_; i++)
  {
    size_t ii = i * nbStates_;
    for (size_t j = 0; j < nbStates_; j++) {
      trans[ii + j] = transitionMatrix_->Pij(j, i);
      if (isnan(trans[ii + j]))
        throw Exception("RescaledHmmLikelihood::computeForward_. NaN transition probability");
      if (trans[ii + j] < 0)
        throw Exception("RescaledHmmLikelihood::computeForward_. Negative transition probability: " + TextTools::toString(trans[ii + j]));
    }
  }

  //Initialisation:
  scales_[0] = 0;
  const vector<double>* emissions = &(*emissionProbabilities_)(0);
  for (size_t j = 0; j < nbStates_; j++)
  {
    size_t jj = j * nbStates_;
    x = 0;
    for (size_t k = 0; k < nbStates_; k++)
    {
      x += trans[k + jj] * transitionMatrix_->getEquilibriumFrequencies()[k];
      //cerr << j << "\t" << k << "\t" << trans[k + jj] << "\t" << transitionMatrix_->getEquilibriumFrequencies()[k] << "\t" << trans[k + jj] * transitionMatrix_->getEquilibriumFrequencies()[k] << "\t" << x << endl;  
    }
    tmp[j] = (*emissions)[j] * x;
    //cerr << "e[j]=" << (*emissions)[j] << "\t" << tmp[j] << endl;
    scales_[0] += tmp[j];
  }
  for (size_t j = 0; j < nbStates_; j++)
  {
    likelihood_[j] = tmp[j] / scales_[0];
  }
  lScales[0] = log(scales_[0]);
 
  //Recursion:
  size_t nextBrkPt = nbSites_; //next break point
  vector<size_t>::const_iterator bpIt = breakPoints_.begin();
  if (bpIt != breakPoints_.end()) nextBrkPt = *bpIt;
  
  double a;
  for (size_t i = 1; i < nbSites_; i++)
  {
    size_t ii = i * nbStates_;
    size_t iip = (i - 1) * nbStates_;
    scales_[i] = 0 ;
    emissions = &(*emissionProbabilities_)(i);
    if (i < nextBrkPt)
    {
      for (size_t j = 0; j < nbStates_; j++)
      {
        size_t jj = j * nbStates_;
        x = 0;
        for (size_t k = 0; k < nbStates_; k++)
        {
          a = trans[jj + k] * likelihood_[iip + k];
          //if (a < 0)
          //{
          //  (*ApplicationTools::warning << "Negative value for likelihood at " << i << ", state " << j << ": " << likelihood_[iip + k] << ", Pij = " << trans[jj + k]).endLine();
          //  a = 0;
          //}
          x += a;
        }
        tmp[j] = (*emissions)[j] * x;
        if (tmp[j] < 0)
        {
          (*ApplicationTools::warning << "Negative probability at " << i << ", state " << j << ": " << (*emissions)[j] << "\t" << x).endLine();
          tmp[j] = 0;
        }
        scales_[i] += tmp[j];
      }
    }
    else //Reset markov chain:
    {
      for (size_t j = 0; j < nbStates_; j++)
      {
        size_t jj = j * nbStates_;
        x = 0;
        for (size_t k = 0; k < nbStates_; k++)
        {
          a = trans[jj + k] * transitionMatrix_->getEquilibriumFrequencies()[k];
          //if (a < 0)
          //{
          //  (*ApplicationTools::warning << "Negative value for likelihood at " << i << ", state " << j << ": ,Pij = " << trans[jj + k]).endLine();
          //  a = 0;
          //}
          x += a;
        }
        tmp[j] = (*emissions)[j] * x;
        //if (tmp[j] < 0)
        //{
        //  (*ApplicationTools::warning << "Negative emission probability at " << i << ", state " << j << ": " << (*emissions)[j]).endLine();
        //  tmp[j] = 0;
        //}
        scales_[i] += tmp[j];
      }
      bpIt++;
      if (bpIt != breakPoints_.end()) nextBrkPt = *bpIt;
      else nextBrkPt = nbSites_;
    }

    for (size_t j = 0; j < nbStates_; j++)
    {
      if (scales_[i] > 0) likelihood_[ii + j] = tmp[j] / scales_[i];
      else                likelihood_[ii + j] = 0;
    }
    lScales[i] = log(scales_[i]);
  }
  greater<double> cmp;
  sort(lScales.begin(), lScales.end(), cmp);
  logLik_ = 0;
  for (size_t i = 0; i < nbSites_; ++i)
  {
    logLik_ += lScales[i];
  }
}

/***************************************************************************************************************************/

void RescaledHmmLikelihood::computeBackward_() const
{
  if (backLikelihood_.size()==0)
    {
      backLikelihood_.resize(nbSites_);
      for (size_t i=0;i<nbSites_;i++)
        backLikelihood_[i].resize(nbStates_);
    }

  double x;

  //Transition probabilities:
  vector<double> trans(nbStates_ * nbStates_);
  for (size_t i = 0; i < nbStates_; i++)
  {
    size_t ii = i * nbStates_;
    for (size_t j = 0; j < nbStates_; j++)
      trans[ii + j] = transitionMatrix_->Pij(i, j);
  }


  //Initialisation:
  const vector<double>* emissions = 0;
  size_t nextBrkPt = 0; //next break point
  vector<size_t>::const_reverse_iterator bpIt = breakPoints_.rbegin();
  if (bpIt != breakPoints_.rend()) nextBrkPt = *bpIt;
  
  for (size_t j = 0; j < nbStates_; j++)
  {
    x = 0;
    backLikelihood_[nbSites_ - 1][j] = 1.;
  }

  //Recursion:
  for (size_t i = nbSites_ - 1; i > 0; i--)
  {
    emissions = &(*emissionProbabilities_)(i);
    if (i > nextBrkPt)
    {
      for (size_t j = 0; j < nbStates_; j++)
      {
        x = 0;
        size_t jj = j * nbStates_;
        for (size_t k = 0; k < nbStates_; k++)
        {
          x += (*emissions)[k] * trans[jj + k] * backLikelihood_[i][k];
        }
        backLikelihood_[i-1][j] = x / scales_[i];
      }    
    }
    else //Reset markov chain
    {
      for (size_t j = 0; j < nbStates_; j++)
      {
        backLikelihood_[i-1][j] = 1.;
      }    
      bpIt++;
      if (bpIt != breakPoints_.rend()) nextBrkPt = *bpIt;
      else nextBrkPt = 0;
    }
  }

  backLikelihoodUpToDate_=true;
}

/***************************************************************************************************************************/

double RescaledHmmLikelihood::getLikelihoodForASite(size_t site) const
{
  Vdouble probs=getHiddenStatesPosteriorProbabilitiesForASite(site);
  double x=0;
  for (size_t i=0;i<nbStates_;i++)
    x+=probs[i]*(*emissionProbabilities_)(site,i);

  return x;
}

Vdouble RescaledHmmLikelihood::getLikelihoodForEachSite() const
{
  std::vector< std::vector<double> > vv;
  getHiddenStatesPosteriorProbabilities(vv);

  Vdouble ret(nbSites_);
  for (size_t i=0;i<nbSites_;i++)
    {
      ret[i]=0;
      for (size_t j=0;j<nbStates_;j++)
        ret[i]+=vv[i][j]*(*emissionProbabilities_)(i,j);
    }

  return ret;
}

/***************************************************************************************************************************/

Vdouble RescaledHmmLikelihood::getHiddenStatesPosteriorProbabilitiesForASite(size_t site) const
{
  if (!backLikelihoodUpToDate_)
    computeBackward_();

  Vdouble probs(nbStates_);
  
  for (size_t j = 0; j < nbStates_; j++)
  {
    probs[j] = likelihood_[site * nbStates_ + j] * backLikelihood_[site][j];
  }

  return probs;
}


void RescaledHmmLikelihood::getHiddenStatesPosteriorProbabilities(std::vector< std::vector<double> >& probs, bool append) const throw (Exception)
{
  size_t offset = append ? probs.size() : 0;
  probs.resize(offset + nbSites_);
  for (size_t i = 0; i < nbSites_; i++)
    {
      probs[offset + i].resize(nbStates_);
    }

  if (!backLikelihoodUpToDate_)
    computeBackward_();
  
  for (size_t i = 0; i < nbSites_; i++)
    {
      size_t ii = i * nbStates_;
      for (size_t j = 0; j < nbStates_; j++)
        {
          probs[offset + i][j] = likelihood_[ii + j] * backLikelihood_[i][j];
        }
    }
}

/***************************************************************************************************************************/

void RescaledHmmLikelihood::computeDForward_() const
{
  //Init arrays:
  if (dLikelihood_.size()==0){
    dLikelihood_.resize(nbSites_);
    for (size_t i=0;i<nbSites_;i++)
      dLikelihood_[i].resize(nbStates_);
  }
  if (dScales_.size()==0)
    dScales_.resize(nbSites_);
  
  double x;
  vector<double> tmp(nbStates_), dTmp(nbStates_);
  vector<double> dLScales(nbSites_);
  
  //Transition probabilities:
  const ColMatrix<double> trans(transitionMatrix_->getPij());

  //Initialisation:
  dScales_[0] = 0;
  const vector<double>* emissions = &(*emissionProbabilities_)(0);
  const vector<double>* dEmissions = &emissionProbabilities_->getDEmissionProbabilities(0);

  for (size_t j = 0; j < nbStates_; j++)
  {
    dTmp[j] = (*dEmissions)[j] * transitionMatrix_->getEquilibriumFrequencies()[j];
    tmp[j] = (*emissions)[j] * transitionMatrix_->getEquilibriumFrequencies()[j];

    dScales_[0] += dTmp[j];
  }

  dLScales[0]=dScales_[0]/scales_[0];

  
  for (size_t j = 0; j < nbStates_; j++)
    dLikelihood_[0][j] = (dTmp[j] * scales_[0] - tmp[j] * dScales_[0]) / pow(scales_[0],2);
 
  //Recursion:

  size_t nextBrkPt = nbSites_; //next break point
  vector<size_t>::const_iterator bpIt = breakPoints_.begin();
  if (bpIt != breakPoints_.end()) nextBrkPt = *bpIt;
  
  for (size_t i = 1; i < nbSites_; i++)
  {
    size_t iip = (i - 1) * nbStates_;

    dScales_[i] = 0 ;

    emissions = &(*emissionProbabilities_)(i);
    dEmissions = &emissionProbabilities_->getDEmissionProbabilities(i);
    
    if (i < nextBrkPt)
    {
      for (size_t j = 0; j < nbStates_; j++)
      {
        x = 0;
        for (size_t k = 0; k < nbStates_; k++)
          x += trans(k,j) * likelihood_[iip + k];

        tmp[j] = (*emissions)[j] * x;
        dTmp[j] = (*dEmissions)[j] * x + (*emissions)[j] * VectorTools::sum(trans.getCol(j) * dLikelihood_[i-1]);
          
        dScales_[i] += dTmp[j];
      }
    }
    else //Reset markov chain:
    {
      for (size_t j = 0; j < nbStates_; j++)
      {
        dTmp[j] = (*dEmissions)[j] * transitionMatrix_->getEquilibriumFrequencies()[j];
        tmp[j] = (*emissions)[j] * transitionMatrix_->getEquilibriumFrequencies()[j];
        
        dScales_[i] += dTmp[j];
      }
      
      bpIt++;
      if (bpIt != breakPoints_.end()) nextBrkPt = *bpIt;
      else nextBrkPt = nbSites_;
    }

    dLScales[i]=dScales_[i]/scales_[i];

    for (size_t j = 0; j < nbStates_; j++)
      dLikelihood_[i][j] = (dTmp[j] * scales_[i] - tmp[j] * dScales_[i]) / pow(scales_[i],2);
  }
  
  greater<double> cmp;
  sort(dLScales.begin(), dLScales.end(), cmp);
  dLogLik_ = 0;
  for (size_t i = 0; i < nbSites_; ++i)
  {
    dLogLik_ += dLScales[i];
  }
}

double RescaledHmmLikelihood::getDLogLikelihoodForASite(size_t site) const
{
  return dScales_[site]/scales_[site];
}

/***************************************************************************************************************************/


void RescaledHmmLikelihood::computeD2Forward_() const
{
  //Init arrays:
  if (d2Likelihood_.size()==0){
    d2Likelihood_.resize(nbSites_);
    for (size_t i=0;i<nbSites_;i++)
      d2Likelihood_[i].resize(nbStates_);
  }
  if (d2Scales_.size()==0)
    d2Scales_.resize(nbSites_);
  
  double x;
  vector<double> tmp(nbStates_), dTmp(nbStates_), d2Tmp(nbStates_);
  vector<double> d2LScales(nbSites_);
  
  //Transition probabilities:
  const ColMatrix<double> trans(transitionMatrix_->getPij());

  //Initialisation:
  d2Scales_[0] = 0;
  const vector<double>* emissions = &(*emissionProbabilities_)(0);
  const vector<double>* dEmissions = &emissionProbabilities_->getDEmissionProbabilities(0);
  const vector<double>* d2Emissions = &emissionProbabilities_->getD2EmissionProbabilities(0);

  for (size_t j = 0; j < nbStates_; j++)
  {
    tmp[j] = (*emissions)[j] * transitionMatrix_->getEquilibriumFrequencies()[j];
    dTmp[j] = (*dEmissions)[j] * transitionMatrix_->getEquilibriumFrequencies()[j];
    d2Tmp[j] = (*d2Emissions)[j] * transitionMatrix_->getEquilibriumFrequencies()[j];

    d2Scales_[0] += d2Tmp[j];
  }

  d2LScales[0]=d2Scales_[0]/scales_[0]-pow(dScales_[0]/scales_[0],2);
  
  for (size_t j = 0; j < nbStates_; j++)
    d2Likelihood_[0][j] = d2Tmp[j] / scales_[0] - (d2Scales_[0] * tmp[j] + 2 * dScales_[0] * dTmp[j]) / pow(scales_[0],2)
      +  2 * pow(dScales_[0],2) * tmp[j] / pow(scales_[0],3);
   
  //Recursion:

  size_t nextBrkPt = nbSites_; //next break point
  vector<size_t>::const_iterator bpIt = breakPoints_.begin();
  if (bpIt != breakPoints_.end()) nextBrkPt = *bpIt;
  
  for (size_t i = 1; i < nbSites_; i++)
  {
    dScales_[i] = 0 ;

    emissions = &(*emissionProbabilities_)(i);
    dEmissions = &emissionProbabilities_->getDEmissionProbabilities(i);
    d2Emissions = &emissionProbabilities_->getD2EmissionProbabilities(i);
    
    if (i < nextBrkPt)
    {
      size_t iip = (i - 1) * nbStates_;

      for (size_t j = 0; j < nbStates_; j++)
      {
        x = 0;
        for (size_t k = 0; k < nbStates_; k++)
          x += trans(k,j) * likelihood_[iip + k];

        tmp[j] = (*emissions)[j] * x;
        dTmp[j] = (*dEmissions)[j] * x + (*emissions)[j] * VectorTools::sum(trans.getCol(j) * dLikelihood_[i-1]);
        d2Tmp[j] = (*d2Emissions)[j] * x + 2 * (*dEmissions)[j] * VectorTools::sum(trans.getCol(j) * dLikelihood_[i-1])
          + (*emissions)[j] * VectorTools::sum(trans.getCol(j) * d2Likelihood_[i-1]);
          
        d2Scales_[i] += d2Tmp[j];
      }
    }
    else //Reset markov chain:
    {
      for (size_t j = 0; j < nbStates_; j++)
      {
        tmp[j] = (*emissions)[j] * transitionMatrix_->getEquilibriumFrequencies()[j];
        dTmp[j] = (*dEmissions)[j] * transitionMatrix_->getEquilibriumFrequencies()[j];
        d2Tmp[j] = (*d2Emissions)[j] * transitionMatrix_->getEquilibriumFrequencies()[j];
        
        d2Scales_[i] += d2Tmp[j];
      }
      
      bpIt++;
      if (bpIt != breakPoints_.end()) nextBrkPt = *bpIt;
      else nextBrkPt = nbSites_;
    }

    d2LScales[i]=d2Scales_[i]/scales_[i]-pow(dScales_[i]/scales_[i],2);
  
    for (size_t j = 0; j < nbStates_; j++)
      d2Likelihood_[i][j] = d2Tmp[j] / scales_[i] - (d2Scales_[i] * tmp[j] + 2 * dScales_[i] * dTmp[j]) / pow(scales_[i],2)
        +  2 * pow(dScales_[i],2) * tmp[j] / pow(scales_[i],3);
  }
  
  greater<double> cmp;
  sort(d2LScales.begin(), d2LScales.end(), cmp);
  dLogLik_ = 0;
  for (size_t i = 0; i < nbSites_; ++i)
  {
    d2LogLik_ += d2LScales[i];
  }
}

/***************************************************************************************************************************/

double RescaledHmmLikelihood::getD2LogLikelihoodForASite(size_t site) const
{
  return d2Scales_[site]/scales_[site]-pow(dScales_[site]/scales_[site],2);
}

