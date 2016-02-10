//
// File: LogsumHmmLikelihood.cpp
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

#include "LogsumHmmLikelihood.h"

// from the STL:
#include <iostream>
#include <algorithm>
using namespace bpp;
using namespace std;

LogsumHmmLikelihood::LogsumHmmLikelihood(
    HmmStateAlphabet* hiddenAlphabet,
    HmmTransitionMatrix* transitionMatrix,
    HmmEmissionProbabilities* emissionProbabilities,
    bool ownsPointers,
    const std::string& prefix) throw (Exception):
  AbstractHmmLikelihood(),
  AbstractParametrizable(prefix),
  hiddenAlphabet_(hiddenAlphabet),
  transitionMatrix_(transitionMatrix),
  emissionProbabilities_(emissionProbabilities),
  ownsPointers_(ownsPointers),
  logLikelihood_(),
  partialLogLikelihoods_(),
  logLik_(),
  dLogLikelihood_(),
  partialDLogLikelihoods_(),
  d2LogLikelihood_(),
  partialD2LogLikelihoods_(),
  backLogLikelihood_(),
  backLogLikelihoodUpToDate_(false),
  breakPoints_(),
  nbStates_(),
  nbSites_()
{
  if (!hiddenAlphabet)        throw Exception("LogsumHmmLikelihood: null pointer passed for HmmStateAlphabet.");
  if (!transitionMatrix)      throw Exception("LogsumHmmLikelihood: null pointer passed for HmmTransitionMatrix.");
  if (!emissionProbabilities) throw Exception("LogsumHmmLikelihood: null pointer passed for HmmEmissionProbabilities.");

  if (!hiddenAlphabet_->worksWith(transitionMatrix_->getHmmStateAlphabet()))
    throw Exception("LogsumHmmLikelihood: HmmTransitionMatrix and HmmEmissionProbabilities should point toward the same HmmStateAlphabet object.");
  if (!hiddenAlphabet_->worksWith(emissionProbabilities_->getHmmStateAlphabet()))
    throw Exception("LogsumHmmLikelihood: HmmTransitionMatrix and HmmEmissionProbabilities should point toward the same HmmStateAlphabet object.");
  nbStates_ = hiddenAlphabet_->getNumberOfStates();
  nbSites_ = emissionProbabilities_->getNumberOfPositions();

  //Manage parameters:
  addParameters_(hiddenAlphabet_->getParameters());
  addParameters_(transitionMatrix_->getParameters());
  addParameters_(emissionProbabilities_->getParameters());

  //Init arrays:
  logLikelihood_.resize(nbSites_ * nbStates_);

  //Compute:
  computeForward_();
}

void LogsumHmmLikelihood::setNamespace(const std::string& nameSpace)
{
  AbstractParametrizable::setNamespace(nameSpace);

  hiddenAlphabet_->setNamespace(nameSpace);
  transitionMatrix_->setNamespace(nameSpace);
  emissionProbabilities_->setNamespace(nameSpace);
}

void LogsumHmmLikelihood::fireParameterChanged(const ParameterList& pl)
{
  dVariable_="";
  d2Variable_="";

  bool alphabetChanged    = hiddenAlphabet_->matchParametersValues(pl);
  bool transitionsChanged = transitionMatrix_->matchParametersValues(pl);
  bool emissionChanged    = emissionProbabilities_->matchParametersValues(pl);
  // these lines are necessary because the transitions and emissions can depend on the alphabet.
  // we could use a StateChangeEvent, but this would result in computing some calculations twice in some cases
  // (when both the alphabet and other parameter changed).
  if (alphabetChanged && !transitionsChanged) transitionMatrix_->setParametersValues(transitionMatrix_->getParameters());
  if (alphabetChanged && !emissionChanged) emissionProbabilities_->setParametersValues(emissionProbabilities_->getParameters());

  backLogLikelihoodUpToDate_=false;
  computeLikelihood();
}

void LogsumHmmLikelihood::computeLikelihood()
{
  computeForward_();
}

/***************************************************************************************************************************/

void LogsumHmmLikelihood::computeForward_()
{
  double x, a;
  vector<double> logTrans(nbStates_ * nbStates_);

  //Transition probabilities:
  for (size_t i = 0; i < nbStates_; i++)
  {
    size_t ii = i * nbStates_;
    for (size_t j = 0; j < nbStates_; j++)
      logTrans[ii + j] = log(transitionMatrix_->Pij(j, i));
  }

  //Initialisation:
  const vector<double>* emissions = &(* emissionProbabilities_)(0);

  for (size_t j = 0; j < nbStates_; j++)
  {
    size_t jj = j * nbStates_;
    x = logTrans[jj] + log(transitionMatrix_->getEquilibriumFrequencies()[0]);

    for (size_t k = 1; k < nbStates_; k++)
    {
      a = logTrans[k + jj] + log(transitionMatrix_->getEquilibriumFrequencies()[k]);
      x = NumTools::logsum(x, a);
    }

    logLikelihood_[j] = log((*emissions)[j]) + x;
  }
 
  //Recursion:
  size_t nextBrkPt = nbSites_; //next break point
  vector<size_t>::const_iterator bpIt = breakPoints_.begin();
  if (bpIt != breakPoints_.end()) nextBrkPt = *bpIt;
  partialLogLikelihoods_.clear();
 
  for (size_t i = 1; i < nbSites_; i++)
  {
    size_t ii = i * nbStates_;
    size_t iip = (i - 1) * nbStates_;
    emissions = &(*emissionProbabilities_)(i);
    if (i < nextBrkPt)
    {
      for (size_t j = 0; j < nbStates_; j++)
      {
        size_t jj = j * nbStates_;
        x = logTrans[jj] + logLikelihood_[iip];
        for (size_t k = 1; k < nbStates_; k++)
        {
          a = logTrans[jj + k] + logLikelihood_[iip + k];
          x = NumTools::logsum(x, a);
        }
        logLikelihood_[ii + j] = log((*emissions)[j]) + x;
      }
    }
    else //Reset markov chain:
    {
      //Termination of previous segment:
      double tmpLog = logLikelihood_[(i - 1) * nbStates_];
      for (size_t k = 1; k < nbStates_; k++)
        tmpLog = NumTools::logsum(tmpLog, logLikelihood_[(i - 1) * nbStates_ + k]);
      partialLogLikelihoods_.push_back(tmpLog);

      for (size_t j = 0; j < nbStates_; j++)
      {
        size_t jj = j * nbStates_;
        x = logTrans[jj] + log(transitionMatrix_->getEquilibriumFrequencies()[0]);
        for (size_t k = 1; k < nbStates_; k++)
        {
          a = logTrans[jj + k] + log(transitionMatrix_->getEquilibriumFrequencies()[k]);
          x = NumTools::logsum(x, a);
        }
        logLikelihood_[ii + j] = log((*emissions)[j]) + x;
      }
      bpIt++;
      if (bpIt != breakPoints_.end()) nextBrkPt = *bpIt;
      else nextBrkPt = nbSites_;
    }
  }

  //Termination:
  double tmpLog = logLikelihood_[(nbSites_ - 1) * nbStates_];
  for (size_t k = 1; k < nbStates_; k++)
    tmpLog = NumTools::logsum(tmpLog, logLikelihood_[(nbSites_ - 1) * nbStates_ + k]);
  partialLogLikelihoods_.push_back(tmpLog);
  
  //Compute likelihood:
  logLik_ = 0;
  vector<double> copy = partialLogLikelihoods_; //We need to keep the original order for posterior decoding.
  sort(copy.begin(), copy.end());
  for (size_t i = copy.size(); i > 0; --i)
    logLik_ += copy[i - 1];
}

/***************************************************************************************************************************/

void LogsumHmmLikelihood::computeBackward_() const
{
  if (backLogLikelihood_.size()==0)
  {
    backLogLikelihood_.resize(nbSites_);
    for (size_t i=0;i<nbSites_;i++)
      backLogLikelihood_[i].resize(nbStates_);
  }
  
  double x;

  //Transition probabilities:
  vector<double> logTrans(nbStates_ * nbStates_);
  for (size_t i = 0; i < nbStates_; i++)
  {
    size_t ii = i * nbStates_;
    for (size_t j = 0; j < nbStates_; j++)
      logTrans[ii + j] = log(transitionMatrix_->Pij(i, j));
  }


  //Initialisation:
  const vector<double>* emissions = 0;
  size_t nextBrkPt = 0; //next break point
  vector<size_t>::const_reverse_iterator bpIt = breakPoints_.rbegin();
  if (bpIt != breakPoints_.rend()) nextBrkPt = *bpIt;
  
  for (size_t k = 0; k < nbStates_; k++)
  {
    backLogLikelihood_[nbSites_ - 1][k] = 0.;
  }

  //Recursion:
  for (size_t i = nbSites_ - 1; i > 0; i--)
  {
    emissions = &(*emissionProbabilities_)(i);
    if (i > nextBrkPt)
    {
      for (size_t j = 0; j < nbStates_; j++)
      {
        size_t jj = j * nbStates_;
        x = log((*emissions)[0]) + logTrans[jj] + backLogLikelihood_[i][0];
        for (size_t k = 1; k < nbStates_; k++)
        {
          x = NumTools::logsum(x, log((*emissions)[k]) + logTrans[jj + k] + backLogLikelihood_[i][k]);
        }
        backLogLikelihood_[i - 1][j] = x;
      }    
    }
    else //Reset markov chain
    {
      for (unsigned int j = 0; j < nbStates_; j++)
      {
        backLogLikelihood_[i - 1][j] = 0.;
      }    
      bpIt++;
      if (bpIt != breakPoints_.rend()) nextBrkPt = *bpIt;
      else nextBrkPt = 0;
    }
  }

  backLogLikelihoodUpToDate_=true;
}


/***************************************************************************************************************************/

double LogsumHmmLikelihood::getLikelihoodForASite(size_t site) const
{
  Vdouble probs=getHiddenStatesPosteriorProbabilitiesForASite(site);
  double x=0;
  for (size_t i=0;i<nbStates_;i++)
    x+=probs[i]*(*emissionProbabilities_)(site,i);

  return x;
}

Vdouble LogsumHmmLikelihood::getLikelihoodForEachSite() const
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

Vdouble LogsumHmmLikelihood::getHiddenStatesPosteriorProbabilitiesForASite(size_t site) const
{
  if (!backLogLikelihoodUpToDate_)
    computeBackward_();

  Vdouble probs(nbStates_);
  
  vector<size_t>::const_iterator bpIt = breakPoints_.begin();
  vector<double>::const_iterator logLikIt = partialLogLikelihoods_.begin();
  while (bpIt != breakPoints_.end())
  {
    if (site>=(*bpIt))
      logLikIt++;
    else
      break;
    bpIt++;
  }

  for (size_t j = 0; j < nbStates_; j++)
  {
    probs[j] = exp(logLikelihood_[site * nbStates_ + j] + backLogLikelihood_[site][j] - *logLikIt);
  }

  return probs;
}

void LogsumHmmLikelihood::getHiddenStatesPosteriorProbabilities(std::vector< std::vector<double> >& probs, bool append) const throw (Exception)
{
  size_t offset = append ? probs.size() : 0;
  probs.resize(offset + nbSites_);
  for (size_t i = 0; i < nbSites_; i++)
  {
    probs[offset + i].resize(nbStates_);
  }

  if (!backLogLikelihoodUpToDate_)
    computeBackward_();
 
  size_t nextBrkPt = nbSites_; //next break point
  vector<size_t>::const_iterator bpIt = breakPoints_.begin();
  if (bpIt != breakPoints_.end()) nextBrkPt = *bpIt;
 
  vector<double>::const_iterator logLikIt = partialLogLikelihoods_.begin();
  for (size_t i = 0; i < nbSites_; i++)
  {
    if (i == nextBrkPt) {
      logLikIt++;
      bpIt++;
      if (bpIt != breakPoints_.end()) nextBrkPt = *bpIt;
      else nextBrkPt = nbSites_;
    }

    size_t ii = i * nbStates_;
    for (size_t j = 0; j < nbStates_; j++)
    {
      probs[offset + i][j] = exp(logLikelihood_[ii + j] + backLogLikelihood_[i][j] - *logLikIt);
    }
  }
}

/***************************************************************************************************************************/

void LogsumHmmLikelihood::computeDForward_() const
{
  //Init arrays:
  if (dLogLikelihood_.size()==0){
    dLogLikelihood_.resize(nbSites_);
    for (size_t i=0;i<nbSites_;i++)
      dLogLikelihood_[i].resize(nbStates_);
  }

  partialDLogLikelihoods_.clear();

  double x;

  vector<double> num(nbStates_);

  //Transition probabilities:
  const ColMatrix<double> trans(transitionMatrix_->getPij());

  //Initialisation:
  const vector<double>* emissions = &(*emissionProbabilities_)(0);
  const vector<double>* dEmissions = &emissionProbabilities_->getDEmissionProbabilities(0);
  
  for (size_t j = 0; j < nbStates_; j++)
    dLogLikelihood_[0][j] = (*dEmissions)[j] / (*emissions)[j];

  //Recursion:
  size_t nextBrkPt = nbSites_; //next break point
  vector<size_t>::const_iterator bpIt = breakPoints_.begin();
  if (bpIt != breakPoints_.end()) nextBrkPt = *bpIt;
  partialDLogLikelihoods_.clear();
 
  for (size_t i = 1; i < nbSites_; i++)
  {
    size_t iip = (i - 1) * nbStates_;

    emissions = &(*emissionProbabilities_)(i);
    dEmissions = &emissionProbabilities_->getDEmissionProbabilities(i);

    if (i < nextBrkPt)
    {
      for (size_t j = 0; j < nbStates_; j++)
      {
        x=(*dEmissions)[j]/(*emissions)[j];

        for (size_t k = 0; k < nbStates_; k++)
        {
          for (size_t kp = 0; kp < nbStates_; kp++)
            num[kp]=logLikelihood_[iip+kp]-logLikelihood_[iip+k];

          x+=dLogLikelihood_[i-1][k]*trans(k,j)/VectorTools::sumExp(num,trans.getCol(j));
        }
        
        dLogLikelihood_[i][j] = x;
      }
    }      
    else //Reset markov chain:
    {
      //Termination of previous segment
      x = 0;
      for (size_t k = 0; k < nbStates_; k++)
      {
        for (size_t kp = 0; kp < nbStates_; kp++)
          num[kp]=logLikelihood_[iip+kp]-logLikelihood_[iip+k];
        
        x += dLogLikelihood_[i-1][k] / VectorTools::sumExp(num);
      }
          
      partialDLogLikelihoods_.push_back(x);

      for (size_t j = 0; j < nbStates_; j++)
        dLogLikelihood_[i][j] = (*dEmissions)[j] / (*emissions)[j];
      
      bpIt++;
      if (bpIt != breakPoints_.end())
        nextBrkPt = *bpIt;
      else
        nextBrkPt = nbSites_;
    }
  }
  
  //Termination:
  x=0;
  for (size_t k = 0; k < nbStates_; k++)
  {
    for (size_t kp = 0; kp < nbStates_; kp++)
      num[kp]=logLikelihood_[nbStates_*(nbSites_-1)+kp]-logLikelihood_[nbStates_*(nbSites_-1)+k];
            
    x += dLogLikelihood_[nbSites_-1][k] / VectorTools::sumExp(num);
  }
          
  partialDLogLikelihoods_.push_back(x);
  
  //Compute dLogLikelihood
  
  dLogLik_ = 0;
  vector<double> copy = partialDLogLikelihoods_; //We need to keep the original order for posterior decoding.
  sort(copy.begin(), copy.end());
  for (size_t i = copy.size(); i > 0; --i)
    dLogLik_ += copy[i - 1];
}

double LogsumHmmLikelihood::getDLogLikelihoodForASite(size_t site) const
{
  return partialDLogLikelihoods_[site];
}

/***************************************************************************************************************************/

void LogsumHmmLikelihood::computeD2Forward_() const
{
  // Make sure that Dlikelihoods are correctly computed
  getFirstOrderDerivative(d2Variable_);
  
  //Init arrays:
  if (d2LogLikelihood_.size()==0){
    d2LogLikelihood_.resize(nbSites_);
    for (size_t i=0;i<nbSites_;i++)
      d2LogLikelihood_[i].resize(nbStates_);
  }

  partialD2LogLikelihoods_.clear();
  
  double x, z, snum;

  vector<double> num(nbStates_);
  
  //Transition probabilities:
  const ColMatrix<double> trans(transitionMatrix_->getPij());
  
  //Initialisation:
  const vector<double>* emissions = &(*emissionProbabilities_)(0);
  const vector<double>* dEmissions = &emissionProbabilities_->getDEmissionProbabilities(0);
  const vector<double>* d2Emissions = &emissionProbabilities_->getD2EmissionProbabilities(0);
  
  for (size_t j = 0; j < nbStates_; j++)
    d2LogLikelihood_[0][j] = (*d2Emissions)[j] / (*emissions)[j] - pow((*dEmissions)[j] / (*emissions)[j],2);

  //Recursion:
  size_t nextBrkPt = nbSites_; //next break point
  vector<size_t>::const_iterator bpIt = breakPoints_.begin();
  if (bpIt != breakPoints_.end()) nextBrkPt = *bpIt;
  partialDLogLikelihoods_.clear();
 
  for (size_t i = 1; i < nbSites_; i++)
  {
    size_t iip = (i - 1) * nbStates_;

    emissions = &(*emissionProbabilities_)(i);
    dEmissions = &emissionProbabilities_->getDEmissionProbabilities(i);
    d2Emissions = &emissionProbabilities_->getD2EmissionProbabilities(i);

    if (i < nextBrkPt)
    {
      for (size_t j = 0; j < nbStates_; j++)
      {
        x=(*d2Emissions)[j] / (*emissions)[j] - pow((*dEmissions)[j] / (*emissions)[j],2);

        for (size_t k = 0; k < nbStates_; k++)
        {
          for (size_t kp = 0; kp < nbStates_; kp++)
            num[kp]=logLikelihood_[iip+kp]-logLikelihood_[iip+k];
          snum=VectorTools::sumExp(num,trans.getCol(j));

          
          z=d2LogLikelihood_[i-1][k]+pow(dLogLikelihood_[i-1][k],2)
            - dLogLikelihood_[i-1][k] * VectorTools::sumExp(num, trans.getCol(j) * dLogLikelihood_[i-1])/snum;

          x += z * trans(k,j) / snum;
        }

        d2LogLikelihood_[i][j] = x;
      }
    }
    else //Reset markov chain:
    {
      x=0;
      
      //Termination of previous segment:
      for (size_t k = 1; k < nbStates_; k++)
      {
        for (size_t kp = 0; kp < nbStates_; kp++)
          num[kp]=logLikelihood_[iip+kp]-logLikelihood_[iip+k];
        
        snum=VectorTools::sumExp(num);
        
        x += (d2LogLikelihood_[i-1][k]+pow(dLogLikelihood_[i-1][k],2)
              - dLogLikelihood_[i-1][k] * VectorTools::sumExp(num, dLogLikelihood_[i-1])/snum)/snum;
      }
      
      partialD2LogLikelihoods_.push_back(x);
      
      for (size_t j = 0; j < nbStates_; j++)
        d2LogLikelihood_[i][j] = (*d2Emissions)[j] / (*emissions)[j] - pow((*dEmissions)[j] / (*emissions)[j],2);

      
      bpIt++;
      if (bpIt != breakPoints_.end()) nextBrkPt = *bpIt;
      else nextBrkPt = nbSites_;
    }
  }  

  //Termination:
  x=0;
  for (size_t k = 0; k < nbStates_; k++)
  {
    for (size_t kp = 0; kp < nbStates_; kp++)
      num[kp]=logLikelihood_[nbStates_*(nbSites_-1)+kp]-logLikelihood_[nbStates_*(nbSites_-1)+k];
    
    snum=VectorTools::sumExp(num);
    
    x += (d2LogLikelihood_[nbSites_-1][k]+pow(dLogLikelihood_[nbSites_-1][k],2)
          - dLogLikelihood_[nbSites_-1][k] * VectorTools::sumExp(num, dLogLikelihood_[nbSites_-1])/snum)/snum;
  }

  partialD2LogLikelihoods_.push_back(x);
  
  //Compute dLogLikelihood
  
  d2LogLik_ = 0;
  vector<double> copy = partialD2LogLikelihoods_; //We need to keep the original order for posterior decoding.
  sort(copy.begin(), copy.end());
  for (size_t i = copy.size(); i > 0; --i)
    d2LogLik_ += copy[i - 1];
}

double LogsumHmmLikelihood::getD2LogLikelihoodForASite(size_t site) const
{
  return partialD2LogLikelihoods_[site];
}
