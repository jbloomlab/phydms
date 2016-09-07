//  
// File: ExperimentallyInformedCodonModel.cpp
// Created by: Jesse Bloom
// Created on: May 2015
//

/*
  This file was created by extensively modifying the 
  CodonDistanceFitnessPhaseFrequenciesSubstitutionModel.cpp file
  distributed with Bio++
*/
 
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
  with loading,  using,  modifying and/or developi_ng or reproducing the
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


#include "ExperimentallyInformedCodonModel.h"
#include <Bpp/Numeric/NumConstants.h>
#include <Bpp/Phyl/Model/Nucleotide/K80.h>
#include <Bpp/Phyl/Model/FrequenciesSet/CodonFrequenciesSet.h>
#include <iostream>


bppextensions::ExperimentallyInformedCodonModel::ExperimentallyInformedCodonModel(
    const bpp::GeneticCode* gCode,
    bpp::FrequenciesSet* preferences,
    const std::string& prefix,
    bool prefsasparams,
    bool divpressure,
    double maxdeltar,
    double mindeltar,
    double deltar, 
    std::string fixationmodel) :
  AbstractParameterAliasable(prefix),
  AbstractCodonSubstitutionModel(gCode, new bpp::K80(dynamic_cast<const bpp::CodonAlphabet*>(gCode->getSourceAlphabet())->getNucleicAlphabet()), prefix),
  AbstractCodonPhaseFrequenciesSubstitutionModel(bpp::CodonFrequenciesSet::getFrequenciesSetForCodons(bpp::CodonFrequenciesSet::F1X4, gCode), prefix),
  prefix_(""),
  preferences_(preferences),
  omega_(1),
  omega2_(0),
  stringencyparameter_(1),
  rateparameter_(1),
  prefsasparams_(prefsasparams),
  deltar_(deltar),
  divpressure_(divpressure),
  fixationmodel_(fixationmodel)
{
  if (dynamic_cast<bpp::CodonFrequenciesSet*>(preferences) == NULL) {
    throw std::runtime_error("Invalid preferences");
  }
  prefix_ = prefix;
  prefName_ = "preferences_" + preferences_->getNamespace();
  preferences_->setNamespace(prefix + prefName_);
  if (prefsasparams_) {
    addParameters_(preferences_->getParameters());
  }
  float const smallnumber = 0.0001; // small number to avoid omega of zero, or omega * (1 + omega2 * delta) that is negative
  addParameter_(new bpp::Parameter(prefix + "omega", omega_, new bpp::IntervalConstraint(smallnumber, 99, true, true), true));
  if (divpressure){
  	if (maxdeltar - mindeltar < smallnumber){
  		throw std::runtime_error("Maximum diversifying selection not sufficiently larger than minimum diversifying selection. Rescale your diversiying selection numbers to be larger.\n");
  		}
  	else if (maxdeltar != mindeltar){
  		if (mindeltar >= 0){
  			addParameter_(new bpp::Parameter(prefix + "omega2", omega2_, new bpp::IntervalConstraint((-1.0/maxdeltar) + smallnumber, 99, true, true), true));
  		}
  		else if (maxdeltar <= 0){
			addParameter_(new bpp::Parameter(prefix + "omega2", omega2_, new bpp::IntervalConstraint(-99,(-1.0/mindeltar) - smallnumber, true, true), true));
  		}
  		else{
  			addParameter_(new bpp::Parameter(prefix + "omega2", omega2_, new bpp::IntervalConstraint((-1.0/maxdeltar) + smallnumber ,(-1.0/mindeltar) - smallnumber, true, true), true));
  			}
  	}
  }
  addParameter_(new bpp::Parameter(prefix + "stringencyparameter", 1, new bpp::IntervalConstraint(0.1, 10.0, true, true), true));
  f_gwF_ = 0.5;
  if (fixationmodel == "gwF") {
    addParameter_(new bpp::Parameter(prefix + "f_gwF", f_gwF_, new bpp::IntervalConstraint(0, 1, true, true), true));
    // we need sum to compute fixation probabilities for gwF method
    scaledprefsum_ = 0.0;
    for (size_t i = 0; i < preferences_->getNumberOfFrequencies(); i++) {
      double pi_i = std::pow(preferences_->getFrequencies()[i], stringencyparameter_);
      scaledprefsum_ += pi_i;
    }
  }
  updateMatrices();
}

bppextensions::ExperimentallyInformedCodonModel::~ExperimentallyInformedCodonModel()
{
  if (preferences_) delete preferences_;
}

std::string bppextensions::ExperimentallyInformedCodonModel::getName() const
{
  return prefix_;
}

void bppextensions::ExperimentallyInformedCodonModel::fireParameterChanged(const bpp::ParameterList& parameters)
{
  AbstractCodonPhaseFrequenciesSubstitutionModel::fireParameterChanged(parameters);
  omega_ = getParameterValue("omega");
  if(divpressure_){
  	omega2_ = getParameterValue("omega2");
  }
  stringencyparameter_ = getParameterValue("stringencyparameter");
  if (hasParameter("rateparameter")) {
      rateparameter_ = getParameterValue("rateparameter");
  }
  if (prefsasparams_) {
      preferences_->matchParametersValues(parameters);
  }
  if (fixationmodel_ == "gwF") {
    f_gwF_ = getParameterValue("f_gwF");
    // we need sum to compute fixation probabilities for gwF method
    scaledprefsum_ = 0.0;
    for (size_t i = 0; i < preferences_->getNumberOfFrequencies(); i++) {
      double pi_i = std::pow(preferences_->getFrequencies()[i], stringencyparameter_);
      scaledprefsum_ += pi_i;
    }
  }
  // this next call MUST be last!
  AbstractCodonSubstitutionModel::fireParameterChanged(parameters);
}

double bppextensions::ExperimentallyInformedCodonModel::getCodonsMulRate(size_t i, size_t j) const
{
  if (getGeneticCode()->areSynonymous(static_cast<int>(i), static_cast<int>(j))) {
    return rateparameter_ * 
      AbstractCodonSubstitutionModel::getCodonsMulRate(i,j)
      * AbstractCodonPhaseFrequenciesSubstitutionModel::getCodonsMulRate(i,j);
  } else {
    double fixationprob;
    double pi_i = std::pow(preferences_->getFrequencies()[i], stringencyparameter_);
    double pi_j = std::pow(preferences_->getFrequencies()[j], stringencyparameter_);
    if (fixationmodel_ == "HalpernBruno") {
      if (pi_j == pi_i) {
        fixationprob = 1;
      } else if (pi_i == 0) {
        fixationprob = 1000.0; // very large value if starting codon has zero preference
      } else if (pi_j == 0) {
        fixationprob = 0;
      } else {
        fixationprob = std::log(pi_j / pi_i) / (1 - (pi_i / pi_j));  // correct version of Halpern and Bruno (1998) equation; note that their paper has a typo
      }
    } else if (fixationmodel_ == "FracTolerated") {
        if (pi_j > pi_i) {
            fixationprob = 1.0;
        } else {
            fixationprob = pi_j / pi_i;
        }
    } else if (fixationmodel_ == "gwF") {
        if (pi_j == 0) {
          fixationprob = 0;
        } else if (pi_i == 0) {
          fixationprob = 1000.0; // very large value if starting codon has zero preference
        } else {
          fixationprob = std::pow(pi_j / scaledprefsum_, 1 - f_gwF_) / std::pow(pi_i / scaledprefsum_, f_gwF_);
        }
    } else {
      throw std::runtime_error("Invalid fixationmodel");
    }
    return omega_ * (1+omega2_*deltar_) * rateparameter_
      * AbstractCodonSubstitutionModel::getCodonsMulRate(i,j)
      * AbstractCodonPhaseFrequenciesSubstitutionModel::getCodonsMulRate(i,j)
      * fixationprob;
  }
}

void bppextensions::ExperimentallyInformedCodonModel::setNamespace(const std::string& st)
{
  AbstractCodonSubstitutionModel::setNamespace(st);
  AbstractParameterAliasable::setNamespace(st);
  AbstractCodonPhaseFrequenciesSubstitutionModel::setNamespace(st); 
  preferences_->setNamespace(st + prefName_);
}

std::map<std::string, double> bppextensions::ExperimentallyInformedCodonModel::getPreferences()
{
    std::map<std::string, double> prefs;
    for (size_t icodon = 0; icodon < preferences_->getNumberOfFrequencies(); icodon++) {
        prefs[preferences_->getAlphabet()->intToChar((int) icodon)] = preferences_->getFrequencies()[icodon];
    }
    return prefs;
}

std::string bppextensions::ExperimentallyInformedCodonModel::getPreferencesNamespace() {
    return preferences_->getNamespace();
}

void bppextensions::ExperimentallyInformedCodonModel::setFreq(std::map<int,double>& frequencies)
{
  throw std::runtime_error("It's not good that ExperimentallyInformedCodonModel->setFreq is being called, because I don't understand what this function does and it is probably not written correctly.\n");     
/*  AbstractCodonPhaseFrequenciesSubstitutionModel::setFreq(frequencies);
  map<int, double> freq1 = AbstractCodonPhaseFrequenciesSubstitutionModel::getFrequenciesSet()->getAlphabetStatesFrequencies();

  map<int, double> freq2;
  double s=0;
  map<int, double>::iterator it;

  for (it=frequencies.begin();it!=frequencies.end();it++)
  {
    freq2[it->first]=(freq1[it->first] != 0 ? it->second/freq1[it->first] : 0);
    s += freq2[it->first];
  }
  
  for (it = freq2.begin(); it != freq2.end(); it++)
    freq2[it->first] /= s;

  updateMatrices(); */
}


