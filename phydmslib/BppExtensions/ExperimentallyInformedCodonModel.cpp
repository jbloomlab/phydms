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


bppextensions::ExperimentallyInformedCodonModel::ExperimentallyInformedCodonModel(
    const bpp::GeneticCode* gCode,
    bpp::FrequenciesSet* preferences,
    const std::string& prefix,
    bool prefsasparams) :
  AbstractParameterAliasable(prefix),
  AbstractCodonSubstitutionModel(gCode, new bpp::K80(dynamic_cast<const bpp::CodonAlphabet*>(gCode->getSourceAlphabet())->getNucleicAlphabet()), prefix),
  AbstractCodonPhaseFrequenciesSubstitutionModel(bpp::CodonFrequenciesSet::getFrequenciesSetForCodons(bpp::CodonFrequenciesSet::F1X4, gCode), prefix),
  prefix_(""),
  preferences_(preferences),
  omega_(1),
  stringencyparameter_(1),
  rateparameter_(1),
  prefsasparams_(prefsasparams)
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
  addParameter_(new bpp::Parameter(prefix + "omega", 1, new bpp::IntervalConstraint(0.001, 99, true, true), true));
  addParameter_(new bpp::Parameter(prefix + "stringencyparameter", 1, new bpp::IntervalConstraint(0.1, 10.0, true, true), true));
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
  stringencyparameter_ = getParameterValue("stringencyparameter");
  if (hasParameter("rateparameter")) {
      rateparameter_ = getParameterValue("rateparameter");
  }
  if (prefsasparams_) {
      preferences_->matchParametersValues(parameters);
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
    if (pi_j == pi_i) {
      fixationprob = 1;
    } else if (pi_i == 0) {
      fixationprob = 1000.0; // very large value of starting codon has zero preference
    } else if (pi_j == 0) {
      fixationprob = 0;
    } else {
      fixationprob = std::log(pi_j / pi_i) / (1 - (pi_i / pi_j));  // correct version of Halpern and Bruno (1998) equation; note that their paper has a typo
    }
    return omega_ * rateparameter_
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


