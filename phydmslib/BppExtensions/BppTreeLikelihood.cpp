/* BppTreeLikelihood.cpp
   Implements C++ class BppTreeLikelihood
   Written by Jesse Bloom
*/

#include <iostream>
#include <fstream>
#include <exception>
#include <string>
#include <sstream>
#include <cmath>
#include <Bpp/Phyl/Model/FrequenciesSet/FrequenciesSet.h>
#include <Bpp/Phyl/Io/BppOFrequenciesSetFormat.h>
#include <Bpp/Phyl/Model/Codon/YN98.h>
#include <Bpp/Phyl/Model/Codon/YNGKP_M1.h>
#include <Bpp/Phyl/Model/Codon/YNGKP_M2.h>
#include <Bpp/Phyl/Model/Codon/YNGKP_M3.h>
#include <Bpp/Phyl/Model/Codon/YNGKP_M7.h>
#include <Bpp/Phyl/Model/Codon/YNGKP_M8.h>
#include "BppTreeLikelihood.h"
#include "ExperimentallyInformedCodonModel.h"

// This function is a patch for the fact that there is a bug in some g++ so that to_string isn't included.
// See: http://stackoverflow.com/questions/12975341/to-string-is-not-a-member-of-std-says-so-g
namespace patch
{
    template < typename T > std::string to_string(const T& n)
    {
        std::ostringstream stm;
        stm << n;
        return stm.str();
    }
}


// constructor
bppextensions::BppTreeLikelihood::BppTreeLikelihood(std::vector<std::string> seqnames, std::vector<std::string> seqs, std::string treefile, std::string modelstring, int infertopology, std::map<int, std::map<std::string, double> > preferences, std::map<std::string, double> fixedmodelparams, std::map<std::string, double> initializemodelparams, int oldlikelihoodmethod, int fixbrlen, int addrateparameter, int prefsasparams, char recursion, int useLog, int ngammarates, int ncats)
{

    // setup some parameters / options
    oldlikmethod = oldlikelihoodmethod != 0;
    verbose = false;
    nrates = ngammarates;
    optimizationparams["optimization.reparametrization"] = "false";
    optimizationparams["optimization.profiler"] = "none";
    optimizationparams["optimization.backup.file"] = "none";
    optimizationparams["optimization.message_handler"] = "none";
    optimizationparams["optimization.verbose"] = "false";
    if (infertopology) {
        optimizationparams["optimization.topology"] = "true";
    } else {
        optimizationparams["optimization.topology"] = "false";
    }
    if (prefsasparams) {
        optimizationparams["optimization.ignore_parameters"] = "*preferences*";
    } else {
        optimizationparams["optimization.ignore_parameters"] = "";
    }

    // error checking on calling variables
    if ((fixbrlen) && (infertopology)) {
        throw std::runtime_error("Cannot use both infertopology and fixbrlen");
    }
    if ((addrateparameter) && (! fixbrlen)) {
        throw std::runtime_error("Cannot use addrateparameter without fixbrlen");
    }

    // Initialize to DNA alphabet of codons using standard genetic code
    ntalphabet = new bpp::DNA(false); // false indicates no gap characters other than -
    alphabet = new bpp::CodonAlphabet(ntalphabet);
    gcode = new bpp::StandardGeneticCode(ntalphabet);

    // Initialize sites
    std::vector<std::string>::size_type nseqs, i;
    nseqs = seqnames.size();
    if (nseqs != seqs.size()) {
        throw std::invalid_argument("seqnames and seqs of different sizes");
    }
    sites = new bpp::VectorSiteContainer(alphabet);
    bpp::BasicSequence *iseq;
    for (i = 0; i < nseqs; i++) {
        iseq = new bpp::BasicSequence(seqnames[i], seqs[i], alphabet);
        sites->addSequence(*dynamic_cast <bpp::Sequence*> (iseq));
        delete iseq;
    }
    bpp::SiteContainerTools::changeGapsToUnknownCharacters(*sites);
    long nsites = sites->getNumberOfSites();
    if (nsites < 1) {
        throw std::runtime_error("No sites are specified");
    }

    // initialize some values
    sharedtreeindex = 0;
    sharedrateindex = 0;
    sharedmodelindex = 0;
    sequenceevolution = 0;
    substitutionprocesscollection = 0;
    tree = 0;
    paramtree = 0;

    // read in tree
    treeReaderWriter = new bpp::Newick(true); // true indicates comments allowed in brackets
    tree = treeReaderWriter->read(treefile);
    if (nseqs != tree->getLeavesNames().size()) {
        throw std::invalid_argument("seqs and treefile specify different size sequence sets");
    }
    if (fixbrlen) {
        if (optimizationparams["optimization.ignore_parameters"].empty()) {
            optimizationparams["optimization.ignore_parameters"] = "*BrLen*"; 
        } else {
            optimizationparams["optimization.ignore_parameters"] = optimizationparams["optimization.ignore_parameters"] + "," + "*BrLen*";
        }
    }

    // set up substitution models
    if ((modelstring.length() >= 12) && (modelstring.substr(0, 6) == "YNGKP_")) {
        if (modelstring.substr(modelstring.length() - 8, 8) == "_empF3X4") {
            std::string thetaparams = "*_Full.theta*";
            if (optimizationparams["optimization.ignore_parameters"].empty()) {
                optimizationparams["optimization.ignore_parameters"] = thetaparams;
            } else {
                optimizationparams["optimization.ignore_parameters"] = optimizationparams["optimization.ignore_parameters"] + "," + thetaparams;
            }
        }
        else if (modelstring.substr(modelstring.length() - 8, 8) != "_fitF3X4") {
            throw std::invalid_argument("Invalid YNGKP frequencies method");
        }
        // These next few lines parallel bpp::PhylogeneticsApplicationTools::getSubstitutionModel to set up models
        bpp::BppOFrequenciesSetFormat freqReader(bpp::BppOFrequenciesSetFormat::ALL, verbose, 1);
        freqReader.setGeneticCode(gcode);
        auto_ptr<bpp::FrequenciesSet> codonFreqs(freqReader.read(alphabet, "F3X4", sites, true));
        // now set up the models
        if (modelstring.substr(6, 2) == "M0") {
            models[sharedmodelindex] = dynamic_cast<bpp::SubstitutionModel*>(new bpp::YN98(gcode, codonFreqs.release()));
        }
        else if (modelstring.substr(6, 2) == "M1") {
            models[sharedmodelindex] = dynamic_cast<bpp::SubstitutionModel*>(new bpp::YNGKP_M1(gcode, codonFreqs.release()));
        }
        else if (modelstring.substr(6, 2) == "M2") {
            models[sharedmodelindex] = dynamic_cast<bpp::SubstitutionModel*>(new bpp::YNGKP_M2(gcode, codonFreqs.release()));
        }
        else if (modelstring.substr(6, 2) == "M3") {
            models[sharedmodelindex] = dynamic_cast<bpp::SubstitutionModel*>(new bpp::YNGKP_M3(gcode, codonFreqs.release(), 3));
        }
        else if (modelstring.substr(6, 2) == "M7") {
            models[sharedmodelindex] = dynamic_cast<bpp::SubstitutionModel*>(new bpp::YNGKP_M7(gcode, codonFreqs.release(), ncats));
        }
        else if (modelstring.substr(6, 2) == "M8") {
            models[sharedmodelindex] = dynamic_cast<bpp::SubstitutionModel*>(new bpp::YNGKP_M8(gcode, codonFreqs.release(), ncats));
        } else {
            throw std::invalid_argument("Invalid model variant of YNGKP");
        }
        if (! models[sharedmodelindex]) {
            throw std::runtime_error("Error casting sharedmodel");
        }
        if (nsites == 1) {
            // add a pseudocount if only one site, as some frequencies are likely to be near zero
            models[sharedmodelindex]->setFreqFromData(*sites, 0.1);
        }
        else {
            models[sharedmodelindex]->setFreqFromData(*sites, 0);
        }
        if (addrateparameter) {
            models[sharedmodelindex]->addRateParameter();
        }
    } 
    else if (modelstring == "ExpCM") {
        if (oldlikmethod) {
            throw runtime_error("Cannot use experimentally informed models with old likelihood method");
        }
        if (nsites != (long) preferences.size()) {
            throw runtime_error("number of sites isn't equal to number of preferences");
        }
        std::map<int, double> init_rprefs;
        std::string codon;
        for (long isite = 1; isite <= nsites; isite++) {
            bpp::FullCodonFrequenciesSet *rprefs = new bpp::FullCodonFrequenciesSet(gcode);
            init_rprefs.clear();
            for (size_t icodon = 0; icodon < rprefs->getNumberOfFrequencies(); icodon++) {
                codon = rprefs->getAlphabet()->intToChar((int) icodon);
                if (preferences[isite].find(codon) == preferences[isite].end()) {
                    throw std::runtime_error("Failed to find codon " + codon + "\n");
                }
                init_rprefs[icodon] = preferences[isite][codon];    
            }
            rprefs->setFrequenciesFromAlphabetStatesFrequencies(init_rprefs);
            models[isite] = dynamic_cast<bpp::SubstitutionModel*>(new bppextensions::ExperimentallyInformedCodonModel(gcode, rprefs, "ExpCM.", prefsasparams != 0));
            if (! models[isite]) {
                throw std::runtime_error("error casting ExperimentallyInformedCodonModel");
            }
            if (addrateparameter) {
                models[isite]->addRateParameter();
            }
        }
    }
    else {
        throw std::invalid_argument("Invalid modelstring");
    }

    // fix and initialize model parameters
    if ((! fixedmodelparams.empty()) || (! initializemodelparams.empty())) {
//        if (oldlikmethod) {
//            throw std::runtime_error("fixedmodelparams and initializemodelparams not guaranteed to work if using oldlikmethod");
//        }
        for (std::map<std::string, double>::iterator itr = fixedmodelparams.begin(); itr != fixedmodelparams.end(); itr++) {
            if (initializemodelparams.find(itr->first) != initializemodelparams.end()) {
                throw std::runtime_error("Cannot specify a key in both fixedmodelparams and initializemodelparams, but you did for " + itr->first);
            }
            for (std::map<size_t, bpp::SubstitutionModel*>::iterator imodel_itr = models.begin(); imodel_itr != models.end(); ++imodel_itr) {
                if (imodel_itr->second->hasParameter(itr->first)) {
                    imodel_itr->second->setParameterValue(itr->first, itr->second);
                } else if (itr->first == "Gamma.alpha") {
                    if (nrates <= 1) {
                        throw std::runtime_error("You shouldn't have a Gamma.alpha parameter as you don't appear to have set the option for gamma distributed rates");
                    }
                } else {
                    // note that Gamma.alpha is set in initialization of GammaDiscreteRateDistribution
                    throw std::runtime_error("Cannot find parameter in fixedmodelparams: " + itr->first);
                }
            }
            if (optimizationparams["optimization.ignore_parameters"].empty()) {
                optimizationparams["optimization.ignore_parameters"] = "*" + itr->first + "*";
            } else {
                optimizationparams["optimization.ignore_parameters"] = optimizationparams["optimization.ignore_parameters"] + "," + "*" + itr->first + "*";
            }
        }
        for (std::map<std::string, double>::iterator itr = initializemodelparams.begin(); itr != initializemodelparams.end(); itr++) {
            for (std::map<size_t, bpp::SubstitutionModel*>::iterator imodel_itr = models.begin(); imodel_itr != models.end(); ++imodel_itr) {
                if (imodel_itr->second->hasParameter(itr->first)) {
                    imodel_itr->second->setParameterValue(itr->first, itr->second);
                } else {
                    throw std::runtime_error("Cannot find parameter in initializemodelparams: " + itr->first);
                }
            }
        }
        for (std::map<size_t, bpp::SubstitutionModel*>::iterator imodel_itr = models.begin(); imodel_itr != models.end(); ++imodel_itr) {
            bpp::AbstractParametrizable* imodel = dynamic_cast<bpp::AbstractParametrizable*>(imodel_itr->second);
            if (! imodel) {
                throw runtime_error("Failed to cast imodel");
            }
            imodel->fireParameterChanged(imodel->getParameters());
        }
    } 

    // set up rate distribution
    if (ngammarates == 1) {
        ratedistribution = new bpp::ConstantRateDistribution(); 
    } else if (ngammarates > 1) {
        if (fixedmodelparams.find("Gamma.alpha") == fixedmodelparams.end()) {
            ratedistribution = new bpp::GammaDiscreteRateDistribution(ngammarates);
        } else {
            ratedistribution = new bpp::GammaDiscreteRateDistribution(ngammarates, fixedmodelparams["Gamma.alpha"]);
        }
    } else {
        throw std::runtime_error("ngammrates must be integer >= 1");
    }

    // compute likelihoods
    oldtreelikelihood = 0;
    phylolikelihood = 0;
    if (! oldlikmethod) {
        if (infertopology) {
            throw std::runtime_error("Cannot use infertopology with new likelihood method");
        }
        // combine all of these into collection, similar to what is done by bpp::PhylogeneticsApplicationTools::getSubstitutionProcessCollection
        substitutionprocesscollection = new bpp::SubstitutionProcessCollection();
        paramtree = new bpp::ParametrizableTree(*tree);
        substitutionprocesscollection->addTree(paramtree, sharedtreeindex); 
        substitutionprocesscollection->addDistribution(ratedistribution, sharedrateindex); 
        std::vector<size_t> processbysite;
        if ((models.size() == 1) && (modelstring != "ExpCM")) { // all models the same
            substitutionprocesscollection->addModel(models[sharedmodelindex], sharedmodelindex);
            map<size_t, vector<int> > mModBr; // keyed by model, values are branches to which model is assigned
            for (int ibranch=0; ibranch<(int) substitutionprocesscollection->getTree(sharedtreeindex).getNumberOfBranches(); ibranch++) {
                mModBr[sharedmodelindex].push_back(ibranch);
            }
            substitutionprocesscollection->addSubstitutionProcess(sharedmodelindex, mModBr, sharedtreeindex, sharedrateindex);
            for (long isite = 1; isite <= nsites; isite++) {
                processbysite.push_back(sharedmodelindex);
            }
            // all model parameters are the same for each site. We add them to constrainedparams so names are returned correctly in ModelParams
            std::vector<std::string> independentmodelparams = models[sharedmodelindex]->getIndependentParameters().getParameterNames();
            for (std::vector<std::string>::size_type i = 0; i != independentmodelparams.size(); ++i) {
                constrainedparams[independentmodelparams[i]] = independentmodelparams[i]; // name to return in ModelParams
            }
        } else { // a different model for each site
            for (long isite = 1; isite <= nsites; isite++) {
                substitutionprocesscollection->addModel(models[isite], isite);
            }
            for (long isite = 1; isite <= nsites; isite++) {
            // assigns model isite (along with shared tree and rate) to process isite. This is mimicking the bpp::PhylogeneticsApplicationTools::addSubstitutionProcessCollectionMember function
                map<size_t, vector<int> > mModBr; // keyed by model, values are branches to which model is assigned
                for (int ibranch=0; ibranch<(int) substitutionprocesscollection->getTree(sharedtreeindex).getNumberOfBranches(); ibranch++) {
                    mModBr[isite].push_back(ibranch);
                }
                substitutionprocesscollection->addSubstitutionProcess(isite, mModBr, sharedtreeindex, sharedrateindex);
            }
            // now alias parameters for substitution models, similar to what is done at end of bpp::PhylogeneticsApplicationTools::getSubstitutionProcessCollection
            // This constrains the parameters to be the same for all models
            // all parameters are aliased to that for model1
            std::vector<std::string> independentmodelparams = models[1]->getIndependentParameters().getParameterNames();
            for (std::vector<std::string>::size_type i = 0; i != independentmodelparams.size(); ++i) {
                std::string parametername = independentmodelparams[i] + "_1";
                if (parametername.find("preferences") != std::string::npos) {
                    continue;
                }
                constrainedparams[independentmodelparams[i]] = independentmodelparams[i]; // name to return in ModelParams
                for (long isite = 2; isite <= nsites; isite++) {
                    std::string parameternametoalias = independentmodelparams[i] + "_" + patch::to_string(isite);
                    substitutionprocesscollection->aliasParameters(parametername, parameternametoalias);
                }
            }
            for (long isite = 1; isite <= nsites; isite++) {
                processbysite.push_back(isite);
            }
        }
        // set sequenceevolution (partitioned sites), this mimics what is done by bpp::PhylogeneticsApplicationTools::getSequenceEvolutions for a Partition
        sequenceevolution = new bpp::PartitionSequenceEvolution(substitutionprocesscollection, processbysite);

        // get phylogenetic likelihood, this mimics what is done by bpp::PhylogeneticsApplicationTools::getPhyloLikelihoods
        if ((recursion != 'S') && (recursion != 'D')) {
            throw std::runtime_error("Invalid value of recursion, should be S or D");
        }
        bpp::PartitionSequenceEvolution* pse = dynamic_cast<bpp::PartitionSequenceEvolution*>(sequenceevolution);
        if (pse == NULL) {
            throw std::runtime_error("Failed to cast to PartitionSequenceEvolution");
        }
        bool patterns = true; // is default, not sure exactly what it means
        phylolikelihood = new bpp::PartitionProcessPhyloLikelihood(*sites, *pse, 0, 0, verbose, patterns);
        phylolikelihood->setUseLog((bool) useLog);

    } else { // use the old likelihood method
        if (models.find(sharedmodelindex) == models.end()) {
            throw runtime_error("Cannot use old likelihood method as there is not a shared model");
        }
        if (dynamic_cast<bpp::MixedSubstitutionModel*>(models[sharedmodelindex]) == 0) {
            if (infertopology) {
                oldtreelikelihood = new bpp::NNIHomogeneousTreeLikelihood(*tree, *sites, models[sharedmodelindex], ratedistribution, true, verbose);
            } else {
                oldtreelikelihood = new bpp::RHomogeneousTreeLikelihood(*tree, *sites, models[sharedmodelindex], ratedistribution, true, verbose, true);
            }
        } else {
            if (infertopology) {
                throw std::invalid_argument("Cannot infer topology with a mixed model");
            } else {
                oldtreelikelihood = new bpp::RHomogeneousMixedTreeLikelihood(*tree, *sites, models[sharedmodelindex], ratedistribution, true, verbose, true);
            }
        }
        oldtreelikelihood->initialize();
    }

    // check that we can compute initial log likelihood, throws an error if this is 0
    LogLikelihood();
}

// destructor
bppextensions::BppTreeLikelihood::~BppTreeLikelihood()
{
    if (oldtreelikelihood) delete oldtreelikelihood;
    if (phylolikelihood) delete phylolikelihood;
    if (substitutionprocesscollection) {
        delete substitutionprocesscollection;
    } else {
        for (std::map<size_t, bpp::SubstitutionModel*>::iterator itr = models.begin(); itr != models.end(); itr++) {
            if (itr->second) delete itr->second;
        }
        if (ratedistribution) delete ratedistribution;
        if (paramtree) delete paramtree;
    }
    if (sites) delete sites;
    if (ntalphabet) delete ntalphabet;
    if (alphabet) delete alphabet;
    if (gcode) delete gcode;
    if (treeReaderWriter) delete treeReaderWriter;
    if (sequenceevolution) delete sequenceevolution;
    if (tree) delete tree;
}

long bppextensions::BppTreeLikelihood::NSeqs()
{
    return (long) sites->getNumberOfSequences();
}

long bppextensions::BppTreeLikelihood::NSites()
{
    return (long) sites->getNumberOfSites();
}

void bppextensions::BppTreeLikelihood::NewickTree(std::string fname)
{
    treeReaderWriter->write(*tree, fname);
}

double bppextensions::BppTreeLikelihood::LogLikelihood()
{
    double logL;
    if (oldlikmethod) {
        logL = oldtreelikelihood->getValue();
    } else {
        logL = phylolikelihood->getValue();
    }
    if (std::isinf(logL)) {
        throw std::runtime_error("Tree likelihood is zero. Perhaps you have branch lengths of zero, stop codons, or a very bad tree topology?");
    }
    if (std::isnan(logL)) {
        throw std::runtime_error("Tree likelihood is nan. Some problem in likelihood computation. Or maybe you have branch lengths of zero, stop codons, or a very bad tree topology?");
    }
    return -logL;
}

void bppextensions::BppTreeLikelihood::OptimizeLikelihood()
{
    if (oldlikmethod) {
        oldtreelikelihood = dynamic_cast<bpp::DiscreteRatesAcrossSitesTreeLikelihood*>(bpp::PhylogeneticsApplicationTools::optimizeParameters(oldtreelikelihood, oldtreelikelihood->getParameters(), optimizationparams, "", true, verbose, 1));
        tree = new bpp::TreeTemplate<bpp::Node>(oldtreelikelihood->getTree());
    } else {
        phylolikelihood = bpp::PhylogeneticsApplicationTools::optimizeParameters(phylolikelihood, phylolikelihood->getParameters(), optimizationparams, "", true, verbose, 1);
        tree = new bpp::TreeTemplate<bpp::Node>(substitutionprocesscollection->getTree(sharedtreeindex).getTree());
    }
}

std::map<std::string, double> bppextensions::BppTreeLikelihood::ModelParams()
{
    std::map<std::string, double> modelparametersmap;
    if (oldlikmethod) {
        if (nrates > 1) {
            modelparametersmap["Gamma.alpha"] = oldtreelikelihood->getParameters().getParameterValue("Gamma.alpha");
        }
        models[sharedmodelindex]->matchParametersValues(oldtreelikelihood->getParameters());
        std::vector<std::string> modelparameters = models[sharedmodelindex]->getIndependentParameters().getParameterNames();
        for (std::vector<std::string>::size_type i = 0; i != modelparameters.size(); ++i) {
            std::string parametername = modelparameters[i];
            modelparametersmap[parametername] = models[sharedmodelindex]->getIndependentParameters().getParameterValue(parametername);
        }
    } else {
        if (nrates > 1) {
            modelparametersmap["Gamma.alpha"] = phylolikelihood->getParameters().getParameterValue("Gamma.alpha_" + patch::to_string(sharedrateindex));
        }
        std::vector<size_t> modelnumbers = substitutionprocesscollection->getModelNumbers();
        for (size_t i = 0; i < modelnumbers.size(); i++) {
            const bpp::SubstitutionModel& imodel = substitutionprocesscollection->getModel(modelnumbers[i]);
            std::vector<std::string> imodelparameters = imodel.getIndependentParameters().getParameterNames();
            for (std::vector<std::string>::size_type j = 0; j != imodelparameters.size(); ++j) {
                std::string parametername = imodelparameters[j];
                if (constrainedparams.find(parametername) != constrainedparams.end()) {
                    if (modelparametersmap.find(parametername) == modelparametersmap.end()) {
                        modelparametersmap[parametername] = imodel.getIndependentParameters().getParameterValue(parametername);
                    } else {
                        if (modelparametersmap[parametername] != imodel.getIndependentParameters().getParameterValue(parametername)) {
                            throw std::runtime_error("Constrained parameter " + parametername + " does not have the same value for all sites");
                        }
                    }
                } else {
                    modelparametersmap["site" + patch::to_string(modelnumbers[i]) + "." + parametername] = imodel.getIndependentParameters().getParameterValue(parametername);
                }
            }
        }
    }
    return modelparametersmap;
}

std::string bppextensions::BppTreeLikelihood::OptimizationIgnoredParameters()
{
    return optimizationparams["optimization.ignore_parameters"];
}

std::map<std::string, double> bppextensions::BppTreeLikelihood::StationaryState(long isite)
{
    if ((isite < 1) || (isite > NSites())) {
        throw std::runtime_error("Invalid site number, must be >= 1 and <= NSites");
    }
    if (oldlikmethod) {
        isite = sharedmodelindex;
    }
    map<std::string, double> stationarystate;
    for (size_t icodon = 0; icodon < models[isite]->getNumberOfStates(); icodon++) { 
        std::string codon = models[isite]->getAlphabetStateAsChar(icodon);
        stationarystate[codon] = models[isite]->freq((int) icodon);
    }
    return stationarystate;
}
             
std::map<std::string, double> bppextensions::BppTreeLikelihood::GetPreferences(long isite)
{
    if ((isite < 1) || (isite > NSites())) {
        throw std::runtime_error("Invalid site number, must be >= 1 and <= NSites");
    }
    if (models.find(isite) == models.end()) {
        throw std::runtime_error("There is not a site with key " + patch::to_string(isite) + ". Are you sure you are using ExpCM with just one model?");
    }
    bppextensions::ExperimentallyInformedCodonModel *model = dynamic_cast<bppextensions::ExperimentallyInformedCodonModel*>(models[isite]);
    if (! model) {
        throw std::runtime_error("You did not use an ExpCM model");
    }
    std::map<std::string, double> aaprefs;
    std::map<std::string, double> codonprefs = model->getPreferences();
    double sum = 0.0;
    for (std::map<std::string, double>::iterator itr = codonprefs.begin(); itr != codonprefs.end(); itr++) 
    {
        if (gcode->isStop(itr->first)) {
            if (itr->second > 1.0e-8) {
                throw std::runtime_error("Stop codon " + itr->first + " has substantially nonzero value of " + patch::to_string(itr->second));
            }
        } else {
            std::string aa = gcode->translate(itr->first);
            if (aaprefs.find(aa) == aaprefs.end()) {
                aaprefs[aa] = itr->second;
                sum += itr->second;
            } else if (std::fabs(aaprefs[aa] - itr->second) > 1.0e-8) {
                throw std::runtime_error("Inconsistent preferences for " + aa + ": " + patch::to_string(aaprefs[aa]) + " and " + patch::to_string(itr->second) + "\n");
            }
        }
    }
    if ((sum <= 0) || (sum > 1)) {
        throw std::runtime_error("sum of aaprefs is invalid value of " + patch::to_string(sum));
    }
    for (std::map<std::string, double>::iterator itr = aaprefs.begin(); itr != aaprefs.end(); itr++) {
        itr->second /= sum;
    }
    return aaprefs;
}

void bppextensions::BppTreeLikelihood::SetPreferences(std::map<std::string, double> aaprefs, long isite) {
    if ((isite < 1) || (isite > NSites())) {
        throw std::runtime_error("Invalid site number, must be >= 1 and <= NSites");
    }
    if (models.find(isite) == models.end()) {
        throw std::runtime_error("There is not a site with key " + patch::to_string(isite) + ". Are you sure you are using ExpCM with just one model?");
    }
    bppextensions::ExperimentallyInformedCodonModel *model = dynamic_cast<bppextensions::ExperimentallyInformedCodonModel*>(models[isite]);
    if (! model) {
        throw std::runtime_error("You did not use an ExpCM model");
    }
    double sum = 0.0;
    for (std::map<std::string, double>::iterator itr = aaprefs.begin(); itr != aaprefs.end(); itr++) {
        sum += itr->second;
    }
    if (std::fabs(sum - 1.0) > 1.e-4) {
        throw std::runtime_error("Sum of amino-acid preferences is not close to one: " + patch::to_string(sum));
    }
    std::map<int, double> initcodonprefs;
    bpp::FullCodonFrequenciesSet *codonprefs = new bpp::FullCodonFrequenciesSet(gcode);
    for (size_t icodon = 0; icodon < codonprefs->getNumberOfFrequencies(); icodon++) {
        if (gcode->isStop(icodon)) {
            initcodonprefs[icodon] = 0.0;
        } else {
            std::string codon = codonprefs->getAlphabet()->intToChar((int) icodon);
            std::string aa = gcode->translate(codon);
            if (aaprefs.find(aa) == aaprefs.end()) {
                throw std::runtime_error("Failed to find preference for amino acid " + aa + "\n");
            }
            initcodonprefs[icodon] = aaprefs[aa];
        }
    }
    codonprefs->setFrequenciesFromAlphabetStatesFrequencies(initcodonprefs);
    codonprefs->setNamespace(model->getPreferencesNamespace());
    bpp::ParameterList preflist = codonprefs->getParameters();
    std::vector<std::string> prefnames = preflist.getParameterNames();
    bpp::ParameterList *preflist_renamed = new bpp::ParameterList();
    for (size_t i = 0; i < prefnames.size(); i++) {
        bpp::Parameter iparam = preflist.getParameter(prefnames[i]);
        // add site number suffix to match phylolikelihood names
        iparam.setName(iparam.getName() + "_" + patch::to_string(isite));
        preflist_renamed->addParameter(iparam);
    }
    std::vector<std::string> prefnames_renamed = preflist_renamed->getParameterNames();
    for (size_t i = 0; i < prefnames_renamed.size(); i++) {
        if (! phylolikelihood->getSubstitutionModelParameters().hasParameter(prefnames_renamed[i])) {
            throw std::runtime_error("Can't find parameter " + prefnames_renamed[i] + ". Did you use an ExpCM with prefsasparams?");
        }
    }
    phylolikelihood->setParameters(*preflist_renamed);
    delete codonprefs;
    delete preflist_renamed;
}

