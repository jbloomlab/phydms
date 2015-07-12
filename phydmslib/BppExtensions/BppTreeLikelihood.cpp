/* BppTreeLikelihood.cpp
   Implements C++ class BppTreeLikelihood
   Written by Jesse Bloom
*/

#include <iostream>
#include <fstream>
#include <exception>
#include <string>
#include <sstream>
#include <Bpp/Phyl/Model/FrequenciesSet/FrequenciesSet.h>
#include <Bpp/Phyl/Io/BppOFrequenciesSetFormat.h>
#include <Bpp/Phyl/Model/Codon/YN98.h>
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
bppextensions::BppTreeLikelihood::BppTreeLikelihood(std::vector<std::string> seqnames, std::vector<std::string> seqs, std::string treefile, std::string modelstring, int infertopology, std::map<int, std::map<std::string, double> > preferences, std::map<std::string, double> fixedmodelparams, int oldlikelihoodmethod, int fixbrlen, int addrateparameter, char recursion)
{

    // setup some parameters / options
    oldlikmethod = oldlikelihoodmethod != 0;
    verbose = false;
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
    optimizationparams["optimization.ignore_parameters"] = "";

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
    size_t sharedrateindex = 0;
    sharedmodelindex = 0;
    sequenceevolution = 0;
    tree = 0;

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
            if (optimizationparams["optimization.ignore_parameters"].empty()) {
                optimizationparams["optimization.ignore_parameters"] = "*theta*"; 
            } else {
                optimizationparams["optimization.ignore_parameters"] = optimizationparams["optimization.ignore_parameters"] + "," + "*theta*";
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
        else if (modelstring.substr(6, 2) == "M7") {
            models[sharedmodelindex] = dynamic_cast<bpp::SubstitutionModel*>(new bpp::YNGKP_M7(gcode, codonFreqs.release(), 3));
        }
        else if (modelstring.substr(6, 2) == "M8") {
            models[sharedmodelindex] = dynamic_cast<bpp::SubstitutionModel*>(new bpp::YNGKP_M8(gcode, codonFreqs.release(), 3));
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
        std::vector<double> init_rprefs;
        std::string codon;
        for (long isite = 1; isite <= nsites; isite++) {
            init_rprefs.clear();
            for (long icodon = 0; icodon < alphabet->getNumberOfTypes() - 1; icodon++) { 
                codon = alphabet->intToChar(icodon);
                init_rprefs.push_back(preferences[isite][codon]);
            }
            init_rprefs[alphabet->getNumberOfTypes() - 1] = 0.0; // this is the ambiguous character code
            bpp::FullCodonFrequenciesSet *rprefs = new bpp::FullCodonFrequenciesSet(gcode, init_rprefs);
            models[isite] = dynamic_cast<bpp::SubstitutionModel*>(new bppextensions::ExperimentallyInformedCodonModel(gcode, rprefs, "ExpCM."));
            delete rprefs;
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

    // fix any model parameters in fixedmodelparams
    if (! fixedmodelparams.empty()) {
        if (oldlikmethod) {
            throw std::runtime_error("fixedmodelparams is not guaranteed to work if using oldlikmethod");
        }
        for (std::map<std::string, double>::iterator itr = fixedmodelparams.begin(); itr != fixedmodelparams.end(); itr++) {
            for (std::map<size_t, bpp::SubstitutionModel*>::iterator imodel_itr = models.begin(); imodel_itr != models.end(); ++imodel_itr) {
                if (imodel_itr->second->hasParameter(itr->first)) {
                    imodel_itr->second->setParameterValue(itr->first, itr->second);
                } else {
                    throw std::runtime_error("Cannot find parameter in fixedmodelparams: " + itr->first);
                }
            }
            if (optimizationparams["optimization.ignore_parameters"].empty()) {
                optimizationparams["optimization.ignore_parameters"] = "*" + itr->first + "*";
            } else {
                optimizationparams["optimization.ignore_parameters"] = optimizationparams["optimization.ignore_parameters"] + "," + "*" + itr->first + "*";
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
    ratedistribution = new bpp::ConstantRateDistribution(); // just a single rate

    // compute likelihoods
    oldtreelikelihood = 0;
    phylolikelihood = 0;
    if (! oldlikmethod) {
        if (infertopology) {
            throw std::runtime_error("Cannot use infertopology with new likelihood method");
        }
        // combine all of these into collection, similar to what is done by bpp::PhylogeneticsApplicationTools::getSubstitutionProcessCollection
        substitutionprocesscollection = new bpp::SubstitutionProcessCollection();
        substitutionprocesscollection->addTree(new bpp::ParametrizableTree(*tree), sharedtreeindex); 
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
        char compression = 'R'; // appears to be Bpp default, not sure what it means
        bpp::PartitionSequenceEvolution* pse = dynamic_cast<bpp::PartitionSequenceEvolution*>(sequenceevolution);
        if (pse == NULL) {
            throw std::runtime_error("Failed to cast to PartitionSequenceEvolution");
        }
        phylolikelihood = new bpp::PartitionPhyloLikelihood(*sites, *pse, recursion, 0, 0, verbose, compression=='R');

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
    for (std::map<size_t, bpp::SubstitutionModel*>::iterator itr = models.begin(); itr != models.end(); itr++) {
        if (itr->second) delete itr->second;
    }
    if (sites) delete sites;
    if (ntalphabet) delete ntalphabet;
    if (alphabet) delete alphabet;
    if (gcode) delete gcode;
    if (treeReaderWriter) delete treeReaderWriter;
    if (ratedistribution) delete ratedistribution;
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
        models[sharedmodelindex]->matchParametersValues(oldtreelikelihood->getParameters());
        std::vector<std::string> modelparameters = models[sharedmodelindex]->getIndependentParameters().getParameterNames();
        for (std::vector<std::string>::size_type i = 0; i != modelparameters.size(); ++i) {
            std::string parametername = modelparameters[i];
            modelparametersmap[parametername] = models[sharedmodelindex]->getIndependentParameters().getParameterValue(parametername);
        }
    } else {
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
        throw std::runtime_error("Invalid site number, must be > 1 and < NSites");
    }
    if (oldlikmethod) {
        isite = sharedmodelindex;
    }
    map<std::string, double> stationarystate;
    for (long icodon = 0; icodon < alphabet->getNumberOfTypes() - 1; icodon++) { 
        std::string codon = alphabet->intToChar(icodon);
        stationarystate[codon] = models[isite]->freq((int) icodon);
    }
    return stationarystate;
}
