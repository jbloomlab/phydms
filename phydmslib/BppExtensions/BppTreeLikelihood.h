/* BppTreeLikelihood.h

   Defines a C++ class that interfaces with Bio++
   to define a tree, model, and likelihood calculation object.
   This class comprehensively performs all of the operations
   needed to interface with Bio++.

   Written by Jesse Bloom
*/

#include <string>
#include <vector>
#include <Bpp/Seq/Alphabet/NucleicAlphabet.h>
#include <Bpp/Seq/Alphabet/DNA.h>
#include <Bpp/Seq/Alphabet/Alphabet.h>
#include <Bpp/Seq/Alphabet/CodonAlphabet.h>
#include <Bpp/Seq/GeneticCode/GeneticCode.h>
#include <Bpp/Seq/GeneticCode/StandardGeneticCode.h>
#include <Bpp/Seq/Container/VectorSiteContainer.h>
#include <Bpp/Seq/Container/SiteContainerTools.h>
#include <Bpp/Seq/Sequence.h>
#include <Bpp/Phyl/App/PhylogeneticsApplicationTools.h>
#include <Bpp/Phyl/Io/IoTree.h>
#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Phyl/Model/RateDistribution/ConstantRateDistribution.h>
#include <Bpp/Phyl/Model/RateDistribution/GammaDiscreteRateDistribution.h>
#include <Bpp/Phyl/Likelihood/DiscreteRatesAcrossSitesTreeLikelihood.h>
#include <Bpp/Phyl/Likelihood/NNIHomogeneousTreeLikelihood.h>
#include <Bpp/Phyl/Likelihood/RHomogeneousTreeLikelihood.h>
#include <Bpp/Phyl/Likelihood/RHomogeneousMixedTreeLikelihood.h>
#include <Bpp/Phyl/NewLikelihood/SequenceEvolution.h>
#include <Bpp/Phyl/NewLikelihood/PhyloLikelihoods/PhyloLikelihood.h>
#include <Bpp/Phyl/NewLikelihood/PartitionSequenceEvolution.h>
#include <Bpp/Phyl/NewLikelihood/PhyloLikelihoods/PartitionProcessPhyloLikelihood.h>

namespace bppextensions {

    /**
     *@brief Interface between bpp and the Python cython wrapper, defined as class that does all likelihood computations.
     *
     * You initialize an object of this class with the sequences, initial tree, substitution model, and various
     * other relevant options. The object then can be manipulated via the interface defined here for all necessary
     * operations.
     */
    class BppTreeLikelihood {

        public:
            /**
             *@brief Constructor
             *
             *@param seqnames A vector of strings giving the sequence names
             *
             *@param seqs A vector of strings giving the sequences (should all be aligned, no stop codons)
             *
             *@param treefile Name of an existing file giving Newick tree with names matching those in seqnames
             *
             *@param modelstring A string defining the model. See the Python wrapper for list of valid values.
             *
             * @param infertopology Do we infer the tree topology by maximum likelihood? 
             *
             * @param preferences Site-specific amino-acid preferences, keyed by integer site (1, 2, ... numbering), then maps keyed by codon and values preference. Value does not matter if modelstring is not "ExpCM".
             *
             * @param fixedmodelparams A map giving model parameters that should be fixed. Key is string name of model parameter, value is what it is fixed to. You will get an error if you specify a parameter that does not exist with that name.
             *
             * @param initializemodelparams A map giving model parameters that should be initialized but not fixed. Key is string name of model parameter, value is what it is initialized to. You will get an error if you specify a parameter that does not exist with that name, or if an entry is duplicated in fixedmodelparams.
             *
             * @param oldlikelihoodmethod Do we use the old Bpp likelihood method rather then the NewLikelihood? Only can be used for non-partitioned data.
             * 
             * @param fixbrlen Do we fix the branch lengths?
             *
             * @param addrateparameter Add a parameter that scales all rates. You can only use this if fixing branch lengths (fixbrlen =! 0)
             *
             * @param prefsasparams Do we define site-specific amino-acid preferences as parameters for ExpCM?
             *
             * @param recursion Recursion method used for likelihood. Can be "S" for simple or "D" for double.
             *
             * @param useLog Do we compute the likelihoods using their logs? Only meaningful when not using oldlikelihoodmethod.
             *
             * @param ngammarates Number of gamma-distributed rate categories. Set to 1 for a single constant rate, otherwise an integer > 1 for gamma distributed rates.
             *
             * @param ncats Number of categories in beta distribution for YNGKP_M7 and YNGPK_M8
             */
            BppTreeLikelihood(std::vector<std::string> seqnames, std::vector<std::string> seqs, std::string treefile, std::string modelstring, int infertopology, std::map<int, std::map<std::string, double> > preferences, std::map<std::string, double> fixedmodelparams, std::map<std::string, double> initializemodelparams, int oldlikelihoodmethod, int fixbrlen, int addrateparameter, int prefsasparams, char recursion, int useLog, int ngammarates, int ncats);

            /**
             *@brief Destructor
             */
            ~BppTreeLikelihood();

            /**
             *@return returns number of sequences
             */
            long NSeqs();

            /**
             *@return returns number of sites
             */
            long NSites();

            /**
             *@brief Writes the current Newick tree to a file

             *@param fname Name of file to which we write the tree
             */
            void NewickTree(std::string fname); 
            
            /**
             *@return Returns current log likelihood
             */
            double LogLikelihood(); 

            /**
             *@brief Optimizes the object by maximum likelihood. May optimize topology, branches, and model parameters depending on how object was initialized.
             */
            void OptimizeLikelihood(); 

            /**
             *@brief Gets current values of model parameters.
             *
             *@return map keyed by parameter names, values current values
             *
             */
            std::map<std::string, double> ModelParams(); 

            /**
             *@brief Returns current stationary state of substitition model for a site
             *
             *@param isite the site for which we get the stationary state in 1, 2, ..., NSites() numbering
             *
             *@return map keyed by codons, values are stationary state
             */
            std::map<std::string, double> StationaryState(long isite);

            /**
             *@brief returns the site-specific preferences. Will work only if object was constructed for an ExpCM.
             *
             * @param isite is the site for which we get the stationary state in 1, 2, ..., NSites() numbering
             *
             *@return map keyed by amino-acids, values as preferences.
             */
             std::map<std::string, double> GetPreferences(long isite);

            /**
             *@brief Sets new values for the site-specific preferences for a site. Will only work if object was constructed for an ExpCM.
             *
             * @param aaprefs is a map keyed by amino acids with values as preferences, we set to these values.
             *
             * @param isite is the site for which we set the preferences in 1, 2, ..., NSites() numbering.
             */
             void SetPreferences(std::map<std::string, double> aaprefs, long isite);

            /**
             *@brief Returns parameters currently being ignored for optimization
             *
             *@return a comma-delimited string of the parameters being ignore, which may include wildcards indicated by *
             *
             */
            std::string OptimizationIgnoredParameters();


        private:
            bool verbose; // verbosity of Bpp functions
            bool oldlikmethod;
            std::map<std::string, std::string> optimizationparams; // specifies parameters for optimization
            bpp::NucleicAlphabet *ntalphabet;
            bpp::CodonAlphabet *alphabet;
            bpp::GeneticCode *gcode;
            bpp::VectorSiteContainer *sites;
            bpp::Tree *tree;
            bpp::ParametrizableTree *paramtree;
            bpp::Newick *treeReaderWriter;
            std::map<size_t, bpp::SubstitutionModel*> models;
            bpp::DiscreteDistribution *ratedistribution;
            bpp::SubstitutionProcessCollection *substitutionprocesscollection;
            bpp::SequenceEvolution *sequenceevolution;
            bpp::PhyloLikelihood *phylolikelihood;
            bpp::DiscreteRatesAcrossSitesTreeLikelihood *oldtreelikelihood;
            size_t sharedmodelindex;
            size_t sharedtreeindex;
            size_t sharedrateindex;
            std::map<std::string, std::string> constrainedparams; // keyed by param, value is what it is constrained to
            int nrates;
    };
};

