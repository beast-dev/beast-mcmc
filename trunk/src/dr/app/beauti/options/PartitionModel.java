package dr.app.beauti.options;

import dr.app.beauti.priorsPanel.PriorType;
import dr.evolution.datatype.DataType;
import dr.evomodel.tree.RateStatistic;
import dr.evomodelxml.BirthDeathModelParser;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 */
public class PartitionModel extends ModelOptions {

    public static final String[] GTR_RATE_NAMES = {"ac", "ag", "at", "cg", "gt"};
    static final String[] GTR_TRANSITIONS = {"A-C", "A-G", "A-T", "C-G", "G-T"};

    public PartitionModel(BeautiOptions options, DataPartition partition) {
        this(options, partition.getName(), partition.getAlignment().getDataType());
    }

    /**
     * A copy constructor
     *
     * @param options the beauti options
     * @param name    the name of the new model
     * @param source  the source model
     */
    public PartitionModel(BeautiOptions options, String name, PartitionModel source) {
        this(options, name, source.dataType);

        nucSubstitutionModel = source.nucSubstitutionModel;
        aaSubstitutionModel = source.aaSubstitutionModel;
        binarySubstitutionModel = source.binarySubstitutionModel;

        frequencyPolicy = source.frequencyPolicy;
        gammaHetero = source.gammaHetero;
        gammaCategories = source.gammaCategories;
        invarHetero = source.invarHetero;
        codonHeteroPattern = source.codonHeteroPattern;
        unlinkedSubstitutionModel = source.unlinkedSubstitutionModel;
        unlinkedHeterogeneityModel = source.unlinkedHeterogeneityModel;
        unlinkedFrequencyModel = source.unlinkedFrequencyModel;
    }

    public PartitionModel(BeautiOptions options, String name, DataType dataType) {

        this.options = options;
        this.name = name;
        this.dataType = dataType;

        double substWeights = 1.0;

        //Substitution model parameters
        createParameter("frequencies", "base frequencies", UNITY_SCALE, 0.25, 0.0, 1.0);
        createParameter("CP1.frequencies", "base frequencies for codon position 1", UNITY_SCALE, 0.25, 0.0, 1.0);
        createParameter("CP2.frequencies", "base frequencies for codon position 2", UNITY_SCALE, 0.25, 0.0, 1.0);
        createParameter("CP1+2.frequencies", "base frequencies for codon positions 1 & 2", UNITY_SCALE, 0.25, 0.0, 1.0);
        createParameter("CP3.frequencies", "base frequencies for codon position 3", UNITY_SCALE, 0.25, 0.0, 1.0);

        createScaleParameter("kappa", "HKY transition-transversion parameter", SUBSTITUTION_PARAMETER_SCALE, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        createScaleParameter("CP1.kappa", "HKY transition-transversion parameter for codon position 1", SUBSTITUTION_PARAMETER_SCALE, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        createScaleParameter("CP2.kappa", "HKY transition-transversion parameter for codon position 2", SUBSTITUTION_PARAMETER_SCALE, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        createScaleParameter("CP1+2.kappa", "HKY transition-transversion parameter for codon positions 1 & 2", SUBSTITUTION_PARAMETER_SCALE, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        createScaleParameter("CP3.kappa", "HKY transition-transversion parameter for codon position 3", SUBSTITUTION_PARAMETER_SCALE, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);

//        createParameter("frequencies", "GTR base frequencies", UNITY_SCALE, 0.25, 0.0, 1.0);
//        createParameter("CP1.frequencies", "GTR base frequencies for codon position 1", UNITY_SCALE, 0.25, 0.0, 1.0);
//        createParameter("CP2.frequencies", "GTR base frequencies for codon position 2", UNITY_SCALE, 0.25, 0.0, 1.0);
//        createParameter("CP1+2.frequencies", "GTR base frequencies for codon positions 1 & 2", UNITY_SCALE, 0.25, 0.0, 1.0);
//        createParameter("CP3.frequencies", "GTR base frequencies for codon position 3", UNITY_SCALE, 0.25, 0.0, 1.0);

        // create the relative rate parameters for the GTR rate matrix
        for (int j = 0; j < 5; j++) {
            createScaleParameter(GTR_RATE_NAMES[j], "GTR " + GTR_TRANSITIONS[j] + " substitution parameter",
                    SUBSTITUTION_PARAMETER_SCALE, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
            for (int i = 1; i <= 3; i++) {

                createScaleParameter(
                        "CP" + i + "." + GTR_RATE_NAMES[j],
                        "GTR " + GTR_TRANSITIONS[j] + " substitution parameter for codon position " + i,
                        SUBSTITUTION_PARAMETER_SCALE, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
            }
            createScaleParameter("CP1+2." + GTR_RATE_NAMES[j],
                    "GTR " + GTR_TRANSITIONS[j] + " substitution parameter for codon positions 1 & 2",
                    SUBSTITUTION_PARAMETER_SCALE, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        }

//        createParameter("frequencies", "Binary Simple frequencies", UNITY_SCALE, 0.5, 0.0, 1.0);
//
//        createParameter("frequencies", "Binary Covarion frequencies of the visible states", UNITY_SCALE, 0.5, 0.0, 1.0);
        createParameter("hfrequencies", "Binary Covarion frequencies of the hidden rates", UNITY_SCALE, 0.5, 0.0, 1.0);
        createParameter("bcov.alpha", "Binary Covarion rate of evolution in slow mode", UNITY_SCALE, 0.5, 0.0, 1.0);
        createParameter("bcov.s", "Binary Covarion rate of flipping between slow and fast modes", SUBSTITUTION_PARAMETER_SCALE, 0.5, 0.0, 100.0);

        createParameter("alpha", "gamma shape parameter", SUBSTITUTION_PARAMETER_SCALE, 0.5, 0.0, 1000.0);
        createParameter("CP1.alpha", "gamma shape parameter for codon position 1", SUBSTITUTION_PARAMETER_SCALE, 0.5, 0.0, 1000.0);
        createParameter("CP2.alpha", "gamma shape parameter for codon position 2", SUBSTITUTION_PARAMETER_SCALE, 0.5, 0.0, 1000.0);
        createParameter("CP1+2.alpha", "gamma shape parameter for codon positions 1 & 2", SUBSTITUTION_PARAMETER_SCALE, 0.5, 0.0, 1000.0);
        createParameter("CP3.alpha", "gamma shape parameter for codon position 3", SUBSTITUTION_PARAMETER_SCALE, 0.5, 0.0, 1000.0);

        createParameter("pInv", "proportion of invariant sites parameter", NONE, 0.5, 0.0, 1.0);
        createParameter("CP1.pInv", "proportion of invariant sites parameter for codon position 1", NONE, 0.5, 0.0, 1.0);
        createParameter("CP2.pInv", "proportion of invariant sites parameter for codon position 2", NONE, 0.5, 0.0, 1.0);
        createParameter("CP1+2.pInv", "proportion of invariant sites parameter for codon positions 1 & 2", NONE, 0.5, 0.0, 1.0);
        createParameter("CP3.pInv", "proportion of invariant sites parameter for codon position 3", NONE, 0.5, 0.0, 1.0);

        createParameter("mu", "relative rate parameter", SUBSTITUTION_PARAMETER_SCALE, 1.0, 0.0, Double.POSITIVE_INFINITY);
        createParameter("CP1.mu", "relative rate parameter for codon position 1", SUBSTITUTION_PARAMETER_SCALE, 1.0, 0.0, Double.POSITIVE_INFINITY);
        createParameter("CP2.mu", "relative rate parameter for codon position 2", SUBSTITUTION_PARAMETER_SCALE, 1.0, 0.0, Double.POSITIVE_INFINITY);
        createParameter("CP1+2.mu", "relative rate parameter for codon positions 1 & 2", SUBSTITUTION_PARAMETER_SCALE, 1.0, 0.0, Double.POSITIVE_INFINITY);
        createParameter("CP3.mu", "relative rate parameter for codon position 3", SUBSTITUTION_PARAMETER_SCALE, 1.0, 0.0, Double.POSITIVE_INFINITY);

        createScaleOperator("kappa", substWeights);
        createScaleOperator("CP1.kappa", substWeights);
        createScaleOperator("CP2.kappa", substWeights);
        createScaleOperator("CP1+2.kappa", substWeights);
        createScaleOperator("CP3.kappa", substWeights);

        createOperator("frequencies", OperatorType.DELTA_EXCHANGE, 0.01, substWeights);
        createOperator("CP1.frequencies", OperatorType.DELTA_EXCHANGE, 0.01, substWeights);
        createOperator("CP2.frequencies", OperatorType.DELTA_EXCHANGE, 0.01, substWeights);
        createOperator("CP1+2.frequencies", OperatorType.DELTA_EXCHANGE, 0.01, substWeights);
        createOperator("CP3.frequencies", OperatorType.DELTA_EXCHANGE, 0.01, substWeights);

        for (String rateName : GTR_RATE_NAMES) {
            createScaleOperator(rateName, substWeights);
            for (int j = 1; j <= 3; j++) {
                createScaleOperator("CP" + j + "." + rateName, substWeights);
            }
            createScaleOperator("CP1+2." + rateName, substWeights);
        }

        createScaleOperator("alpha", substWeights);
        for (int i = 1; i <= 3; i++) {
            createScaleOperator("CP" + i + ".alpha", substWeights);
        }
        createScaleOperator("CP1+2.alpha", substWeights);

        createScaleOperator("pInv", substWeights);
        for (int i = 1; i <= 3; i++) {
            createScaleOperator("CP" + i + ".pInv", substWeights);
        }
        createScaleOperator("CP1+2.pInv", substWeights);

        createScaleOperator("bcov.alpha", substWeights);
        createScaleOperator("bcov.s", substWeights);
        //createOperator("frequencies", OperatorType.DELTA_EXCHANGE, 0.01, substWeights);
        createOperator("hfrequencies", OperatorType.DELTA_EXCHANGE, 0.01, substWeights);
        
    }    
    
    // use override method getOperator(String name) in PartitionModel containing prefix
    public void initParametersAndOperatorsForEachPartitionModel () {
        double demoWeights = 3.0;
        double branchWeights = 30.0;
        double treeWeights = 15.0;
        double rateWeights = 3.0;

        createParameter("tree", "The tree");
        createParameter("treeModel.internalNodeHeights", "internal node heights of the tree (except the root)");
        createParameter("treeModel.allInternalNodeHeights", "internal node heights of the tree");
        createParameter("treeModel.rootHeight", "root height of the tree", true, 1.0, 0.0, Double.POSITIVE_INFINITY);

        // A vector of relative rates across all partitions...
        createParameter("allMus", "All the relative rates");

        createParameter("clock.rate", "substitution rate", SUBSTITUTION_RATE_SCALE, 1.0, 0.0, Double.POSITIVE_INFINITY);
        createParameter(ClockType.UCED_MEAN, "uncorrelated exponential relaxed clock mean", SUBSTITUTION_RATE_SCALE, 1.0, 0.0, Double.POSITIVE_INFINITY);
        createParameter(ClockType.UCLD_MEAN, "uncorrelated lognormal relaxed clock mean", SUBSTITUTION_RATE_SCALE, 1.0, 0.0, Double.POSITIVE_INFINITY);
        createParameter(ClockType.UCLD_STDEV, "uncorrelated lognormal relaxed clock stdev", LOG_STDEV_SCALE, 0.1, 0.0, Double.POSITIVE_INFINITY);
        createParameter("branchRates.categories", "relaxed clock branch rate categories");
        createParameter(ClockType.LOCAL_CLOCK + "." + "rates", "random local clock rates", SUBSTITUTION_RATE_SCALE, 1.0, 0.0, Double.POSITIVE_INFINITY);
        createParameter(ClockType.LOCAL_CLOCK + "." + "changes", "random local clock rate change indicator");

        {
            final Parameter p = createParameter("treeModel.rootRate", "autocorrelated lognormal relaxed clock root rate", ROOT_RATE_SCALE, 1.0, 0.0, Double.POSITIVE_INFINITY);
            p.priorType = PriorType.GAMMA_PRIOR;
            p.gammaAlpha = 1;
            p.gammaBeta = 0.0001;
        }
        createParameter("treeModel.nodeRates", "autocorrelated lognormal relaxed clock non-root rates", SUBSTITUTION_RATE_SCALE, 1.0, 0.0, Double.POSITIVE_INFINITY);
        createParameter("treeModel.allRates", "autocorrelated lognormal relaxed clock all rates", SUBSTITUTION_RATE_SCALE, 1.0, 0.0, Double.POSITIVE_INFINITY);
        {
            final Parameter p = createParameter("branchRates.var", "autocorrelated lognormal relaxed clock rate variance ", LOG_VAR_SCALE, 0.1, 0.0, Double.POSITIVE_INFINITY);
            p.priorType = PriorType.GAMMA_PRIOR;
            p.gammaAlpha = 1;
            p.gammaBeta = 0.0001;
        }

        createScaleParameter("constant.popSize", "coalescent population size parameter", TIME_SCALE, 1.0, 0.0, Double.POSITIVE_INFINITY);

        createScaleParameter("exponential.popSize", "coalescent population size parameter", TIME_SCALE, 1.0, 0.0, Double.POSITIVE_INFINITY);
        createParameter("exponential.growthRate", "coalescent growth rate parameter", GROWTH_RATE_SCALE, 0.0, Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY);
        createParameter("exponential.doublingTime", "coalescent doubling time parameter", TIME_SCALE, 0.5, 0.0, Double.POSITIVE_INFINITY);
        createScaleParameter("logistic.popSize", "coalescent population size parameter", TIME_SCALE, 1.0, 0.0, Double.POSITIVE_INFINITY);
        createParameter("logistic.growthRate", "coalescent logistic growth rate parameter", GROWTH_RATE_SCALE, 0.001, 0.0, Double.POSITIVE_INFINITY);
        createParameter("logistic.doublingTime", "coalescent doubling time parameter", TIME_SCALE, 0.5, 0.0, Double.POSITIVE_INFINITY);
        createParameter("logistic.t50", "logistic shape parameter", T50_SCALE, 0.1, 0.0, Double.POSITIVE_INFINITY);
        createScaleParameter("expansion.popSize", "coalescent population size parameter", TIME_SCALE, 1.0, 0.0, Double.POSITIVE_INFINITY);
        createParameter("expansion.growthRate", "coalescent logistic growth rate parameter", GROWTH_RATE_SCALE, 0.001, 0.0, Double.POSITIVE_INFINITY);
        createParameter("expansion.doublingTime", "coalescent doubling time parameter", TIME_SCALE, 0.5, 0.0, Double.POSITIVE_INFINITY);
        createParameter("expansion.ancestralProportion", "ancestral population proportion", NONE, 0.1, 0.0, 1.0);
        createParameter("skyline.popSize", "Bayesian Skyline population sizes", TIME_SCALE, 1.0, 0.0, Double.POSITIVE_INFINITY);
        createParameter("skyline.groupSize", "Bayesian Skyline group sizes");

        createParameter("skyride.popSize", "GMRF Bayesian skyride population sizes", TIME_SCALE, 1.0, Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY);
        createParameter("skyride.groupSize", "GMRF Bayesian skyride group sizes (for backward compatibility)");
        {
            final Parameter p = createParameter("skyride.precision", "GMRF Bayesian skyride precision", NONE, 1.0, 0.0, Double.POSITIVE_INFINITY);
            p.priorType = PriorType.GAMMA_PRIOR;
            p.gammaAlpha = 0.001;
            p.gammaBeta = 1000;
            p.priorFixed = true;
        }

        createParameter("demographic.popSize", "Extended Bayesian Skyline population sizes", TIME_SCALE, 1.0, 0.0, Double.POSITIVE_INFINITY);
        createParameter("demographic.indicators", "Extended Bayesian Skyline population switch");
        createScaleParameter("demographic.populationMean", "Extended Bayesian Skyline population prior mean", TIME_SCALE, 1, 0, Double.POSITIVE_INFINITY);
        {
            final Parameter p = createStatistic("demographic.populationSizeChanges", "Average number of population change points", true);
            p.priorType = PriorType.POISSON_PRIOR;
            p.poissonMean = Math.log(2);
        }
        createParameter("yule.birthRate", "Yule speciation process birth rate", BIRTH_RATE_SCALE, 1.0, 0.0, Double.POSITIVE_INFINITY);

        createParameter(BirthDeathModelParser.BIRTHDIFF_RATE_PARAM_NAME, "Birth-Death speciation process rate", BIRTH_RATE_SCALE, 1.0, 0.0, Double.POSITIVE_INFINITY);
        createParameter(BirthDeathModelParser.RELATIVE_DEATH_RATE_PARAM_NAME, "Death/Birth speciation process relative death rate", BIRTH_RATE_SCALE, 0.5, 0.0, 1.0);

        createScaleOperator("constant.popSize", demoWeights);
        createScaleOperator("exponential.popSize", demoWeights);
        createOperator("exponential.growthRate", OperatorType.RANDOM_WALK, 1.0, demoWeights);
        createScaleOperator("exponential.doublingTime", demoWeights);
        createScaleOperator("logistic.popSize", demoWeights);
        createScaleOperator("logistic.growthRate", demoWeights);
        createScaleOperator("logistic.doublingTime", demoWeights);
        createScaleOperator("logistic.t50", demoWeights);
        createScaleOperator("expansion.popSize", demoWeights);
        createScaleOperator("expansion.growthRate", demoWeights);
        createScaleOperator("expansion.doublingTime", demoWeights);
        createScaleOperator("expansion.ancestralProportion", demoWeights);
        createScaleOperator("skyline.popSize", demoWeights * 5);
        createOperator("skyline.groupSize", OperatorType.INTEGER_DELTA_EXCHANGE, 1.0, demoWeights * 2);

        createOperator("demographic.populationMean", OperatorType.SCALE, 0.9, demoWeights);
        createOperator("demographic.indicators", OperatorType.BITFLIP, 1, 2 * treeWeights);
        // hack pass distribution in name
        createOperator("demographic.popSize", "demographic.populationMeanDist", "", "demographic.popSize",
                "demographic.indicators", OperatorType.SAMPLE_NONACTIVE, 1, 5 * demoWeights);
        createOperator("demographic.scaleActive", "demographic.scaleActive", "", "demographic.popSize",
                "demographic.indicators", OperatorType.SCALE_WITH_INDICATORS, 0.5, 2 * demoWeights);

        createOperator("gmrfGibbsOperator", "gmrfGibbsOperator", "Gibbs sampler for GMRF", "skyride.popSize",
                "skyride.precision", OperatorType.GMRF_GIBBS_OPERATOR, 2, 2);

        createScaleOperator("yule.birthRate", demoWeights);

        createScaleOperator(BirthDeathModelParser.BIRTHDIFF_RATE_PARAM_NAME, demoWeights);
        createScaleOperator(BirthDeathModelParser.RELATIVE_DEATH_RATE_PARAM_NAME, demoWeights);

        // These are statistics which could have priors on...
        createStatistic("meanRate", "The mean rate of evolution over the whole tree", 0.0, Double.POSITIVE_INFINITY);
        createStatistic(RateStatistic.COEFFICIENT_OF_VARIATION, "The variation in rate of evolution over the whole tree",
                0.0, Double.POSITIVE_INFINITY);
        createStatistic("covariance",
                "The covariance in rates of evolution on each lineage with their ancestral lineages",
                Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY);

        // This only works if the partitions are of the same size...
//        createOperator("centeredMu", "Relative rates",
//                "Scales codon position rates relative to each other maintaining mean", "allMus",
//                OperatorType.CENTERED_SCALE, 0.75, rateWeights);
        createOperator("deltaMu", "Relative rates",
                "Changes partition relative rates relative to each other maintaining their mean", "allMus",
                OperatorType.DELTA_EXCHANGE, 0.75, rateWeights);

        createScaleOperator("clock.rate", rateWeights);
        createScaleOperator(ClockType.UCED_MEAN, rateWeights);
        createScaleOperator(ClockType.UCLD_MEAN, rateWeights);
        createScaleOperator(ClockType.UCLD_STDEV, rateWeights);

        createOperator("scaleRootRate", "treeModel.rootRate",
                "Scales root rate", "treeModel.rootRate",
                OperatorType.SCALE, 0.75, rateWeights);
        createOperator("scaleOneRate", "treeModel.nodeRates",
                "Scales one non-root rate", "treeModel.nodeRates",
                OperatorType.SCALE, 0.75, branchWeights);
        createOperator("scaleAllRates", "treeModel.allRates",
                "Scales all rates simultaneously", "treeModel.allRates",
                OperatorType.SCALE_ALL, 0.75, rateWeights);
        createOperator("scaleAllRatesIndependently", "treeModel.nodeRates",
                "Scales all non-root rates independently", "treeModel.nodeRates",
                OperatorType.SCALE_INDEPENDENTLY, 0.75, rateWeights);

        createOperator("upDownAllRatesHeights", "All rates and heights",
                "Scales all rates inversely to node heights of the tree", "treeModel.allRates",
                "treeModel.allInternalNodeHeights", OperatorType.UP_DOWN, 0.75, branchWeights);
        createScaleOperator("branchRates.var", rateWeights);

        createOperator("swapBranchRateCategories", "branchRates.categories",
                "Performs a swap of branch rate categories", "branchRates.categories",
                OperatorType.SWAP, 1, branchWeights / 3);
        createOperator("randomWalkBranchRateCategories", "branchRates.categories",
                "Performs an integer random walk of branch rate categories", "branchRates.categories",
                OperatorType.INTEGER_RANDOM_WALK, 1, branchWeights / 3);
        createOperator("unformBranchRateCategories", "branchRates.categories",
                "Performs an integer uniform draw of branch rate categories", "branchRates.categories",
                OperatorType.INTEGER_UNIFORM, 1, branchWeights / 3);

        createScaleOperator(ClockType.LOCAL_CLOCK + "." + "rates", treeWeights);
        createOperator(ClockType.LOCAL_CLOCK + "." + "changes", OperatorType.BITFLIP, 1, treeWeights);
        createOperator("treeBitMove", "Tree", "Swaps the rates and change locations of local clocks", "tree",
                OperatorType.TREE_BIT_MOVE, -1.0, treeWeights);

        createScaleOperator("treeModel.rootHeight", demoWeights);
        createOperator("uniformHeights", "Internal node heights", "Draws new internal node heights uniformally",
                "treeModel.internalNodeHeights", OperatorType.UNIFORM, -1, branchWeights);

        createOperator("upDownRateHeights", "Substitution rate and heights",
                "Scales substitution rates inversely to node heights of the tree", "clock.rate",
                "treeModel.allInternalNodeHeights", OperatorType.UP_DOWN, 0.75, rateWeights);
        createOperator("upDownUCEDMeanHeights", "UCED mean and heights",
                "Scales UCED mean inversely to node heights of the tree", ClockType.UCED_MEAN,
                "treeModel.allInternalNodeHeights", OperatorType.UP_DOWN, 0.75, rateWeights);
        createOperator("upDownUCLDMeanHeights", "UCLD mean and heights",
                "Scales UCLD mean inversely to node heights of the tree", ClockType.UCLD_MEAN,
                "treeModel.allInternalNodeHeights", OperatorType.UP_DOWN, 0.75, rateWeights);

        createOperator("subtreeSlide", "Tree", "Performs the subtree-slide rearrangement of the tree", "tree",
                OperatorType.SUBTREE_SLIDE, 1.0, treeWeights);
        createOperator("narrowExchange", "Tree", "Performs local rearrangements of the tree", "tree",
                OperatorType.NARROW_EXCHANGE, -1, treeWeights);
        createOperator("wideExchange", "Tree", "Performs global rearrangements of the tree", "tree",
                OperatorType.WIDE_EXCHANGE, -1, demoWeights);
        createOperator("wilsonBalding", "Tree", "Performs the Wilson-Balding rearrangement of the tree", "tree",
                OperatorType.WILSON_BALDING, -1, demoWeights);
    }
    
    // use override method getParameter(String name) in PartitionModel containing prefix	
    public void selectParameters(List<Parameter> params) {
        if (options.hasData()) {

            // if not fixed then do mutation rate move and up/down move
            boolean fixed = options.isFixedSubstitutionRate();
            Parameter rateParam;

            switch (options.clockType) {
                case STRICT_CLOCK:
                    rateParam = getParameter("clock.rate");
                    if (!fixed) params.add(rateParam);
                    break;

                case UNCORRELATED_EXPONENTIAL:
                    rateParam = getParameter(ClockType.UCED_MEAN);
                    if (!fixed) params.add(rateParam);
                    break;

                case UNCORRELATED_LOGNORMAL:
                    rateParam = getParameter(ClockType.UCLD_MEAN);
                    if (!fixed) params.add(rateParam);
                    params.add(getParameter(ClockType.UCLD_STDEV));
                    break;

                case AUTOCORRELATED_LOGNORMAL:
                    rateParam = getParameter("treeModel.rootRate");
                    if (!fixed) params.add(rateParam);
                    params.add(getParameter("branchRates.var"));
                    break;

                case RANDOM_LOCAL_CLOCK:
                    rateParam = getParameter("clock.rate");
                    if (!fixed) params.add(rateParam);
                    break;

                default:
                    throw new IllegalArgumentException("Unknown clock model");
            }

            /*if (clockType == ClockType.STRICT_CLOCK || clockType == ClockType.RANDOM_LOCAL_CLOCK) {
				rateParam = getParameter("clock.rate");
				if (!fixed) params.add(rateParam);
			} else {
				if (clockType == ClockType.UNCORRELATED_EXPONENTIAL) {
					rateParam = getParameter("uced.mean");
					if (!fixed) params.add(rateParam);
				} else if (clockType == ClockType.UNCORRELATED_LOGNORMAL) {
					rateParam = getParameter("ucld.mean");
					if (!fixed) params.add(rateParam);
					params.add(getParameter("ucld.stdev"));
				} else {
					throw new IllegalArgumentException("Unknown clock model");
				}
			}*/

            rateParam.isFixed = fixed;

        }

        if (options.nodeHeightPrior == TreePrior.CONSTANT) {
            params.add(getParameter("constant.popSize"));
        } else if (options.nodeHeightPrior == TreePrior.EXPONENTIAL) {
            params.add(getParameter("exponential.popSize"));
            if (options.parameterization == GROWTH_RATE) {
                params.add(getParameter("exponential.growthRate"));
            } else {
                params.add(getParameter("exponential.doublingTime"));
            }
        } else if (options.nodeHeightPrior == TreePrior.LOGISTIC) {
            params.add(getParameter("logistic.popSize"));
            if (options.parameterization == GROWTH_RATE) {
                params.add(getParameter("logistic.growthRate"));
            } else {
                params.add(getParameter("logistic.doublingTime"));
            }
            params.add(getParameter("logistic.t50"));
        } else if (options.nodeHeightPrior == TreePrior.EXPANSION) {
            params.add(getParameter("expansion.popSize"));
            if (options.parameterization == GROWTH_RATE) {
                params.add(getParameter("expansion.growthRate"));
            } else {
                params.add(getParameter("expansion.doublingTime"));
            }
            params.add(getParameter("expansion.ancestralProportion"));
        } else if (options.nodeHeightPrior == TreePrior.SKYLINE) {
            params.add(getParameter("skyline.popSize"));
        } else if (options.nodeHeightPrior == TreePrior.EXTENDED_SKYLINE) {
            params.add(getParameter("demographic.populationSizeChanges"));
            params.add(getParameter("demographic.populationMean"));
        } else if (options.nodeHeightPrior == TreePrior.GMRF_SKYRIDE) {
//            params.add(getParameter("skyride.popSize"));
            params.add(getParameter("skyride.precision"));
        } else if (options.nodeHeightPrior == TreePrior.YULE) {
            params.add(getParameter("yule.birthRate"));
        } else if (options.nodeHeightPrior == TreePrior.BIRTH_DEATH) {
            params.add(getParameter(BirthDeathModelParser.BIRTHDIFF_RATE_PARAM_NAME));
            params.add(getParameter(BirthDeathModelParser.RELATIVE_DEATH_RATE_PARAM_NAME));
        }

        params.add(getParameter("treeModel.rootHeight"));
    }
    
    // use override method getParameter(String name) and getOperator(String name) in PartitionModel containing prefix	
    public void selectOperators(List<Operator> ops) { 
        if (options.hasData()) {

            if (!options.isFixedSubstitutionRate()) {
                switch (options.clockType) {
                    case STRICT_CLOCK:
                        ops.add(getOperator("clock.rate"));
                        ops.add(getOperator("upDownRateHeights"));
                        break;

                    case UNCORRELATED_EXPONENTIAL:
                        ops.add(getOperator(ClockType.UCED_MEAN));
                        ops.add(getOperator("upDownUCEDMeanHeights"));
                        ops.add(getOperator("swapBranchRateCategories"));
                        ops.add(getOperator("randomWalkBranchRateCategories"));
                        ops.add(getOperator("unformBranchRateCategories"));
                        break;

                    case UNCORRELATED_LOGNORMAL:
                        ops.add(getOperator(ClockType.UCLD_MEAN));
                        ops.add(getOperator(ClockType.UCLD_STDEV));
                        ops.add(getOperator("upDownUCLDMeanHeights"));
                        ops.add(getOperator("swapBranchRateCategories"));
                        ops.add(getOperator("randomWalkBranchRateCategories"));
                        ops.add(getOperator("unformBranchRateCategories"));
                        break;

                    case AUTOCORRELATED_LOGNORMAL:
                        ops.add(getOperator("scaleRootRate"));
                        ops.add(getOperator("scaleOneRate"));
                        ops.add(getOperator("scaleAllRates"));
                        ops.add(getOperator("scaleAllRatesIndependently"));
                        ops.add(getOperator("upDownAllRatesHeights"));
                        ops.add(getOperator("branchRates.var"));
                        break;

                    case RANDOM_LOCAL_CLOCK:
                        ops.add(getOperator("clock.rate"));
                        ops.add(getOperator("upDownRateHeights"));
                        ops.add(getOperator(ClockType.LOCAL_CLOCK + "." + "rates"));
                        ops.add(getOperator(ClockType.LOCAL_CLOCK + "." + "changes"));
                        ops.add(getOperator("treeBitMove"));
                        break;

                    default:
                        throw new IllegalArgumentException("Unknown clock model");
                }
            } else {
                switch (options.clockType) {
                    case STRICT_CLOCK:
                        // no parameter to operator on
                        break;

                    case UNCORRELATED_EXPONENTIAL:
                        ops.add(getOperator("swapBranchRateCategories"));
                        ops.add(getOperator("randomWalkBranchRateCategories"));
                        ops.add(getOperator("unformBranchRateCategories"));
                        break;

                    case UNCORRELATED_LOGNORMAL:
                        ops.add(getOperator(ClockType.UCLD_STDEV));
                        ops.add(getOperator("swapBranchRateCategories"));
                        ops.add(getOperator("randomWalkBranchRateCategories"));
                        ops.add(getOperator("unformBranchRateCategories"));
                        break;

                    case AUTOCORRELATED_LOGNORMAL:
                        ops.add(getOperator("scaleOneRate"));
                        ops.add(getOperator("scaleAllRatesIndependently"));
                        ops.add(getOperator("branchRates.var"));
                        break;

                    case RANDOM_LOCAL_CLOCK:
                        ops.add(getOperator(ClockType.LOCAL_CLOCK + "." + "rates"));
                        ops.add(getOperator(ClockType.LOCAL_CLOCK + "." + "changes"));
                        ops.add(getOperator("treeBitMove"));
                        break;

                    default:
                        throw new IllegalArgumentException("Unknown clock model");
                }
            }
        }

        if (options.nodeHeightPrior == TreePrior.CONSTANT) {
            ops.add(getOperator("constant.popSize"));
        } else if (options.nodeHeightPrior == TreePrior.EXPONENTIAL) {
            ops.add(getOperator("exponential.popSize"));
            if (options.parameterization == GROWTH_RATE) {
                ops.add(getOperator("exponential.growthRate"));
            } else {
                ops.add(getOperator("exponential.doublingTime"));
            }
        } else if (options.nodeHeightPrior == TreePrior.LOGISTIC) {
            ops.add(getOperator("logistic.popSize"));
            if (options.parameterization == GROWTH_RATE) {
                ops.add(getOperator("logistic.growthRate"));
            } else {
                ops.add(getOperator("logistic.doublingTime"));
            }
            ops.add(getOperator("logistic.t50"));
        } else if (options.nodeHeightPrior == TreePrior.EXPANSION) {
            ops.add(getOperator("expansion.popSize"));
            if (options.parameterization == GROWTH_RATE) {
                ops.add(getOperator("expansion.growthRate"));
            } else {
                ops.add(getOperator("expansion.doublingTime"));
            }
            ops.add(getOperator("expansion.ancestralProportion"));
        } else if (options.nodeHeightPrior == TreePrior.SKYLINE) {
            ops.add(getOperator("skyline.popSize"));
            ops.add(getOperator("skyline.groupSize"));
        } else if (options.nodeHeightPrior == TreePrior.GMRF_SKYRIDE) {
            ops.add(getOperator("gmrfGibbsOperator"));
        } else if (options.nodeHeightPrior == TreePrior.EXTENDED_SKYLINE) {
            ops.add(getOperator("demographic.populationMean"));
            ops.add(getOperator("demographic.popSize"));
            ops.add(getOperator("demographic.indicators"));
            ops.add(getOperator("demographic.scaleActive"));
        } else if (options.nodeHeightPrior == TreePrior.YULE) {
            ops.add(getOperator("yule.birthRate"));
        } else if (options.nodeHeightPrior == TreePrior.BIRTH_DEATH) {
            ops.add(getOperator(BirthDeathModelParser.BIRTHDIFF_RATE_PARAM_NAME));
            ops.add(getOperator(BirthDeathModelParser.RELATIVE_DEATH_RATE_PARAM_NAME));
        }

        ops.add(getOperator("treeModel.rootHeight"));
        ops.add(getOperator("uniformHeights"));

        // if not a fixed tree then sample tree space
        if (!options.fixedTree) {
            ops.add(getOperator("subtreeSlide"));
            ops.add(getOperator("narrowExchange"));
            ops.add(getOperator("wideExchange"));
            ops.add(getOperator("wilsonBalding"));
        }

    }
    
    
    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }

    public List<Operator> getOperators() {
        List<Operator> operators = new ArrayList<Operator>();

        switch (dataType.getType()) {
            case DataType.NUCLEOTIDES:

                if (getCodonPartitionCount() > 1 && unlinkedSubstitutionModel) {
                    if (codonHeteroPattern.equals("123")) {
                        switch (nucSubstitutionModel) {
                            case HKY:
                                operators.add(getOperator("CP1.kappa"));
                                operators.add(getOperator("CP2.kappa"));
                                operators.add(getOperator("CP3.kappa"));
                                break;

                            case GTR:
                                for (int i = 1; i <= 3; i++) {
                                    for (String rateName : GTR_RATE_NAMES) {
                                        operators.add(getOperator("CP" + i + "." + rateName));
                                    }
                                }
                                break;

                            default:
                                throw new IllegalArgumentException("Unknown nucleotides substitution model");
                        }

                        if (frequencyPolicy == FrequencyPolicy.ESTIMATED) {
                            if (getCodonPartitionCount() > 1 && unlinkedSubstitutionModel) {
                                operators.add(getOperator("CP1.frequencies"));
                                operators.add(getOperator("CP2.frequencies"));
                                operators.add(getOperator("CP3.frequencies"));
                            } else {
                                operators.add(getOperator("frequencies"));
                            }
                        }
                    } else if (codonHeteroPattern.equals("112")) {
                        switch (nucSubstitutionModel) {
                            case HKY:
                                operators.add(getOperator("CP1+2.kappa"));
                                operators.add(getOperator("CP3.kappa"));
                                break;

                            case GTR:
                                for (String rateName : GTR_RATE_NAMES) {
                                    operators.add(getOperator("CP1+2." + rateName));
                                }
                                for (String rateName : GTR_RATE_NAMES) {
                                    operators.add(getOperator("CP3." + rateName));
                                }
                                break;

                            default:
                                throw new IllegalArgumentException("Unknown nucleotides substitution model");
                        }
                        if (frequencyPolicy == FrequencyPolicy.ESTIMATED) {
                            if (getCodonPartitionCount() > 1 && unlinkedSubstitutionModel) {
                                operators.add(getOperator("CP1+2.frequencies"));
                                operators.add(getOperator("CP3.frequencies"));
                            } else {
                                operators.add(getOperator("frequencies"));
                            }
                        }

                    } else {
                        throw new IllegalArgumentException("codonHeteroPattern must be one of '111', '112' or '123'");
                    }
                } else { // no codon partitioning
                    switch (nucSubstitutionModel) {
                        case HKY:
                            operators.add(getOperator("kappa"));
                            break;

                        case GTR:
                            for (String rateName : GTR_RATE_NAMES) {
                                operators.add(getOperator(rateName));
                            }
                            break;

                        default:
                            throw new IllegalArgumentException("Unknown nucleotides substitution model");
                    }
                    if (frequencyPolicy == FrequencyPolicy.ESTIMATED) {
                        operators.add(getOperator("frequencies"));
                    }
                }
                break;

            case DataType.AMINO_ACIDS:
                break;

            case DataType.TWO_STATES:
            case DataType.COVARION:
                switch (binarySubstitutionModel) {
                    case BIN_SIMPLE:
                        break;

                    case BIN_COVARION:
                        operators.add(getOperator("bcov.alpha"));
                        operators.add(getOperator("bcov.s"));
                        operators.add(getOperator("bcov.frequencies"));
                        operators.add(getOperator("bcov.hfrequencies"));
                        break;

                    default:
                        throw new IllegalArgumentException("Unknown binary substitution model");
                }
                break;

            default:
                throw new IllegalArgumentException("Unknown data type");
        }

        // if gamma do shape move
        if (gammaHetero) {
            if (getCodonPartitionCount() > 1 && unlinkedHeterogeneityModel) {
                if (codonHeteroPattern.equals("123")) {
                    operators.add(getOperator("CP1.alpha"));
                    operators.add(getOperator("CP2.alpha"));
                    operators.add(getOperator("CP3.alpha"));
                } else if (codonHeteroPattern.equals("112")) {
                    operators.add(getOperator("CP1+2.alpha"));
                    operators.add(getOperator("CP3.alpha"));
                } else {
                    throw new IllegalArgumentException("codonHeteroPattern must be one of '111', '112' or '123'");
                }
            } else {
                operators.add(getOperator("alpha"));
            }
        }
        // if pinv do pinv move
        if (invarHetero) {
            if (getCodonPartitionCount() > 1 && unlinkedHeterogeneityModel) {
                if (codonHeteroPattern.equals("123")) {
                    operators.add(getOperator("CP1.pInv"));
                    operators.add(getOperator("CP2.alpha"));
                    operators.add(getOperator("CP3.pInv"));
                } else if (codonHeteroPattern.equals("112")) {
                    operators.add(getOperator("CP1+2.pInv"));
                    operators.add(getOperator("CP3.pInv"));
                } else {
                    throw new IllegalArgumentException("codonHeteroPattern must be one of '111', '112' or '123'");
                }
            } else {
                operators.add(getOperator("pInv"));
            }
        }

        return operators;
    }

    /**
     * @param includeRelativeRates true if relative rate parameters should be added
     * @return a list of parameters that are required
     */
    List<Parameter> getParameters(boolean includeRelativeRates) {

        List<Parameter> params = new ArrayList<Parameter>();

        switch (dataType.getType()) {
            case DataType.NUCLEOTIDES:
                if (getCodonPartitionCount() > 1 && unlinkedSubstitutionModel) {
                    if (codonHeteroPattern.equals("123")) {
                        switch (nucSubstitutionModel) {
                            case HKY:
                                params.add(getParameter("CP1.kappa"));
                                params.add(getParameter("CP2.kappa"));
                                params.add(getParameter("CP3.kappa"));
                                break;
                            case GTR:
                                for (int i = 1; i <= getCodonPartitionCount(); i++) {
                                    for (String rateName : GTR_RATE_NAMES) {
                                        params.add(getParameter("CP" + i + "." + rateName));
                                    }
                                }
                                break;

                            default:
                                throw new IllegalArgumentException("Unknown nucleotides substitution model");
                        }
                        params.add(getParameter("CP1.mu"));
                        params.add(getParameter("CP2.mu"));
                        params.add(getParameter("CP3.mu"));
                    } else if (codonHeteroPattern.equals("112")) {
                        switch (nucSubstitutionModel) {
                            case HKY:
                                params.add(getParameter("CP1+2.kappa"));
                                params.add(getParameter("CP3.kappa"));
                                break;
                            case GTR:
                                for (String rateName : GTR_RATE_NAMES) {
                                    params.add(getParameter("CP1+2." + rateName));
                                }
                                for (String rateName : GTR_RATE_NAMES) {
                                    params.add(getParameter("CP3." + rateName));
                                }
                                break;

                            default:
                                throw new IllegalArgumentException("Unknown nucleotides substitution model");
                        }
                        params.add(getParameter("CP1+2.mu"));
                        params.add(getParameter("CP3.mu"));
                    } else {
                        throw new IllegalArgumentException("codonHeteroPattern must be one of '111', '112' or '123'");
                    }
                } else { // no codon partitioning
                    switch (nucSubstitutionModel) {
                        case HKY:
                            params.add(getParameter("kappa"));
                            break;
                        case GTR:
                            for (String rateName : GTR_RATE_NAMES) {
                                params.add(getParameter(rateName));
                            }
                            break;

                        default:
                            throw new IllegalArgumentException("Unknown nucleotides substitution model");
                    }
                    if (includeRelativeRates) {
                        params.add(getParameter("mu"));
                    }
                }
                break;

            case DataType.AMINO_ACIDS:
                if (includeRelativeRates) {
                    params.add(getParameter("mu"));
                }
                break;

            case DataType.TWO_STATES:
            case DataType.COVARION:
                switch (binarySubstitutionModel) {
                    case BIN_SIMPLE:
                        break;

                    case BIN_COVARION:
                        params.add(getParameter("bcov.alpha"));
                        params.add(getParameter("bcov.s"));
                        break;

                    default:
                        throw new IllegalArgumentException("Unknown binary substitution model");
                }
                if (includeRelativeRates) {
                    params.add(getParameter("mu"));
                }
                break;

            default:
                throw new IllegalArgumentException("Unknown data type");
        }

        // if gamma do shape move
        if (gammaHetero) {
            if (getCodonPartitionCount() > 1 && unlinkedHeterogeneityModel) {
                if (codonHeteroPattern.equals("123")) {
                    params.add(getParameter("CP1.alpha"));
                    params.add(getParameter("CP2.alpha"));
                    params.add(getParameter("CP3.alpha"));
                } else if (codonHeteroPattern.equals("112")) {
                    params.add(getParameter("CP1+2.alpha"));
                    params.add(getParameter("CP3.alpha"));
                } else {
                    throw new IllegalArgumentException("codonHeteroPattern must be one of '111', '112' or '123'");
                }
            } else {
                params.add(getParameter("alpha"));
            }
        }
        // if pinv do pinv move
        if (invarHetero) {
            if (getCodonPartitionCount() > 1 && unlinkedHeterogeneityModel) {
                if (codonHeteroPattern.equals("123")) {
                    params.add(getParameter("CP1.pInv"));
                    params.add(getParameter("CP2.pInv"));
                    params.add(getParameter("CP3.pInv"));
                } else if (codonHeteroPattern.equals("112")) {
                    params.add(getParameter("CP1+2.pInv"));
                    params.add(getParameter("CP3.pInv"));
                } else {
                    throw new IllegalArgumentException("codonHeteroPattern must be one of '111', '112' or '123'");
                }
            } else {
                params.add(getParameter("pInv"));
            }
        }
        if (frequencyPolicy == FrequencyPolicy.ESTIMATED) {
            if (getCodonPartitionCount() > 1 && unlinkedHeterogeneityModel) {
                if (codonHeteroPattern.equals("123")) {
                    params.add(getParameter("CP1.frequencies"));
                    params.add(getParameter("CP2.frequencies"));
                    params.add(getParameter("CP3.frequencies"));
                } else if (codonHeteroPattern.equals("112")) {
                    params.add(getParameter("CP1+2.frequencies"));
                    params.add(getParameter("CP3.frequencies"));
                } else {
                    throw new IllegalArgumentException("codonHeteroPattern must be one of '111', '112' or '123'");
                }
            } else {
                params.add(getParameter("frequencies"));
            }

        }

        return params;
    }

    public Parameter getParameter(String name) {

        if (name.startsWith(getName())) {
            name = name.substring(getName().length() + 1);
        }
        Parameter parameter = parameters.get(name);

        if (parameter == null) {
            throw new IllegalArgumentException("Parameter with name, " + name + ", is unknown");
        }

        parameter.setPrefix(getPrefix());

        return parameter;
    }

    public Operator getOperator(String name) {

        Operator operator = operators.get(name);

        if (operator == null) throw new IllegalArgumentException("Operator with name, " + name + ", is unknown");

        operator.setPrefix(getName());

        return operator;
    }

    public int getCodonPartitionCount() {
        if (codonHeteroPattern == null || codonHeteroPattern.equals("111")) {
            return 1;
        }
        if (codonHeteroPattern.equals("123")) {
            return 3;
        }
        if (codonHeteroPattern.equals("112")) {
            return 2;
        }
        throw new IllegalArgumentException("codonHeteroPattern must be one of '111', '112' or '123'");
    }

    public void addWeightsForPartition(DataPartition partition, int[] weights, int offset) {
        int n = partition.getSiteCount();
        int codonCount = n / 3;
        int remainder = n % 3;
        if (codonHeteroPattern == null || codonHeteroPattern.equals("111")) {
            weights[offset] += n;
            return;
        }
        if (codonHeteroPattern.equals("123")) {
            weights[offset] += codonCount + (remainder > 0 ? 1 : 0);
            weights[offset + 1] += codonCount + (remainder > 1 ? 1 : 0);
            weights[offset + 2] += codonCount;
            return;
        }
        if (codonHeteroPattern.equals("112")) {
            weights[offset] += codonCount * 2 + remainder; // positions 1 + 2
            weights[offset + 1] += codonCount; // position 3
            return;
        }
        throw new IllegalArgumentException("codonHeteroPattern must be one of '111', '112' or '123'");
    }

    public String toString() {
        return getName();
    }

    public NucModelType getNucSubstitutionModel() {
        return nucSubstitutionModel;
    }

    public void setNucSubstitutionModel(NucModelType nucSubstitutionModel) {
        this.nucSubstitutionModel = nucSubstitutionModel;
    }

    public AminoAcidModelType getAaSubstitutionModel() {
        return aaSubstitutionModel;
    }

    public void setAaSubstitutionModel(AminoAcidModelType aaSubstitutionModel) {
        this.aaSubstitutionModel = aaSubstitutionModel;
    }

    public int getBinarySubstitutionModel() {
        return binarySubstitutionModel;
    }

    public void setBinarySubstitutionModel(int binarySubstitutionModel) {
        this.binarySubstitutionModel = binarySubstitutionModel;
    }

    public FrequencyPolicy getFrequencyPolicy() {
        return frequencyPolicy;
    }

    public void setFrequencyPolicy(FrequencyPolicy frequencyPolicy) {
        this.frequencyPolicy = frequencyPolicy;
    }

    public boolean isGammaHetero() {
        return gammaHetero;
    }

    public void setGammaHetero(boolean gammaHetero) {
        this.gammaHetero = gammaHetero;
    }

    public int getGammaCategories() {
        return gammaCategories;
    }

    public void setGammaCategories(int gammaCategories) {
        this.gammaCategories = gammaCategories;
    }

    public boolean isInvarHetero() {
        return invarHetero;
    }

    public void setInvarHetero(boolean invarHetero) {
        this.invarHetero = invarHetero;
    }

    public String getCodonHeteroPattern() {
        return codonHeteroPattern;
    }

    public void setCodonHeteroPattern(String codonHeteroPattern) {
        this.codonHeteroPattern = codonHeteroPattern;
    }

    /**
     * @return true if the rate matrix parameters are unlinked across codon positions
     */
    public boolean isUnlinkedSubstitutionModel() {
        return unlinkedSubstitutionModel;
    }

    public void setUnlinkedSubstitutionModel(boolean unlinkedSubstitutionModel) {
        this.unlinkedSubstitutionModel = unlinkedSubstitutionModel;
    }

    public boolean isUnlinkedHeterogeneityModel() {
        return unlinkedHeterogeneityModel;
    }

    public void setUnlinkedHeterogeneityModel(boolean unlinkedHeterogeneityModel) {
        this.unlinkedHeterogeneityModel = unlinkedHeterogeneityModel;
    }

    public boolean isUnlinkedFrequencyModel() {
        return unlinkedFrequencyModel;
    }

    public void setUnlinkedFrequencyModel(boolean unlinkedFrequencyModel) {
        this.unlinkedFrequencyModel = unlinkedFrequencyModel;
    }

    public DataType getDataType() {
        return dataType;
    }

    public void setDataType(DataType dataType) {
        this.dataType = dataType;
    }

    public boolean isDolloModel() {
        return dolloModel;
    }

    public void setDolloModel(boolean dolloModel) {
        this.dolloModel = dolloModel;
    }

    public String getPrefix() {
        String prefix = "";
        if (options.getActivePartitionModels().size() > 1) {
            // There is more than one active partition model
            prefix += getName() + ".";
        }
        return prefix;
    }

    public String getPrefix(int codonPartitionNumber) {
        String prefix = "";
        if (options.getActivePartitionModels().size() > 1) {
            // There is more than one active partition model
            prefix += getName() + ".";
        }
        if (getCodonPartitionCount() > 1 && codonPartitionNumber > 0) {
            if (getCodonHeteroPattern().equals("123")) {
                prefix += "CP" + codonPartitionNumber + ".";
            } else if (getCodonHeteroPattern().equals("112")) {
                if (codonPartitionNumber == 1) {
                    prefix += "CP1+2.";
                } else {
                    prefix += "CP3.";
                }
            } else {
                throw new IllegalArgumentException("unsupported codon hetero pattern");
            }

        }
        return prefix;
    }

    // Instance variables

    private final BeautiOptions options;

    private NucModelType nucSubstitutionModel = NucModelType.HKY;
    private AminoAcidModelType aaSubstitutionModel = AminoAcidModelType.BLOSUM_62;
    private int binarySubstitutionModel = BeautiOptions.BIN_SIMPLE;

    private FrequencyPolicy frequencyPolicy = FrequencyPolicy.ESTIMATED;
    private boolean gammaHetero = false;
    private int gammaCategories = 4;
    private boolean invarHetero = false;
    private String codonHeteroPattern = null;
    private boolean unlinkedSubstitutionModel = false;
    private boolean unlinkedHeterogeneityModel = false;

    private boolean unlinkedFrequencyModel = false;

    private boolean dolloModel = false;

    public DataType dataType;
    public String name;

}