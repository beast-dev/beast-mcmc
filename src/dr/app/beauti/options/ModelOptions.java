package dr.app.beauti.options;

import dr.app.beauti.priorsPanel.PriorType;
import dr.evolution.util.TaxonList;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 */
public abstract class ModelOptions {

    HashMap<String, Parameter> parameters = new HashMap<String, Parameter>();
    HashMap<TaxonList, Parameter> statistics = new HashMap<TaxonList, Parameter>();
    HashMap<String, Operator> operators = new HashMap<String, Operator>();

    public static final String version = "1.5";
    public static final int YEARS = 0;
    public static final int MONTHS = 1;
    public static final int DAYS = 2;
    public static final int FORWARDS = 0;
    public static final int BACKWARDS = 1;
    public static final int NONE = -1;

    public static final int BIN_SIMPLE = 0;
    public static final int BIN_COVARION = 1;

    public static final int GROWTH_RATE = 0;
    public static final int DOUBLING_TIME = 1;
    public static final int CONSTANT_SKYLINE = 0;
    public static final int LINEAR_SKYLINE = 1;

    public static final int SKYRIDE_UNIFORM_SMOOTHING = 0;
    public static final int SKYRIDE_TIME_AWARE_SMOOTHING = 1;

    public static final int TIME_SCALE = 0;
    public static final int GROWTH_RATE_SCALE = 1;
    public static final int BIRTH_RATE_SCALE = 2;
    public static final int SUBSTITUTION_RATE_SCALE = 3;
    public static final int LOG_STDEV_SCALE = 4;
    public static final int SUBSTITUTION_PARAMETER_SCALE = 5;
    public static final int T50_SCALE = 6;
    public static final int UNITY_SCALE = 7;
    public static final int ROOT_RATE_SCALE = 8;
    public static final int LOG_VAR_SCALE = 9;

    public static final String[] GTR_RATE_NAMES = {"ac", "ag", "at", "cg", "gt"};
    protected static final String[] GTR_TRANSITIONS = {"A-C", "A-G", "A-T", "C-G", "G-T"};

    public static final double demoTuning = 0.75;
    public static final double demoWeights = 3.0;


    /**
     * Initialise parameters and operators, the prefix will be added by parameter.setPrefix(getPrefix()) or operator.setPrefix(getName()).
     */
    public void initAllParametersAndOperators() {

        initGlobalParaAndOpers();

//        initClockModelParaAndOpers();

//        initTreeModelParaAndOpers();        
//        initTreePriorParaAndOpers();        
    }

    public void initGlobalParaAndOpers() {
        double rateWeights = 3.0;

        // A vector of relative rates across all partitions...
        createParameter("allMus", "All the relative rates");

        // This only works if the partitions are of the same size...
//      createOperator("centeredMu", "Relative rates",
//              "Scales codon position rates relative to each other maintaining mean", "allMus",
//              OperatorType.CENTERED_SCALE, 0.75, rateWeights);
        createOperator("deltaMu", "Relative rates",
                "Changes partition relative rates relative to each other maintaining their mean", "allMus",
                OperatorType.DELTA_EXCHANGE, 0.75, rateWeights);
    }

//    public void initClockModelParaAndOpers() {
//    	double rateWeights = 3.0; 
//    	
//    	createParameter("clock.rate", "substitution rate", SUBSTITUTION_RATE_SCALE, 1.0, 0.0, Double.POSITIVE_INFINITY);
//        createParameter(ClockType.UCED_MEAN, "uncorrelated exponential relaxed clock mean", SUBSTITUTION_RATE_SCALE, 1.0, 0.0, Double.POSITIVE_INFINITY);
//        createParameter(ClockType.UCLD_MEAN, "uncorrelated lognormal relaxed clock mean", SUBSTITUTION_RATE_SCALE, 1.0, 0.0, Double.POSITIVE_INFINITY);
//        createParameter(ClockType.UCLD_STDEV, "uncorrelated lognormal relaxed clock stdev", LOG_STDEV_SCALE, 0.1, 0.0, Double.POSITIVE_INFINITY);
//
//        {
//            final Parameter p = createParameter("branchRates.var", "autocorrelated lognormal relaxed clock rate variance ", LOG_VAR_SCALE, 0.1, 0.0, Double.POSITIVE_INFINITY);
//            p.priorType = PriorType.GAMMA_PRIOR;
//            p.gammaAlpha = 1;
//            p.gammaBeta = 0.0001;
//        }
//
//        createScaleOperator("clock.rate", demoTuning, rateWeights);
//        createScaleOperator(ClockType.UCED_MEAN, demoTuning, rateWeights);
//        createScaleOperator(ClockType.UCLD_MEAN, demoTuning, rateWeights);
//        createScaleOperator(ClockType.UCLD_STDEV, demoTuning, rateWeights);
//        
//        createScaleOperator("branchRates.var", demoTuning, rateWeights);
//
//        createOperator("upDownRateHeights", "Substitution rate and heights",
//                "Scales substitution rates inversely to node heights of the tree", "clock.rate",
//                "treeModel.allInternalNodeHeights", OperatorType.UP_DOWN, 0.75, rateWeights);
//        createOperator("upDownUCEDMeanHeights", "UCED mean and heights",
//                "Scales UCED mean inversely to node heights of the tree", ClockType.UCED_MEAN,
//                "treeModel.allInternalNodeHeights", OperatorType.UP_DOWN, 0.75, rateWeights);
//        createOperator("upDownUCLDMeanHeights", "UCLD mean and heights",
//                "Scales UCLD mean inversely to node heights of the tree", ClockType.UCLD_MEAN,
//                "treeModel.allInternalNodeHeights", OperatorType.UP_DOWN, 0.75, rateWeights);
//    }

//    public void initTreeModelParaAndOpers() {
//        double branchWeights = 30.0;
//        double treeWeights = 15.0;
//        double rateWeights = 3.0;
//        
//        createParameter("tree", "The tree");
//        createParameter("treeModel.internalNodeHeights", "internal node heights of the tree (except the root)");
//        createParameter("treeModel.allInternalNodeHeights", "internal node heights of the tree");
//        createParameter("treeModel.rootHeight", "root height of the tree", true, 1.0, 0.0, Double.POSITIVE_INFINITY);
//
//        createParameter("branchRates.categories", "relaxed clock branch rate categories");
//        createParameter(ClockType.LOCAL_CLOCK + "." + "rates", "random local clock rates", SUBSTITUTION_RATE_SCALE, 1.0, 0.0, Double.POSITIVE_INFINITY);
//        createParameter(ClockType.LOCAL_CLOCK + "." + "changes", "random local clock rate change indicator");
//
//        {
//            final Parameter p = createParameter("treeModel.rootRate", "autocorrelated lognormal relaxed clock root rate", ROOT_RATE_SCALE, 1.0, 0.0, Double.POSITIVE_INFINITY);
//            p.priorType = PriorType.GAMMA_PRIOR;
//            p.gammaAlpha = 1;
//            p.gammaBeta = 0.0001;
//        }
//        createParameter("treeModel.nodeRates", "autocorrelated lognormal relaxed clock non-root rates", SUBSTITUTION_RATE_SCALE, 1.0, 0.0, Double.POSITIVE_INFINITY);
//        createParameter("treeModel.allRates", "autocorrelated lognormal relaxed clock all rates", SUBSTITUTION_RATE_SCALE, 1.0, 0.0, Double.POSITIVE_INFINITY);
//
//        // These are statistics which could have priors on...
//        createStatistic("meanRate", "The mean rate of evolution over the whole tree", 0.0, Double.POSITIVE_INFINITY);
//        createStatistic(RateStatistic.COEFFICIENT_OF_VARIATION, "The variation in rate of evolution over the whole tree",
//                0.0, Double.POSITIVE_INFINITY);
//        createStatistic("covariance",
//                "The covariance in rates of evolution on each lineage with their ancestral lineages",
//                Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY);
//
//        createOperator("scaleRootRate", "treeModel.rootRate",
//                "Scales root rate", "treeModel.rootRate",
//                OperatorType.SCALE, 0.75, rateWeights);
//        createOperator("scaleOneRate", "treeModel.nodeRates",
//                "Scales one non-root rate", "treeModel.nodeRates",
//                OperatorType.SCALE, 0.75, branchWeights);
//        createOperator("scaleAllRates", "treeModel.allRates",
//                "Scales all rates simultaneously", "treeModel.allRates",
//                OperatorType.SCALE_ALL, 0.75, rateWeights);
//        createOperator("scaleAllRatesIndependently", "treeModel.nodeRates",
//                "Scales all non-root rates independently", "treeModel.nodeRates",
//                OperatorType.SCALE_INDEPENDENTLY, 0.75, rateWeights);
//
//        createOperator("upDownAllRatesHeights", "All rates and heights",
//                "Scales all rates inversely to node heights of the tree", "treeModel.allRates",
//                "treeModel.allInternalNodeHeights", OperatorType.UP_DOWN, 0.75, branchWeights);
//        createScaleOperator("branchRates.var", demoTuning, rateWeights);
//
//        createOperator("swapBranchRateCategories", "branchRates.categories",
//                "Performs a swap of branch rate categories", "branchRates.categories",
//                OperatorType.SWAP, 1, branchWeights / 3);
//        createOperator("randomWalkBranchRateCategories", "branchRates.categories",
//                "Performs an integer random walk of branch rate categories", "branchRates.categories",
//                OperatorType.INTEGER_RANDOM_WALK, 1, branchWeights / 3);
//        createOperator("unformBranchRateCategories", "branchRates.categories",
//                "Performs an integer uniform draw of branch rate categories", "branchRates.categories",
//                OperatorType.INTEGER_UNIFORM, 1, branchWeights / 3);
//        
//        createScaleOperator(ClockType.LOCAL_CLOCK + "." + "rates", demoTuning, treeWeights);
//        createOperator(ClockType.LOCAL_CLOCK + "." + "changes", OperatorType.BITFLIP, 1, treeWeights);
//        createOperator("treeBitMove", "Tree", "Swaps the rates and change locations of local clocks", "tree",
//                OperatorType.TREE_BIT_MOVE, -1.0, treeWeights);
//
//        createScaleOperator("treeModel.rootHeight", demoTuning, demoWeights);
//        createOperator("uniformHeights", "Internal node heights", "Draws new internal node heights uniformally",
//                "treeModel.internalNodeHeights", OperatorType.UNIFORM, -1, branchWeights);
//
//        createOperator("upDownRateHeights", "Substitution rate and heights",
//                "Scales substitution rates inversely to node heights of the tree", "clock.rate",
//                "treeModel.allInternalNodeHeights", OperatorType.UP_DOWN, 0.75, rateWeights);
//        createOperator("upDownUCEDMeanHeights", "UCED mean and heights",
//                "Scales UCED mean inversely to node heights of the tree", ClockType.UCED_MEAN,
//                "treeModel.allInternalNodeHeights", OperatorType.UP_DOWN, 0.75, rateWeights);
//        createOperator("upDownUCLDMeanHeights", "UCLD mean and heights",
//                "Scales UCLD mean inversely to node heights of the tree", ClockType.UCLD_MEAN,
//                "treeModel.allInternalNodeHeights", OperatorType.UP_DOWN, 0.75, rateWeights);
//
//        createOperator("subtreeSlide", "Tree", "Performs the subtree-slide rearrangement of the tree", "tree",
//                OperatorType.SUBTREE_SLIDE, 1.0, treeWeights);
//        createOperator("narrowExchange", "Tree", "Performs local rearrangements of the tree", "tree",
//                OperatorType.NARROW_EXCHANGE, -1, treeWeights);
//        createOperator("wideExchange", "Tree", "Performs global rearrangements of the tree", "tree",
//                OperatorType.WIDE_EXCHANGE, -1, demoWeights);
//        createOperator("wilsonBalding", "Tree", "Performs the Wilson-Balding rearrangement of the tree", "tree",
//                OperatorType.WILSON_BALDING, -1, demoWeights);
//    }

//    public void initTreePriorParaAndOpers() {
//    	double treeWeights = 15.0;
//    	    	
//        createScaleParameter("constant.popSize", "coalescent population size parameter", TIME_SCALE, 1.0, 0.0, Double.POSITIVE_INFINITY);
//
//        createScaleParameter("exponential.popSize", "coalescent population size parameter", TIME_SCALE, 1.0, 0.0, Double.POSITIVE_INFINITY);
//        createParameter("exponential.growthRate", "coalescent growth rate parameter", GROWTH_RATE_SCALE, 0.0, Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY);
//        createParameter("exponential.doublingTime", "coalescent doubling time parameter", TIME_SCALE, 0.5, 0.0, Double.POSITIVE_INFINITY);
//        createScaleParameter("logistic.popSize", "coalescent population size parameter", TIME_SCALE, 1.0, 0.0, Double.POSITIVE_INFINITY);
//        createParameter("logistic.growthRate", "coalescent logistic growth rate parameter", GROWTH_RATE_SCALE, 0.001, 0.0, Double.POSITIVE_INFINITY);
//        createParameter("logistic.doublingTime", "coalescent doubling time parameter", TIME_SCALE, 0.5, 0.0, Double.POSITIVE_INFINITY);
//        createParameter("logistic.t50", "logistic shape parameter", T50_SCALE, 0.1, 0.0, Double.POSITIVE_INFINITY);
//        createScaleParameter("expansion.popSize", "coalescent population size parameter", TIME_SCALE, 1.0, 0.0, Double.POSITIVE_INFINITY);
//        createParameter("expansion.growthRate", "coalescent logistic growth rate parameter", GROWTH_RATE_SCALE, 0.001, 0.0, Double.POSITIVE_INFINITY);
//        createParameter("expansion.doublingTime", "coalescent doubling time parameter", TIME_SCALE, 0.5, 0.0, Double.POSITIVE_INFINITY);
//        createParameter("expansion.ancestralProportion", "ancestral population proportion", NONE, 0.1, 0.0, 1.0);
//        createParameter("skyline.popSize", "Bayesian Skyline population sizes", TIME_SCALE, 1.0, 0.0, Double.POSITIVE_INFINITY);
//        createParameter("skyline.groupSize", "Bayesian Skyline group sizes");
//
//        createParameter("skyride.popSize", "GMRF Bayesian skyride population sizes", TIME_SCALE, 1.0, Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY);
//        createParameter("skyride.groupSize", "GMRF Bayesian skyride group sizes (for backward compatibility)");
//        {
//            final Parameter p = createParameter("skyride.precision", "GMRF Bayesian skyride precision", NONE, 1.0, 0.0, Double.POSITIVE_INFINITY);
//            p.priorType = PriorType.GAMMA_PRIOR;
//            p.gammaAlpha = 0.001;
//            p.gammaBeta = 1000;
//            p.priorFixed = true;
//        }
//
//        createParameter("demographic.popSize", "Extended Bayesian Skyline population sizes", TIME_SCALE, 1.0, 0.0, Double.POSITIVE_INFINITY);
//        createParameter("demographic.indicators", "Extended Bayesian Skyline population switch");
//        createScaleParameter("demographic.populationMean", "Extended Bayesian Skyline population prior mean", TIME_SCALE, 1, 0, Double.POSITIVE_INFINITY);
//        {
//            final Parameter p = createStatistic("demographic.populationSizeChanges", "Average number of population change points", true);
//            p.priorType = PriorType.POISSON_PRIOR;
//            p.poissonMean = Math.log(2);
//        }
//        createParameter("yule.birthRate", "Yule speciation process birth rate", BIRTH_RATE_SCALE, 1.0, 0.0, Double.POSITIVE_INFINITY);
//
//        createParameter(BirthDeathModelParser.BIRTHDIFF_RATE_PARAM_NAME, "Birth-Death speciation process rate", BIRTH_RATE_SCALE, 1.0, 0.0, Double.POSITIVE_INFINITY);
//        createParameter(BirthDeathModelParser.RELATIVE_DEATH_RATE_PARAM_NAME, "Death/Birth speciation process relative death rate", BIRTH_RATE_SCALE, 0.5, 0.0, 1.0);
//
//        createScaleOperator("constant.popSize", demoTuning, demoWeights);
//        createScaleOperator("exponential.popSize", demoTuning, demoWeights);
//        createOperator("exponential.growthRate", OperatorType.RANDOM_WALK, 1.0, demoWeights);
//        createScaleOperator("exponential.doublingTime", demoTuning, demoWeights);
//        createScaleOperator("logistic.popSize", demoTuning, demoWeights);
//        createScaleOperator("logistic.growthRate", demoTuning, demoWeights);
//        createScaleOperator("logistic.doublingTime", demoTuning, demoWeights);
//        createScaleOperator("logistic.t50", demoTuning, demoWeights);
//        createScaleOperator("expansion.popSize", demoTuning, demoWeights);
//        createScaleOperator("expansion.growthRate", demoTuning, demoWeights);
//        createScaleOperator("expansion.doublingTime", demoTuning, demoWeights);
//        createScaleOperator("expansion.ancestralProportion", demoTuning, demoWeights);
//        createScaleOperator("skyline.popSize", demoTuning, demoWeights * 5);
//        createOperator("skyline.groupSize", OperatorType.INTEGER_DELTA_EXCHANGE, 1.0, demoWeights * 2);
//
//        createOperator("demographic.populationMean", OperatorType.SCALE, 0.9, demoWeights);
//        createOperator("demographic.indicators", OperatorType.BITFLIP, 1, 2 * treeWeights);
//        
//        // hack pass distribution in name
//        createOperator("demographic.popSize", "demographic.populationMeanDist", "", "demographic.popSize",
//                "demographic.indicators", OperatorType.SAMPLE_NONACTIVE, 1, 5 * demoWeights);
//        createOperator("demographic.scaleActive", "demographic.scaleActive", "", "demographic.popSize",
//                "demographic.indicators", OperatorType.SCALE_WITH_INDICATORS, 0.5, 2 * demoWeights);
//
//        createOperator("gmrfGibbsOperator", "gmrfGibbsOperator", "Gibbs sampler for GMRF", "skyride.popSize",
//                "skyride.precision", OperatorType.GMRF_GIBBS_OPERATOR, 2, 2);
//
//        createScaleOperator("yule.birthRate", demoTuning, demoWeights);
//
//        createScaleOperator(BirthDeathModelParser.BIRTHDIFF_RATE_PARAM_NAME, demoTuning, demoWeights);
//        createScaleOperator(BirthDeathModelParser.RELATIVE_DEATH_RATE_PARAM_NAME, demoTuning, demoWeights);
//
//    }

    public void createOperator(String parameterName, OperatorType type, double tuning, double weight) {
        Parameter parameter = getParameter(parameterName);
        operators.put(parameterName, new Operator(parameterName, "", parameter, type, tuning, weight));
    }

    public void createOperator(String key, String name, String description, String parameterName, OperatorType type, double tuning, double weight) {
        Parameter parameter = getParameter(parameterName);
        operators.put(key, new Operator(name, description, parameter, type, tuning, weight));
    }

    public void createOperator(String key, String name, String description, Parameter parameter1, Parameter parameter2, OperatorType type, double tuning, double weight) {
//        Parameter parameter1 = getParameter(parameterName1);
//        Parameter parameter2 = getParameter(parameterName2);
        operators.put(key, new Operator(name, description, parameter1, parameter2, type, tuning, weight));
    }

    public void createScaleOperator(String parameterName, double tuning, double weight) {
        Parameter parameter = getParameter(parameterName);
        operators.put(parameterName, new Operator(parameterName, "", parameter, OperatorType.SCALE, tuning, weight));
    }

    public void createScaleAllOperator(String parameterName, double weight) {
        Parameter parameter = getParameter(parameterName);
        operators.put(parameterName, new Operator(parameterName, "", parameter, OperatorType.SCALE_ALL, 0.75, weight));
    }

    public Parameter createParameter(String name, String description) {
        final Parameter parameter = new Parameter(name, description);
        parameters.put(name, parameter);
        return parameter;
    }

    public Parameter createParameter(String name, String description, int scale, double value, double lower, double upper) {
        final Parameter parameter = new Parameter(name, description, scale, value, lower, upper);
        parameters.put(name, parameter);
        return parameter;
    }

    public void createParameter(String name, String description, boolean isNodeHeight, double value, double lower, double upper) {
        parameters.put(name, new Parameter(name, description, isNodeHeight, value, lower, upper));
    }

    public void createScaleParameter(String name, String description, int scale, double value, double lower, double upper) {
        Parameter p = createParameter(name, description, scale, value, lower, upper);
        p.priorType = PriorType.JEFFREYS_PRIOR;
    }

    public Parameter createStatistic(String name, String description, boolean isDiscrete) {
        final Parameter parameter = new Parameter(name, description, isDiscrete);
        parameters.put(name, parameter);
        return parameter;
    }

    public void createStatistic(String name, String description, double lower, double upper) {
        parameters.put(name, new Parameter(name, description, lower, upper));
    }

    public Parameter getParameter(String name) {
        Parameter parameter = parameters.get(name);
        if (parameter == null) {
            for (String key : parameters.keySet()) {
                System.err.println(key);
            }
            throw new IllegalArgumentException("Parameter with name, " + name + ", is unknown");
        }
        return parameter;
    }

    public Parameter getStatistic(TaxonList taxonList) {
        Parameter parameter = statistics.get(taxonList);
        if (parameter == null) {
            for (TaxonList key : statistics.keySet()) {
                System.err.println("Taxon list: " + key.getId());
            }
            throw new IllegalArgumentException("Statistic for taxon list, " + taxonList.getId() + ", is unknown");
        }
        return parameter;
    }

    public Operator getOperator(String name) {
        Operator operator = operators.get(name);
        if (operator == null) throw new IllegalArgumentException("Operator with name, " + name + ", is unknown");
        return operator;
    }

    abstract public String getPrefix();

    protected void addComponent(ComponentOptions component) {
        components.add(component);
        component.createParameters(this);
    }

    public ComponentOptions getComponentOptions(Class<?> theClass) {
        for (ComponentOptions component : components) {
            if (theClass.isAssignableFrom(component.getClass())) {
                return component;
            }
        }

        return null;
    }

    protected void selectComponentParameters(ModelOptions options, List<Parameter> params) {
        for (ComponentOptions component : components) {
            component.selectParameters(options, params);
        }
    }

    protected void selectComponentStatistics(ModelOptions options, List<Parameter> stats) {
        for (ComponentOptions component : components) {
            component.selectStatistics(options, stats);
        }
    }

    protected void selectComponentOperators(ModelOptions options, List<Operator> ops) {
        for (ComponentOptions component : components) {
            component.selectOperators(options, ops);
        }
    }

    private final List<ComponentOptions> components = new ArrayList<ComponentOptions>();

}
