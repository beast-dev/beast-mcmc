package dr.app.beauti.options;

import dr.app.beauti.priorsPanel.PriorType;
import dr.evolution.tree.Tree;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Andrew Rambaut
 * @author Walter Xie
 * @version $Id$
 */
public class PartitionTreeModel extends ModelOptions {

    // Instance variables

    private final BeautiOptions options;
    private String name;
    private PartitionTreePrior treePrior;
    private List<PartitionData> allPartitionData;

    public Parameter localClockRateChangesStatistic = null;
    public Parameter localClockRatesStatistic = null;

    private StartingTreeType startingTreeType = StartingTreeType.RANDOM;
    private Tree userStartingTree = null;

    private boolean fixedTree = false;

    public PartitionTreeModel(BeautiOptions options, PartitionData partition) {
        this.options = options;
        this.name = partition.getName();

        allPartitionData = new ArrayList<PartitionData>();
        addPartitionData(partition);

        initTreeModelParaAndOpers();
    }

    /**
     * A copy constructor
     *
     * @param options the beauti options
     * @param name    the name of the new model
     * @param source  the source model
     */
    public PartitionTreeModel(BeautiOptions options, String name, PartitionTreeModel source) {
        this.options = options;
        this.name = name;

        this.allPartitionData = source.allPartitionData;

        this.startingTreeType = source.startingTreeType;
        this.userStartingTree = source.userStartingTree;

        initTreeModelParaAndOpers();
    }

//    public PartitionTreeModel(BeautiOptions options, String name) {
//        this.options = options;
//        this.name = name;
//    }

    public void initTreeModelParaAndOpers() {
        double branchWeights = 30.0;
        double treeWeights = 15.0;
        double rateWeights = 3.0;

        createParameter("tree", "The tree");
        createParameter("treeModel.internalNodeHeights", "internal node heights of the tree (except the root)");
        createParameter("treeModel.allInternalNodeHeights", "internal node heights of the tree");
        createParameter("treeModel.rootHeight", "root height of the tree", true, 1.0, 0.0, Double.POSITIVE_INFINITY);

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
                "Scales all rates inversely to node heights of the tree", super.getParameter("treeModel.allRates"),
                super.getParameter("treeModel.allInternalNodeHeights"), OperatorType.UP_DOWN, 0.75, branchWeights);

        createOperator("swapBranchRateCategories", "branchRates.categories",
                "Performs a swap of branch rate categories", "branchRates.categories",
                OperatorType.SWAP, 1, branchWeights / 3);
        createOperator("randomWalkBranchRateCategories", "branchRates.categories",
                "Performs an integer random walk of branch rate categories", "branchRates.categories",
                OperatorType.INTEGER_RANDOM_WALK, 1, branchWeights / 3);
        createOperator("unformBranchRateCategories", "branchRates.categories",
                "Performs an integer uniform draw of branch rate categories", "branchRates.categories",
                OperatorType.INTEGER_UNIFORM, 1, branchWeights / 3);

        createScaleOperator(ClockType.LOCAL_CLOCK + "." + "rates", demoTuning, treeWeights);
        createOperator(ClockType.LOCAL_CLOCK + "." + "changes", OperatorType.BITFLIP, 1, treeWeights);
        createOperator("treeBitMove", "Tree", "Swaps the rates and change locations of local clocks", "tree",
                OperatorType.TREE_BIT_MOVE, -1.0, treeWeights);

        createScaleOperator("treeModel.rootHeight", demoTuning, demoWeights);
        createOperator("uniformHeights", "Internal node heights", "Draws new internal node heights uniformally",
                "treeModel.internalNodeHeights", OperatorType.UNIFORM, -1, branchWeights);

        createOperator("subtreeSlide", "Tree", "Performs the subtree-slide rearrangement of the tree", "tree",
                OperatorType.SUBTREE_SLIDE, 1.0, treeWeights);
        createOperator("narrowExchange", "Tree", "Performs local rearrangements of the tree", "tree",
                OperatorType.NARROW_EXCHANGE, -1, treeWeights);
        createOperator("wideExchange", "Tree", "Performs global rearrangements of the tree", "tree",
                OperatorType.WIDE_EXCHANGE, -1, demoWeights);
        createOperator("wilsonBalding", "Tree", "Performs the Wilson-Balding rearrangement of the tree", "tree",
                OperatorType.WILSON_BALDING, -1, demoWeights);
    }

    /**
     * return a list of parameters that are required
     *
     * @param params the parameter list
     */
    public void selectParameters(List<Parameter> params) {

        params.add(getParameter("treeModel.rootHeight"));

        selectStatistics(params);
    }

    /**
     * return a list of operators that are required
     *
     * @param ops the operator list
     */
    public void selectOperators(List<Operator> ops) {

        // if not a fixed tree then sample tree space
        if (!fixedTree) {
            ops.add(getOperator("subtreeSlide"));
            ops.add(getOperator("narrowExchange"));
            ops.add(getOperator("wideExchange"));
            ops.add(getOperator("wilsonBalding"));
        }

    }

    // use override method getParameter(String name) and getOperator(String name) in PartitionModel containing prefix
    private void selectStatistics(List<Parameter> params) {

//        if (options.taxonSets != null) {
//            for (Taxa taxonSet : options.taxonSets) {
//                Parameter statistic = statistics.get(taxonSet);
//                if (statistic == null) {
//                    statistic = new Parameter(taxonSet, "tMRCA for taxon set ");
//                    statistics.put(taxonSet, statistic);
//                }
//                params.add(statistic);
//            }
//        } else {
//            System.err.println("TaxonSets are null");
//        }
        //TODO ?
        for (PartitionClockModel clock : options.getPartitionClockModels(getAllPartitionData())) {
            if (clock.getClockType() == ClockType.RANDOM_LOCAL_CLOCK) {
                if (localClockRateChangesStatistic == null) {
                    localClockRateChangesStatistic = new Parameter("rateChanges", "number of random local clocks", true);
                    localClockRateChangesStatistic.priorType = PriorType.POISSON_PRIOR;
                    localClockRateChangesStatistic.poissonMean = 1.0;
                    localClockRateChangesStatistic.poissonOffset = 0.0;
                }
                if (localClockRatesStatistic == null) {
                    localClockRatesStatistic = new Parameter(ClockType.LOCAL_CLOCK + "." + "rates", "random local clock rates", false);

                    localClockRatesStatistic.priorType = PriorType.GAMMA_PRIOR;
                    localClockRatesStatistic.gammaAlpha = 0.5;
                    localClockRatesStatistic.gammaBeta = 2.0;
                }

                localClockRateChangesStatistic.setPrefix(getPrefix());
                params.add(localClockRateChangesStatistic);
                localClockRatesStatistic.setPrefix(getPrefix());
                params.add(localClockRatesStatistic);
            }

//	        if (clock.getClockType() != ClockType.STRICT_CLOCK) {
//	            params.add(getParameter("meanRate"));
//	            params.add(getParameter(RateStatistic.COEFFICIENT_OF_VARIATION));
//	            params.add(getParameter("covariance"));
//	        }
        }
    }

    /////////////////////////////////////////////////////////////

    public List<PartitionData> getAllPartitionData() {
        return allPartitionData;
    }

    public void clearAllPartitionData() {
        this.allPartitionData.clear();
    }

    public void addPartitionData(PartitionData partition) {
        allPartitionData.add(partition);
    }

    public boolean removePartitionData(PartitionData partition) {
        return allPartitionData.remove(partition);
    }

    public PartitionTreePrior getPartitionTreePrior() {
        return treePrior;
    }

    public void setPartitionTreePrior(PartitionTreePrior treePrior) {
        this.treePrior = treePrior;
    }

    public StartingTreeType getStartingTreeType() {
        return startingTreeType;
    }

    public void setStartingTreeType(StartingTreeType startingTreeType) {
        this.startingTreeType = startingTreeType;
    }

    public Tree getUserStartingTree() {
        return userStartingTree;
    }

    public void setUserStartingTree(Tree userStartingTree) {
        this.userStartingTree = userStartingTree;
    }

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }

    public String toString() {
        return getName();
    }

    public Parameter getParameter(String name) {

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

        operator.setPrefix(getPrefix());

        return operator;
    }


    public String getPrefix() {
        String prefix = "";
        if (options.getPartitionTreeModels().size() > 1) { //|| options.isSpeciesAnalysis()
            // There is more than one active partition model
            prefix += getName() + ".";
        }
        return prefix;
    }


}
