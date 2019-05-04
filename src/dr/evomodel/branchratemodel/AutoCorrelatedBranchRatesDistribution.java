package dr.evomodel.branchratemodel;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.tree.TreeModel;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.*;
import dr.math.distributions.NormalDistribution;
import dr.util.Citable;
import dr.util.Citation;

import java.util.List;

/**
 * @author Marc A. Suchard
 * @author Xiang Ji
 * @author Philippe Lemey
 */
public class AutoCorrelatedBranchRatesDistribution extends AbstractModelLikelihood
        implements GradientWrtParameterProvider, Citable {

    private final ArbitraryBranchRates branchRateModel;
    private final Parameter precisionParameter;
    private final BranchVarianceScaling scaling;

    private final Tree tree;
    private final Parameter rateParameter;

    private boolean likelihoodKnown = false;
    private boolean savedLikelihoodKnown;

    private double logLikelihood;
    private double savedLogLikelihood;

    private final int dim;
    private double[] increments;
    private double[] savedIncrements;

    private static final boolean TEST_INCREMENTS = true;

    public AutoCorrelatedBranchRatesDistribution(String name,
                                                 ArbitraryBranchRates branchRateModel,
                                                 Parameter precision,
                                                 BranchVarianceScaling scaling) {
        super(name);
        this.branchRateModel = branchRateModel;
        this.precisionParameter = precision;
        this.scaling = scaling;

        this.tree = branchRateModel.getTree();
        this.rateParameter = branchRateModel.getRateParameter();

        addModel(branchRateModel);
        addVariable(precision);

        if (tree instanceof TreeModel) {
            addModel((TreeModel) tree);
        }

        this.dim = branchRateModel.getRateParameter().getDimension();
        this.increments = new double[dim];
        this.savedIncrements = new double[dim];
    }

    @Override
    public Likelihood getLikelihood() {
        return this;
    }

    @Override
    public Parameter getParameter() {
        return rateParameter;
    }

    @Override
    public int getDimension() {
        return rateParameter.getDimension();
    }

    @Override
    public double[] getGradientLogDensity() {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        likelihoodKnown = false;
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        likelihoodKnown = false;
    }

    @Override
    protected void storeState() {
        savedLikelihoodKnown = likelihoodKnown;
        savedLogLikelihood = logLikelihood;

        System.arraycopy(increments, 0, savedIncrements, 0, dim);
    }

    @Override
    protected void restoreState() {
        likelihoodKnown = savedLikelihoodKnown;
        logLikelihood = savedLogLikelihood;

        double[] tmp = savedIncrements;
        savedIncrements = increments;
        increments = tmp;
    }

    @Override
    protected void acceptState() {

    }

    @Override
    public Model getModel() {
        return this;
    }

    @Override
    public double getLogLikelihood() {
        if (!likelihoodKnown) {
            logLikelihood = calculateLogLikelihood();
        }
        return logLikelihood;
    }

    @Override
    public void makeDirty() {
        likelihoodKnown = false;
    }

    @Override
    public Citation.Category getCategory() {
        return null;
    }

    @Override
    public String getDescription() {
        return null;
    }

    @Override
    public List<Citation> getCitations() {
        return null;
    }

    private double calculateLogLikelihood() {
        double assumedUntransformedRootRate = 0.0;
        NodeRef root = tree.getRoot();
        double preOrderComputation = recursePreOrder(tree.getChild(root, 0), assumedUntransformedRootRate)
                + recursePreOrder(tree.getChild(root, 1), assumedUntransformedRootRate);

        if (TEST_INCREMENTS) {
            double logLike = 0;
            for (int i = 0; i < dim; ++i) {
                logLike += getLogLikelihoodOfBranch(increments[i]);
            }

            System.err.println("pre: " + preOrderComputation + " inc: " + logLike);
        }

        return preOrderComputation;
    }

    private double recursePreOrder(NodeRef node, double untransformedParentRate) {

        final double untransformedRate = branchRateModel.getUntransformedBranchRate(tree, node);
        final double branchLength = tree.getBranchLength(node);
        final double rateIncrement = scaling.rescaleIncrement(
                untransformedRate - untransformedParentRate, branchLength);

        increments[branchRateModel.getParameterIndexFromNode(node)] = rateIncrement;

        double likelihoodContribution = getLogLikelihoodOfBranch(rateIncrement);

        if (!tree.isExternal(node)) {
            likelihoodContribution +=
                    recursePreOrder(tree.getChild(node, 0), untransformedRate) +
                            recursePreOrder(tree.getChild(node, 1), untransformedRate);
        }

        return likelihoodContribution;
    }

    private double getLogLikelihoodOfBranch(double rateIncrement) {
        final double sd = 1.0 / Math.sqrt(precisionParameter.getParameterValue(0));
        return NormalDistribution.logPdf(rateIncrement, 0.0, sd);
    }

    public enum BranchVarianceScaling {

        NONE("none") {
            @Override
            double rescaleIncrement(double increment, double branchLength) {
                return increment;
            }
        },

        BY_TIME("byTime") {
            @Override
            double rescaleIncrement(double increment, double branchLength) {
                return increment / Math.sqrt(branchLength);
            }
        };

        BranchVarianceScaling(String name) {
            this.name = name;
        }

        private final String name;

        abstract double rescaleIncrement(double increment, double branchLength);

        public String getName() {
            return name;
        }

        public static BranchVarianceScaling parse(String name) {
            for (BranchVarianceScaling scaling : BranchVarianceScaling.values()) {
                if (scaling.getName().equalsIgnoreCase(name)) {
                    return scaling;
                }
            }
            return null;
        }
    }
}
