package dr.evomodel.branchratemodel;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.tree.TreeModel;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.*;
import dr.math.distributions.NormalDistribution;
import dr.util.Author;
import dr.util.Citable;
import dr.util.Citation;

import java.util.Collections;
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

    private boolean incrementsKnown = false;
    private boolean savedIncrementsKnown;

    private boolean likelihoodKnown = false;
    private boolean savedLikelihoodKnown;

    private double logLikelihood;
    private double savedLogLikelihood;

    private final int dim;
    private double[] increments;
    private double[] savedIncrements;

    private Parameter incrementsProxy;

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
        incrementsKnown = false;
        likelihoodKnown = false;
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        likelihoodKnown = false;
    }

    @Override
    protected void storeState() {
        savedIncrementsKnown = incrementsKnown;
        System.arraycopy(increments, 0, savedIncrements, 0, dim);

        savedLikelihoodKnown = likelihoodKnown;
        savedLogLikelihood = logLikelihood;
    }

    @Override
    protected void restoreState() {
        incrementsKnown = savedIncrementsKnown;
        double[] tmp = savedIncrements;
        savedIncrements = increments;
        increments = tmp;

        likelihoodKnown = savedLikelihoodKnown;
        logLikelihood = savedLogLikelihood;
    }

    @Override
    protected void acceptState() {

    }

    public double getIncrement(int index) {
        checkIncrements();
        return increments[index];
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
        return Collections.singletonList(CITATION);
    }

    public static Citation CITATION = new Citation(
            new Author[]{
            },
            Citation.Status.IN_PREPARATION
    );

    private void checkIncrements() {
        if (!incrementsKnown) {
            recursePreOrder(tree.getRoot(), 0.0);
            incrementsKnown = true;
        }
    }

    private double calculateLogLikelihood() {
        checkIncrements();

        double logLike = 0;
        for (int i = 0; i < dim; ++i) {
            logLike += getLogLikelihoodOfBranch(increments[i]);
        }

        return logLike;
    }

    private void recursePreOrder(NodeRef node, double untransformedParentRate) {

        if (!tree.isRoot(node)) {
            final double untransformedRate = branchRateModel.getUntransformedBranchRate(tree, node);
            final double branchLength = tree.getBranchLength(node);
            final double rateIncrement = scaling.rescaleIncrement(
                    untransformedRate - untransformedParentRate, branchLength);

            increments[branchRateModel.getParameterIndexFromNode(node)] = rateIncrement;

            untransformedParentRate = untransformedRate;
        }

        if (!tree.isExternal(node)) {
            recursePreOrder(tree.getChild(node, 0), untransformedParentRate);
            recursePreOrder(tree.getChild(node, 1), untransformedParentRate);
        }
    }

    private double getLogLikelihoodOfBranch(double rateIncrement) {
        // TODO Delegate to a DistributionLikelihood using proxy
        final double sd = 1.0 / Math.sqrt(precisionParameter.getParameterValue(0));
        return NormalDistribution.logPdf(rateIncrement, 0.0, sd);
    }

    Parameter getIncrementsAsParameter() {

        if (incrementsProxy == null) {

            incrementsProxy = new Parameter.Proxy("incrementsProxy", increments.length) {
                @Override
                public double getParameterValue(int dim) {
                    checkIncrements();
                    return increments[dim];
                }

                @Override
                public void setParameterValue(int dim, double value) {
                    throw new RuntimeException("Cannot set increments directly");
                }

                @Override
                public void setParameterValueQuietly(int dim, double value) {
                    throw new RuntimeException("Cannot set increments directly");
                }

                @Override
                public void setParameterValueNotifyChangedAll(int dim, double value) {
                    throw new RuntimeException("Cannot set increments directly");
                }
            };
        }

        return incrementsProxy;
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
