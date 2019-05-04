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
    }

    @Override
    protected void restoreState() {
        likelihoodKnown = savedLikelihoodKnown;
        logLikelihood = savedLogLikelihood;
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
        return recursePreOrder(tree.getChild(root, 0), assumedUntransformedRootRate)
                + recursePreOrder(tree.getChild(root, 1), assumedUntransformedRootRate);
    }

    private double recursePreOrder(NodeRef node, double untransformedParentRate) {

        double untransformedRate = branchRateModel.getUntransformedBranchRate(tree, node);
        double likelihoodContribution = getLogLikelihoodOfBranch(untransformedParentRate, untransformedRate,
                tree.getBranchLength(node));

        if (!tree.isExternal(node)) {
            likelihoodContribution +=
                    recursePreOrder(tree.getChild(node, 0), untransformedRate) +
                            recursePreOrder(tree.getChild(node, 1), untransformedRate);
        }

        return likelihoodContribution;
    }

    private double getLogLikelihoodOfBranch(double parentValue, double value, double branchLength) {
        double sd = scaling.getSD(precisionParameter.getParameterValue(0), branchLength);
        return NormalDistribution.logPdf(value, parentValue, sd);
    }

    enum BranchVarianceScaling {

        NONE {
            @Override
            double getSD(double precision, double branchLength) {
                return Math.sqrt(1.0 / precision);
            }
        },

        BY_TIME {
            @Override
            double getSD(double precision, double branchLength) {
                return Math.sqrt(branchLength / precision);
            }
        };

        abstract double getSD(double precision, double branchLength);
    }
}
