package dr.evomodel.operators;

import dr.evomodel.tree.TreeModel;
import dr.evomodel.continuous.AbstractMultivariateTraitLikelihood;
import dr.evomodel.continuous.IntegratedMultivariateTraitLikelihood;
import dr.evomodel.branchratemodel.ArbitraryBranchRates;
import dr.inference.model.MatrixParameter;
import dr.inference.operators.OperatorFailedException;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.SimpleMCMCOperator;
import dr.inference.operators.GibbsOperator;
import dr.inference.distribution.GammaDistributionModel;
import dr.inference.distribution.DistributionLikelihood;
import dr.evolution.tree.NodeRef;
import dr.math.distributions.GammaDistribution;
import dr.xml.*;

import java.util.logging.Logger;

/**
 * Gibbs samples each of AbritraryBranchRates when their prior is a gamma distribution
 *
 * @author Marc A. Suchard
 */
public class TraitRateGibbsOperator extends SimpleMCMCOperator implements GibbsOperator {

    public static final String GIBBS_OPERATOR = "traitRateGibbsOperator";

    private final TreeModel treeModel;
    private final MatrixParameter precisionMatrixParameter;
    private final AbstractMultivariateTraitLikelihood traitModel;
    private final GammaDistributionModel ratePriorModel;
    private final GammaDistribution ratePrior;
    private final ArbitraryBranchRates branchRateModel;
    private final int dim;
    private final String traitName;

    public TraitRateGibbsOperator(AbstractMultivariateTraitLikelihood traitModel,
                                  ArbitraryBranchRates branchRateModel,
                                  GammaDistributionModel ratePriorModel,
                                  GammaDistribution ratePrior) {
        super();
        this.traitModel = traitModel;
        this.treeModel = traitModel.getTreeModel();
        this.precisionMatrixParameter = (MatrixParameter) traitModel.getDiffusionModel().getPrecisionParameter();
        this.traitName = traitModel.getTraitName();
        this.branchRateModel = branchRateModel;
        this.ratePriorModel = ratePriorModel;
        this.ratePrior = ratePrior;
        this.dim = treeModel.getMultivariateNodeTrait(treeModel.getRoot(), traitName).length;

        boolean hasDistributionModel = ratePriorModel == null;
        boolean hasDistribution = ratePrior == null;

        if (traitModel instanceof IntegratedMultivariateTraitLikelihood) {
            throw new RuntimeException("Only implemented for a SampledMultivariateTraitLikelihood");
        }

        if ((hasDistribution && hasDistributionModel) || (!hasDistribution && !hasDistributionModel)) {
            throw new RuntimeException("Can only provide one prior density in TraitRateGibbsOperation");
        }

        if (!branchRateModel.usingReciprocal()) {
            throw new RuntimeException("ArbitraryBranchRates in TraitRateGibbsOperatior must use reciprocal rates");
        }
                
        Logger.getLogger("dr.evomodel").info("Using Gibbs operator and trait rates");
    }

    public int getStepCount() {
        return 1;
    }

    private void sampleRateForNode(NodeRef child, double[][] precision, double priorShape, double priorRate) {
        
        NodeRef parent = treeModel.getParent(child);

        final double[] trait       = treeModel.getMultivariateNodeTrait(child,  traitName);
        final double[] parentTrait = treeModel.getMultivariateNodeTrait(parent, traitName);

        final double precisionScalar = branchRateModel.getBranchRate(treeModel, child) /
                                       traitModel.getRescaledBranchLength(child);

        for (int i = 0; i < dim; i++) {
            trait[i] -= parentTrait[i];
        }

        double SSE = 0;

        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
                SSE += trait[i] * precision[i][j] * trait[j];
            }
        }

        final double gammaShape = priorShape + 0.5 * dim;
        final double gammaRate  = priorRate  + 0.5 * SSE * precisionScalar;

        final double newValue = GammaDistribution.nextGamma(gammaShape, 1.0 / gammaRate);

        // Store the reciprocal value as the rate (\propto variance)
        branchRateModel.setBranchRate(treeModel, child, 1.0 / newValue);
    }

    public double doOperation() throws OperatorFailedException {

        double[][] precision = precisionMatrixParameter.getParameterAsMatrix();

        double priorShape;
        double priorRate;

        if (ratePriorModel != null) {
            priorShape = ratePriorModel.getShape();
            priorRate  = 1.0 / ratePriorModel.getScale();
        } else {
            priorShape = ratePrior.getShape();
            priorRate = 1.0 / ratePrior.getScale();
        }

        for (int i = 0; i < treeModel.getNodeCount(); i++) {
            NodeRef node = treeModel.getNode(i);
            if (node != treeModel.getRoot()) {
                sampleRateForNode(node, precision, priorShape, priorRate);
            }
        }
        return 0;
    }

    public String getPerformanceSuggestion() {
        return null;
    }

    public String getOperatorName() {
        return GIBBS_OPERATOR;
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return GIBBS_OPERATOR;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            double weight = xo.getDoubleAttribute(WEIGHT);
            AbstractMultivariateTraitLikelihood traitLikelihood =
                    (AbstractMultivariateTraitLikelihood) xo.getChild(AbstractMultivariateTraitLikelihood.class);
            ArbitraryBranchRates branchRates = (ArbitraryBranchRates) xo.getChild(ArbitraryBranchRates.class);

            DistributionLikelihood priorLikelihood = (DistributionLikelihood) xo.getChild(DistributionLikelihood.class);

            GammaDistributionModel gammaPriorModel = null;
            GammaDistribution gammaPrior = null;

            if (priorLikelihood.getDistribution() instanceof GammaDistributionModel) {
                gammaPriorModel = (GammaDistributionModel) priorLikelihood.getDistribution();
            } else if (priorLikelihood.getDistribution() instanceof GammaDistribution) {
                gammaPrior = (GammaDistribution) priorLikelihood.getDistribution();
            } else {
                throw new XMLParseException("Currently only works with a GammaDistributionModel or GammaDistribution");
            }

            if (!branchRates.usingReciprocal()) {
                throw new XMLParseException(
                        "Gibbs sampling of rates only works with reciprocal rates under an ArbitraryBranchRates model");
            }

            TraitRateGibbsOperator operator = new TraitRateGibbsOperator(traitLikelihood, branchRates, gammaPriorModel, gammaPrior);
            operator.setWeight(weight);

            return operator;
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "This element returns a multivariate Gibbs operator on traits for possible all nodes.";
        }

        public Class getReturnType() {
            return MCMCOperator.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                AttributeRule.newDoubleRule(WEIGHT),
                new ElementRule(AbstractMultivariateTraitLikelihood.class),
                new ElementRule(ArbitraryBranchRates.class),
                new ElementRule(DistributionLikelihood.class),
        };

    };
}
