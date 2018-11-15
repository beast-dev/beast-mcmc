package dr.evomodel.treedatalikelihood.continuous;

import dr.evolution.tree.TreeTrait;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.inference.distribution.MultivariateDistributionLikelihood;
import dr.inference.model.CompoundParameter;
import dr.inference.model.MatrixParameterInterface;
import dr.inference.model.Parameter;
import dr.inference.operators.GibbsOperator;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.distributions.WishartDistribution;
import dr.math.distributions.WishartStatistics;
import dr.math.interfaces.ConjugateWishartStatisticsProvider;
import dr.xml.*;
import dr.math.matrixAlgebra.Matrix;

import static dr.evomodel.treedatalikelihood.preorder.AbstractRealizedContinuousTraitDelegate.REALIZED_TIP_TRAIT;

public class RepeatedMeasuresPrecisionGibbsOperator extends SimpleMCMCOperator implements GibbsOperator {

    public static final String OPERATOR_NAME = "repeatedMeasuresPrecisionGibbsOperator";

    private static WishartStatistics priorDistribution;
    private RepeatedMeasuresTraitDataModel dataModel;
    private TreeDataLikelihood dataLikelihood;


    public RepeatedMeasuresPrecisionGibbsOperator(
            RepeatedMeasuresTraitDataModel dataModel,
            TreeDataLikelihood dataLikelihood,
            WishartStatistics priorDistribution,
            double weight) {

        this.priorDistribution = priorDistribution;
        this.dataModel = dataModel;
        this.dataLikelihood = dataLikelihood;


    }


    @Override
    public String getPerformanceSuggestion() {
        return null;
    }

    @Override
    public String getOperatorName() {
        return OPERATOR_NAME;
    }

    @Override
    public double doOperation() {
        int dim = dataModel.dimTrait;
        double df = priorDistribution.getDF() + dataLikelihood.getTree().getExternalNodeCount();
        //TODO: Draw missing data values
        double[][] scaleMatrix = getPosteriorScaleMatrix();

        double[][] draw = WishartDistribution.nextWishart(df, scaleMatrix);

        for (int i = 0; i < dim; i++) {

            Parameter column = dataModel.getPrecisionMatrix().getParameter(i);

            for (int j = 0; j < dim; j++)

                column.setParameterValueQuietly(j, draw[j][i]);
        }

        dataModel.getPrecisionMatrix().fireParameterChangedEvent();

        return 0;
    }

    private double[][] getPosteriorScaleMatrix() {

        int dimTrait = dataModel.dimTrait;
        int N = dataLikelihood.getTree().getExternalNodeCount();
        double[][] scaleMatrix = new Matrix(priorDistribution.getScaleMatrix()).inverse().toComponents();

        CompoundParameter traitParameter = dataModel.getParameter();

        TreeTrait tipTrait = dataLikelihood.getTreeTrait(REALIZED_TIP_TRAIT + "." + dataModel.getTraitName());
        double[] tipTraits = (double[]) tipTrait.getTrait(dataLikelihood.getTree(), null);

        for (int i = 0; i < dimTrait; i++) {
            for (int j = 0; j < dimTrait; j++) {

                double value = 0;

                for (int k = 0; k < N; k++) {

                    Parameter taxonParameter = traitParameter.getParameter(k);

                    value = value + (taxonParameter.getParameterValue(i) - tipTraits[k * dimTrait + i]) *
                            (taxonParameter.getParameterValue(j) - tipTraits[k * dimTrait + j]);
                }

                scaleMatrix[i][j] += value;
            }
        }

        scaleMatrix = new Matrix(scaleMatrix).inverse().toComponents();

        return scaleMatrix;
    }


    private static final String RM_GIBBS_OPERATOR = "repeatedMeasuresPrecisionGibbsOperator";

    public static AbstractXMLObjectParser PARSER = new AbstractXMLObjectParser() {
        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            RepeatedMeasuresTraitDataModel dataModel = (RepeatedMeasuresTraitDataModel) xo.getChild(RepeatedMeasuresTraitDataModel.class);
            TreeDataLikelihood dataLikelihood = (TreeDataLikelihood) xo.getChild(TreeDataLikelihood.class);
//            ConjugateWishartStatisticsProvider ws = (ConjugateWishartStatisticsProvider) xo.getChild(ConjugateWishartStatisticsProvider.class);
            MultivariateDistributionLikelihood prior = (MultivariateDistributionLikelihood) xo.getChild(MultivariateDistributionLikelihood.class);
            assert (prior.getDistribution() instanceof WishartDistribution);

            double weight = xo.getDoubleAttribute(WEIGHT);


            return new RepeatedMeasuresPrecisionGibbsOperator(dataModel, dataLikelihood, (WishartStatistics) prior.getDistribution(),
                    weight);
        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                AttributeRule.newDoubleRule(WEIGHT),
                new ElementRule(MultivariateDistributionLikelihood.class),
                new ElementRule(RepeatedMeasuresTraitDataModel.class),
                new ElementRule(TreeDataLikelihood.class)
        };

        @Override
        public String getParserDescription() {
            return null;
        }

        @Override
        public Class getReturnType() {
            return null;
        }

        @Override
        public String getParserName() {
            return RM_GIBBS_OPERATOR;
        }
    };


}
