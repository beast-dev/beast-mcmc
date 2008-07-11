package dr.inference.operators;

import dr.inference.distribution.LinearRegression;
import dr.inference.distribution.MultivariateDistributionLikelihood;
import dr.inference.model.Parameter;
import dr.math.MultivariateDistribution;
import dr.math.MultivariateNormalDistribution;
import dr.math.matrixAlgebra.SymmetricMatrix;
import dr.xml.*;


/**
 * @author Marc Suchard
 */
public class RegressionGibbsOperator extends SimpleMCMCOperator implements GibbsOperator {

    public static final String GIBBS_OPERATOR = "regressionGibbsOperator";
    //	private TreeModel treeModel;
    //	private MatrixParameter precisionMatrixParameter;
    private double[][] effectPrecisionParameter;
    private LinearRegression linearModel;
    private Parameter effect;
    private Parameter outcome;
    private MultivariateDistribution effectPrior;
    private int dim;
    private int effectNumber;
    private int N;
    private int numEffects;
    private double[][] X;
//	private String traitName;

//	private int numberObservations;
//	private int weight;

    public RegressionGibbsOperator(LinearRegression linearModel, Parameter effect, MultivariateDistributionLikelihood effectPrior) {
        super();
        this.linearModel = linearModel;
        this.effect = effect;
        effectNumber = linearModel.getEffectNumber(effect);
        this.effectPrior = effectPrior.getDistribution();
        dim = effect.getDimension();
        outcome = linearModel.getDependentVariable();
        N = outcome.getDimension();
        numEffects = linearModel.getNumberOfEffects();
        X = linearModel.getX(effectNumber);

//		effectPrecisionParameter = effectPrior.getScaleMatrix();
//		numEffects = linearModel.getNumberOfEffects();
    }


    public int getStepCount() {
        return 1;
    }


    public double doOperation() throws OperatorFailedException {

//		int effectNum = MathUtils.nextInt(numEffects);
//	    Parameter effectParameter  = linearModel.getEffect(effectNum);
//		Parameter outcomeParameter = linearModel.getDependentVariable();

        double[] W = outcome.getParameterValues(); // outcome, fresh copy
        double[] P = linearModel.getScale();  // outcome precision, fresh copy
        double[] Beta = effect.getParameterValues(); // effect, fresh copy

//		final int N = outcome.getDimension();

        for (int k = 0; k < numEffects; k++) {
            if (k != effectNumber) {
                double[] thisXBeta = linearModel.getXBeta(k);
                for (int i = 0; i < N; i++)
                    W[i] -= thisXBeta[i];

            }
        }

//		double[][] X = linearModel.getX(effectNumber); // pointer only
//		final int dim = X[0].length;

        double[] priorBetaMean = effectPrior.getMean();
        double[][] priorBetaScale = effectPrior.getScaleMatrix();

        double[][] XtP = new double[dim][N];
        for (int j = 0; j < dim; j++) {
            for (int i = 0; i < N; i++)
                XtP[j][i] = X[i][j] * P[i];
        }

        double[][] XtPX = new double[dim][dim];
        for (int i = 0; i < dim; i++) {
            for (int j = i; j < dim; j++) {// symmetric
                for (int k = 0; k < N; k++)
                    XtPX[i][j] += XtP[i][k] * X[k][j];
                XtPX[j][i] = XtPX[i][j]; // symmetric
            }
        }

        double[][] XtPX_plus_P0 = new double[dim][dim];
        for (int i = 0; i < dim; i++) {
            for (int j = i; j < dim; j++) // symmetric
                XtPX_plus_P0[j][i] = XtPX_plus_P0[i][j] = XtPX[i][j] + priorBetaScale[i][j];
        }

        double[] XtPW = new double[dim];
        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < N; j++)
                XtPW[i] += XtP[i][j] * W[j];
        }

        double[] P0Mean0 = new double[dim];
        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++)
                P0Mean0[i] += priorBetaScale[i][j] * priorBetaMean[j];
        }

        double[] unscaledMean = new double[dim];
        for (int i = 0; i < dim; i++)
            unscaledMean[i] = P0Mean0[i] + XtPW[i];

        double[][] variance = new SymmetricMatrix(XtPX_plus_P0).inverse().toComponents();

        double[] scaledMean = new double[dim];
        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++)
                scaledMean[i] += variance[i][j] * unscaledMean[j];
        }

//		System.err.println("Mean:\n"+new Vector(scaledMean));
//		System.err.println("Var :\n"+new Matrix(variance));

        double[] draw = MultivariateNormalDistribution.nextMultivariateNormalVariance(
                scaledMean, variance);

        for (int i = 0; i < dim; i++)
            effect.setParameterValue(i, draw[i]);

        return 0;
    }

    public String getPerformanceSuggestion() {
        return null;
    }

    public String getOperatorName() {
        return GIBBS_OPERATOR;
    }

    public static dr.xml.XMLObjectParser PARSER = new dr.xml.AbstractXMLObjectParser() {

        public String getParserName() {
            return GIBBS_OPERATOR;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            double weight = xo.getDoubleAttribute(WEIGHT);

            LinearRegression linearModel = (LinearRegression) xo.getChild(LinearRegression.class);
            Parameter effect = (Parameter) xo.getChild(Parameter.class);
            MultivariateDistributionLikelihood prior = (MultivariateDistributionLikelihood) xo.getChild(MultivariateDistributionLikelihood.class);
//			System.err.println("prior: "+prior.getId());
            if (prior.getDistribution().getType().compareTo(MultivariateNormalDistribution.TYPE) != 0)
                throw new XMLParseException("Only a multivariate normal prior is conjugate");
//			TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);

//			MultivariateDiffusionModel diffusionModel = (MultivariateDiffusionModel) xo.getChild(MultivariateDiffusionModel.class);

            RegressionGibbsOperator operator = new RegressionGibbsOperator(linearModel, effect, prior);
            operator.setWeight(weight);
            return operator;
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "This element returns a multivariate Gibbs operator on an internal node trait.";
        }

        public Class getReturnType() {
            return MCMCOperator.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                AttributeRule.newDoubleRule(WEIGHT),
//				new ElementRule(TreeModel.class),
                new ElementRule(Parameter.class),
//				new ElementRule(MultivariateDiffusionModel.class)
                new ElementRule(MultivariateDistributionLikelihood.class),
                new ElementRule(LinearRegression.class)
        };

    };

}
