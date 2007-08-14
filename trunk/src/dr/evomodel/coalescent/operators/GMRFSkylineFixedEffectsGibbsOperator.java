package dr.evomodel.coalescent.operators;

import dr.evomodel.coalescent.GMRFSkylineLikelihood;
import dr.inference.distribution.MultivariateDistributionLikelihood;
import dr.inference.model.Parameter;
import dr.inference.operators.GibbsOperator;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.OperatorFailedException;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.MultivariateDistribution;
import dr.math.MultivariateNormalDistribution;
import dr.xml.*;
import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.DenseVector;

/**
 * A Gibbs operator to update the population size parameters under a Gaussian Markov random field prior
 *
 * @author Erik Bloomquist
 * @author Vladimir Minin
 * @author Marc Suchard
 * @version $Id: GMRFSkylineFixedEffectsGibbsOperator.java,v 1.5 2007/03/20 11:26:49 msuchard Exp $
 */

public class GMRFSkylineFixedEffectsGibbsOperator extends SimpleMCMCOperator implements GibbsOperator {


    public static final String GMRF_GIBBS_OPERATOR = "gmrfFixedEffectsGibbsOperator";

    private Parameter param;
    private GMRFSkylineLikelihood gmrfLikelihood;
    private MultivariateDistribution prior;

    private DenseVector mean;
    private DenseMatrix precision;

    private int weight = 1;

    private int fieldLength;
    private int dim;

    public GMRFSkylineFixedEffectsGibbsOperator(Parameter param,
                                                GMRFSkylineLikelihood gmrfLikelihood, MultivariateDistribution prior, int weight) {
        this.param = param;
        this.gmrfLikelihood = gmrfLikelihood;
        this.prior = prior;
        mean = new DenseVector(prior.getMean());
        precision = new DenseMatrix(prior.getScaleMatrix());

        this.fieldLength = gmrfLikelihood.getPopSizeParameter().getDimension();
        this.dim = param.getDimension();

        this.gmrfLikelihood = gmrfLikelihood;
        this.weight = weight;
    }

    public double doOperation() throws OperatorFailedException {

        // todo Erik to fill-in operation

        return 0;
    }

    public int getStepCount() {
        return 0;
    }


    public String getPerformanceSuggestion() {
        return null;
    }

    //MCMCOperator INTERFACE

    public final String getOperatorName() {
        return GMRF_GIBBS_OPERATOR;
    }


    public static dr.xml.XMLObjectParser PARSER = new dr.xml.AbstractXMLObjectParser() {

        public String getParserName() {
            return GMRF_GIBBS_OPERATOR;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {


            int weight = xo.getIntegerAttribute(WEIGHT);

            GMRFSkylineLikelihood gmrfLikelihood = (GMRFSkylineLikelihood) xo.getChild(GMRFSkylineLikelihood.class);

            MultivariateDistributionLikelihood likelihood = (MultivariateDistributionLikelihood) xo.getChild(MultivariateDistributionLikelihood.class);

            MultivariateDistribution prior = likelihood.getDistribution();

            if (prior.getType().compareTo(MultivariateNormalDistribution.TYPE) != 0)
                throw new XMLParseException("Only a multivariate normal distribution is conjugate for the regression coefficients in a GMRF");

            Parameter param = (Parameter) xo.getChild(Parameter.class);

            return new GMRFSkylineFixedEffectsGibbsOperator(param,
                    gmrfLikelihood, prior, weight);

        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "This element returns a Gibbs operator for regression coefficients in a GMRF.";
        }

        public Class getReturnType() {
            return MCMCOperator.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
//            return null;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                AttributeRule.newIntegerRule(WEIGHT),
                new ElementRule(MultivariateDistributionLikelihood.class),
                new ElementRule(Parameter.class),
                new ElementRule(GMRFSkylineLikelihood.class)
        };

    };

}
