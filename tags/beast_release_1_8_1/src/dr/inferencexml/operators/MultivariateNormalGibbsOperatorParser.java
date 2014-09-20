package dr.inferencexml.operators;

import dr.inference.distribution.MultivariateDistributionLikelihood;
import dr.inference.operators.MultivariateNormalGibbsOperator;
import dr.math.matrixAlgebra.IllegalDimension;
import dr.xml.*;

/**
 * Created with IntelliJ IDEA.
 * User: max
 * Date: 2/21/14
 * Time: 2:20 PM
 * To change this template use File | Settings | File Templates.
 */
public class MultivariateNormalGibbsOperatorParser extends AbstractXMLObjectParser {
    public static final String MVN_GIBBS_SAMPLER="MultivariateNormalGibbsOperator";
    public static final String PRIOR="prior";
    public static final String LIKELIHOOD="likelihood";
    public static final String WEIGHT="weight";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

       MultivariateDistributionLikelihood prior= (MultivariateDistributionLikelihood) xo.getChild(PRIOR).getChild(MultivariateDistributionLikelihood.class);
       MultivariateDistributionLikelihood likelihood= (MultivariateDistributionLikelihood) xo.getChild(LIKELIHOOD).getChild(MultivariateDistributionLikelihood.class);
//       CompoundParameter data = (CompoundParameter) xo.getChild(CompoundParameter.class);
       String weightTemp= (String) xo.getAttribute(WEIGHT);
       Double weight=Double.parseDouble(weightTemp);

        //TODO check that it gives the right likelihood and the MVN distributions are conformable
//       if (!(prior.getDistribution() instanceof MultivariateNormalDistributionModel)) {
//            throw new XMLParseException("Only a Wishart distribution is conjugate for Gibbs sampling");
//       }
//
//       // Make sure precMatrix is square and dim(precMatrix) = dim(parameter)
//       if (precMatrix.getColumnDimension() != precMatrix.getRowDimension()) {
//           throw new XMLParseException("The variance matrix is not square or of wrong dimension");
//       }


        try {
            return new MultivariateNormalGibbsOperator(likelihood, prior, weight);  //To change body of implemented methods use File | Settings | File Templates.
        } catch (IllegalDimension illegalDimension) {
            illegalDimension.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }
        return null;
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new ElementRule(PRIOR, new XMLSyntaxRule[]{new ElementRule(MultivariateDistributionLikelihood.class)}),
            new ElementRule(LIKELIHOOD, new XMLSyntaxRule[]{new ElementRule(MultivariateDistributionLikelihood.class)}),
//            new ElementRule(CompoundParameter.class),
            AttributeRule.newDoubleRule(WEIGHT),
    };

    @Override
    public String getParserDescription() {
        return "Multivariate Normal Gibbs Sampler";  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public Class getReturnType() {
        return MultivariateNormalGibbsOperator.class;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public String getParserName() {
        return MVN_GIBBS_SAMPLER;  //To change body of implemented methods use File | Settings | File Templates.
    }
}
