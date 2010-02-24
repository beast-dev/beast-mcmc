package dr.evomodelxml.coalescent;

import dr.evomodel.coalescent.BMPriorLikelihood;
import dr.evomodel.coalescent.VariableDemographicModel;
import dr.inference.distribution.ParametricDistributionModel;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 */
public class BMPriorLikelihoodParser extends AbstractXMLObjectParser { //TODO Parser is still in a bad state
    public static final String BM = "BrownianMotion";
    public static final String LOG_SPACE = "logUnits";
    //static final String NORMALIZE = "normalize";

    //static final String DATA = "data";
    // static final String TIMES = "times";

    // public static final String MEAN = "mean";
    public static final String SIGMA = "sigma";
    //  public static final String LAMBDA = "lambda";

    public String getParserDescription() {
        return "";
    }

    public Class getReturnType() {
        return BMPriorLikelihood.class;
    }

    public String getParserName() {
        return BM;
    }

    private Parameter getParam(XMLObject xo, String name) throws XMLParseException {
        final XMLObject object = xo.getChild(name);
        // optional
        if (object == null) {
            return null;
        }
        final Object child = object.getChild(0);
        if (child instanceof Parameter) {
            return (Parameter) child;
        }

        double x = object.getDoubleChild(0);
        return new Parameter.Default(x);
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        // Parameter mean = getParam(xo, MEAN);
        Parameter sigma = getParam(xo, SIGMA);
        //  Parameter lambda = getParam(xo, LAMBDA);

        final boolean logSpace = xo.getAttribute(LOG_SPACE, false);
        //   final boolean normalize = xo.getAttribute(NORMALIZE, false);

        VariableDemographicModel m = (VariableDemographicModel) xo.getChild(VariableDemographicModel.class);

        ParametricDistributionModel popMeanPrior = (ParametricDistributionModel) xo.getChild(ParametricDistributionModel.class);
        return new BMPriorLikelihood(sigma, m, logSpace, popMeanPrior);
    }

    public XMLSyntaxRule[] getSyntaxRules() {


        return new XMLSyntaxRule[]{
                AttributeRule.newBooleanRule(LOG_SPACE, true),

                new XORRule(
                        new ElementRule(SIGMA, Double.class),
                        new ElementRule(SIGMA, Parameter.class)
                ),

                // you can't have a XOR (b AND c), yikes
                // make all optional and check in parser
                new ElementRule(VariableDemographicModel.class),

//                    new ElementRule(DATA, new XMLSyntaxRule[]{new ElementRule(Statistic.class)}, true),
//                    new ElementRule(TIMES, new XMLSyntaxRule[]{new ElementRule(Statistic.class)}, true)

        };
    }
}
