package dr.inferencexml.distribution;

import dr.inference.distribution.InverseGaussianDistributionModel;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * Reads an inverse gaussian distribution model from a DOM Document element.
 */
public class InverseGaussianDistributionModelParser extends AbstractXMLObjectParser {

    public static final String INVERSEGAUSSIAN_DISTRIBUTION_MODEL = "inverseGaussianDistributionModel";
    public static final String MEAN = "mean";
    public static final String STDEV = "stdev";
    public static final String SHAPE = "shape";
    public static final String OFFSET = "offset";

    public String getParserName() {
        return INVERSEGAUSSIAN_DISTRIBUTION_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Parameter meanParam;

        double offset = xo.getAttribute(OFFSET, 0.0);

        XMLObject cxo = xo.getChild(MEAN);
        if (cxo.getChild(0) instanceof Parameter) {
            meanParam = (Parameter) cxo.getChild(Parameter.class);
        } else {
            meanParam = new Parameter.Default(cxo.getDoubleChild(0));
        }

        if(xo.hasChildNamed(STDEV) && xo.hasChildNamed(SHAPE)) {
            throw new RuntimeException("XML has both standard deviation and shape for Inverse Gaussian distribution");
        }
        else if(xo.hasChildNamed(STDEV)) {
            Parameter stdevParam;
            cxo = xo.getChild(STDEV);
            if (cxo.getChild(0) instanceof Parameter) {
                stdevParam = (Parameter) cxo.getChild(Parameter.class);
            } else {
                stdevParam = new Parameter.Default(cxo.getDoubleChild(0));
            }
            return new InverseGaussianDistributionModel(meanParam, stdevParam, offset, false);
        }
        else if(xo.hasChildNamed(SHAPE)) {
            Parameter shapeParam;
            cxo = xo.getChild(SHAPE);
            if (cxo.getChild(0) instanceof Parameter) {
                shapeParam = (Parameter) cxo.getChild(Parameter.class);
            } else {
                shapeParam = new Parameter.Default(cxo.getDoubleChild(0));
            }
            return new InverseGaussianDistributionModel(meanParam, shapeParam, offset, true);
        }
        else {
            throw new RuntimeException("XML has neither standard deviation nor shape for Inverse Gaussian distribution");
        }
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newDoubleRule(OFFSET, true),
            new ElementRule(MEAN,
                    new XMLSyntaxRule[]{
                            new XORRule(
                                    new ElementRule(Parameter.class),
                                    new ElementRule(Double.class)
                            )}
            , false),
            new ElementRule(STDEV,
                    new XMLSyntaxRule[]{
                            new XORRule(
                                    new ElementRule(Parameter.class),
                                    new ElementRule(Double.class)
                            )}
            , true),

            new ElementRule(SHAPE,
                new XMLSyntaxRule[]{
                        new XORRule(
                                new ElementRule(Parameter.class),
                                new ElementRule(Double.class)
                        )}
            , true)
    };

    public String getParserDescription() {
        return "Describes a inverse gaussian distribution with a given mean and shape (or standard deviation) " +
                "that can be used in a distributionLikelihood element";
    }

    public Class getReturnType() {
        return InverseGaussianDistributionModel.class;
    }

}
