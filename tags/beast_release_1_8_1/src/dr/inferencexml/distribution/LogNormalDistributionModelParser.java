package dr.inferencexml.distribution;

import dr.inference.distribution.LogNormalDistributionModel;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * Reads a normal distribution model from a DOM Document element.
 */
public class LogNormalDistributionModelParser extends AbstractXMLObjectParser {

    public static final String LOGNORMAL_DISTRIBUTION_MODEL = "logNormalDistributionModel";
    public static final String MEAN = "mean";
    public static final String STDEV = "stdev";
    public static final String PRECISION = "precision";
    public static final String OFFSET = "offset";
    public static final String MEAN_IN_REAL_SPACE = "meanInRealSpace";
    public static final String STDEV_IN_REAL_SPACE = "stdevInRealSpace";

    public String getParserName() {
        return LOGNORMAL_DISTRIBUTION_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        Parameter meanParam;

        final double offset = xo.getAttribute(OFFSET, 0.0);

        final boolean meanInRealSpace = xo.getAttribute(MEAN_IN_REAL_SPACE, false);
        final boolean stdevInRealSpace = xo.getAttribute(STDEV_IN_REAL_SPACE, false);
        if(!meanInRealSpace && stdevInRealSpace) {
            throw new RuntimeException("Cannot parameterise Lognormal model with M and Stdev");
        }


        {
            final XMLObject cxo = xo.getChild(MEAN);
            if (cxo.getChild(0) instanceof Parameter) {
                meanParam = (Parameter) cxo.getChild(Parameter.class);
            } else {
                meanParam = new Parameter.Default(cxo.getDoubleChild(0));
            }
        }

        {
            final XMLObject cxo = xo.getChild(PRECISION);
            if (cxo != null) {
                Parameter precParam;
                if (cxo.getChild(0) instanceof Parameter) {
                    precParam = (Parameter) cxo.getChild(Parameter.class);
                } else {
                    precParam = new Parameter.Default(cxo.getDoubleChild(0));
                }
                return new LogNormalDistributionModel(meanParam, precParam, offset, meanInRealSpace,stdevInRealSpace, false);
            }
        }
        {
            final XMLObject cxo = xo.getChild(STDEV);
            Parameter stdevParam;
            if (cxo.getChild(0) instanceof Parameter) {
                stdevParam = (Parameter) cxo.getChild(Parameter.class);
            } else {
                stdevParam = new Parameter.Default(cxo.getDoubleChild(0));
            }

            return new LogNormalDistributionModel(meanParam, stdevParam, offset, meanInRealSpace, stdevInRealSpace);
        }
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newBooleanRule(MEAN_IN_REAL_SPACE, true),
            AttributeRule.newBooleanRule(STDEV_IN_REAL_SPACE, true),
            AttributeRule.newDoubleRule(OFFSET, true),
            new ElementRule(MEAN,
                    new XMLSyntaxRule[]{
                            new XORRule(
                                    new ElementRule(Parameter.class),
                                    new ElementRule(Double.class)
                            )}
            ),
            new XORRule(
                    new ElementRule(STDEV,
                            new XMLSyntaxRule[]{
                                    new XORRule(
                                            new ElementRule(Parameter.class),
                                            new ElementRule(Double.class)
                                    )}
                    ),
                    new ElementRule(PRECISION,
                            new XMLSyntaxRule[]{
                                    new XORRule(
                                            new ElementRule(Parameter.class),
                                            new ElementRule(Double.class)
                                    )}

                    ))
    };

    public String getParserDescription() {
        return "Describes a normal distribution with a given mean and standard deviation " +
                "that can be used in a distributionLikelihood element";
    }

    public Class getReturnType() {
        return LogNormalDistributionModel.class;
    }
}
