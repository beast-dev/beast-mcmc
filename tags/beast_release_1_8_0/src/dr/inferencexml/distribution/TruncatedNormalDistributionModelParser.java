package dr.inferencexml.distribution;

import dr.inference.distribution.TruncatedNormalDistributionModel;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * Reads a normal distribution model from a DOM Document element.
 */
public class TruncatedNormalDistributionModelParser extends AbstractXMLObjectParser {

    public static final String TRUNCATED_NORMAL_DISTRIBUTION_MODEL = "truncatedNormalDistributionModel";
    public static final String MEAN = "mean";
    public static final String STDEV = "stdev";
    public static final String MINIMUM = "minimum";
    public static final String MAXIMUM = "maximum";
    public static final String PREC = "precision";

    public String getParserName() {
        return TRUNCATED_NORMAL_DISTRIBUTION_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Parameter meanParam;
        Parameter stdevParam;
        Parameter minParam;
        Parameter maxParam;
        Parameter precParam;

        XMLObject cxo = xo.getChild(MEAN);
        if (cxo.getChild(0) instanceof Parameter) {
            meanParam = (Parameter) cxo.getChild(Parameter.class);
        } else {
            meanParam = new Parameter.Default(cxo.getDoubleChild(0));
        }

        if(xo.hasChildNamed(MINIMUM)){
            cxo = xo.getChild(MINIMUM);
            if(cxo.getChild(0) instanceof Parameter) {
                minParam = (Parameter) cxo.getChild(Parameter.class);
            } else {
                minParam = new Parameter.Default(cxo.getDoubleChild(0));
            }
        } else {
            minParam = new Parameter.Default(Double.NEGATIVE_INFINITY);
        }

        if(xo.hasChildNamed(MAXIMUM)){
            cxo = xo.getChild(MAXIMUM);
            if(cxo.getChild(0) instanceof Parameter) {
                maxParam = (Parameter) cxo.getChild(Parameter.class);
            } else {
                maxParam = new Parameter.Default(cxo.getDoubleChild(0));
            }
        } else {
            maxParam = new Parameter.Default(Double.POSITIVE_INFINITY);
        }

        if (xo.getChild(STDEV) != null) {

            cxo = xo.getChild(STDEV);
            if (cxo.getChild(0) instanceof Parameter) {
                stdevParam = (Parameter) cxo.getChild(Parameter.class);
            } else {
                stdevParam = new Parameter.Default(cxo.getDoubleChild(0));
            }

            return new TruncatedNormalDistributionModel(meanParam, stdevParam, minParam, maxParam);
        }

        cxo = xo.getChild(PREC);
        if (cxo.getChild(0) instanceof Parameter) {
            precParam = (Parameter) cxo.getChild(Parameter.class);
        } else {
            precParam = new Parameter.Default(cxo.getDoubleChild(0));
        }
        return new TruncatedNormalDistributionModel(meanParam, precParam, minParam, maxParam, true);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
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
                    new ElementRule(PREC,
                            new XMLSyntaxRule[]{
                                    new XORRule(
                                            new ElementRule(Parameter.class),
                                            new ElementRule(Double.class)
                                    )}
                    )
            ),
            new OrRule(
                    new ElementRule(MINIMUM,
                            new XMLSyntaxRule[]{
                                    new XORRule(
                                            new ElementRule(Parameter.class),
                                            new ElementRule(Double.class)
                                    )}
                    ),
                    new ElementRule(MAXIMUM,
                            new XMLSyntaxRule[]{
                                    new XORRule(
                                            new ElementRule(Parameter.class),
                                            new ElementRule(Double.class)
                                    )
                            }
                    )
            )
    };

    public String getParserDescription() {
        return "Describes a truncated normal distribution with a given mean, standard deviation and minimum or" +
                "maximum (or both) values that can be used in a distributionLikelihood element";
    }

    public Class getReturnType() {
        return TruncatedNormalDistributionModel.class;
    }

}
