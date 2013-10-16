package dr.inferencexml.distribution;

import dr.inference.distribution.NormalDistributionModel;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * Reads a normal distribution model from a DOM Document element.
 */
public class NormalDistributionModelParser extends AbstractXMLObjectParser {

    public static final String NORMAL_DISTRIBUTION_MODEL = "normalDistributionModel";
    public static final String MEAN = "mean";
    public static final String STDEV = "stdev";
    public static final String PREC = "precision";

    public String getParserName() {
        return NORMAL_DISTRIBUTION_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Parameter meanParam;
        Parameter stdevParam;
        Parameter precParam;

        XMLObject cxo = xo.getChild(MEAN);
        if (cxo.getChild(0) instanceof Parameter) {
            meanParam = (Parameter) cxo.getChild(Parameter.class);
        } else {
            meanParam = new Parameter.Default(cxo.getDoubleChild(0));
        }

        if (xo.getChild(STDEV) != null) {

            cxo = xo.getChild(STDEV);
            if (cxo.getChild(0) instanceof Parameter) {
                stdevParam = (Parameter) cxo.getChild(Parameter.class);
            } else {
                stdevParam = new Parameter.Default(cxo.getDoubleChild(0));
            }

            return new NormalDistributionModel(meanParam, stdevParam);
        }

        cxo = xo.getChild(PREC);
        if (cxo.getChild(0) instanceof Parameter) {
            precParam = (Parameter) cxo.getChild(Parameter.class);
        } else {
            precParam = new Parameter.Default(cxo.getDoubleChild(0));
        }
        return new NormalDistributionModel(meanParam, precParam, true);
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
            )
    };

    public String getParserDescription() {
        return "Describes a normal distribution with a given mean and standard deviation " +
                "that can be used in a distributionLikelihood element";
    }

    public Class getReturnType() {
        return NormalDistributionModel.class;
    }

}
