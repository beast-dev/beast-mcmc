package dr.evomodel.coalescent.operators;

import dr.inference.model.Parameter;
import dr.inference.operators.*;
import dr.xml.*;

/**
 * A Gibbs operator to update the population size parameters under a Gaussian Markov random field prior
 *
 * @author Vladimir Minin
 * @author Marc Suchard
 * @version $Id: GMRFSkylineGibbsOperator.java,v 1.5 2007/03/20 11:26:49 msuchard Exp $
 */

public class GMRFSkylineGibbsOperator extends SimpleMCMCOperator implements GibbsOperator, CoercableMCMCOperator {

    public static final String GMRF_GIBBS_OPERATOR = "gmrfGibbsOperator";
    public static final String SCALE_FACTOR = "scaleFactor";
    public static final String POPULATION_PARAMETER = "populationSizes";
    public static final String PRECISION_PARAMETER = "precisionParameter";

    private Parameter precisionParameter;
    private Parameter populationSizeParameter;
    private double scaleFactor = 0.5;
    private int mode = CoercableMCMCOperator.DEFAULT;
    private int weight = 1;

    public GMRFSkylineGibbsOperator(Parameter populationSizeParameter, Parameter precisionParameter, double scaleFactor,
                                    int weight, int mode) {
        this.populationSizeParameter = populationSizeParameter;
        this.precisionParameter = precisionParameter;
        this.scaleFactor = scaleFactor;
        this.weight = weight;
        this.mode = mode;
    }

    public double doOperation() throws OperatorFailedException {
        return 0;
    }

    public int getStepCount() {
        return 0;
    }


    //MCMCOperator INTERFACE
    public final String getOperatorName() {
        return GMRF_GIBBS_OPERATOR;
    }

    public double getCoercableParameter() {
        return Math.log(1.0 / scaleFactor - 1.0);
    }

    public void setCoercableParameter(double value) {
        scaleFactor = 1.0 / (Math.exp(value) + 1.0);
    }

    public double getRawParameter() {
        return scaleFactor;
    }

    public int getMode() {
        return mode;
    }

    public double getScaleFactor() {
        return scaleFactor;
    }

    public double getTargetAcceptanceProbability() {
        return 0.234;
    }

    public double getMinimumAcceptanceLevel() {
        return 0.1;
    }

    public double getMaximumAcceptanceLevel() {
        return 0.4;
    }

    public double getMinimumGoodAcceptanceLevel() {
        return 0.20;
    }

    public double getMaximumGoodAcceptanceLevel() {
        return 0.30;
    }

    public int getWeight() {
        return weight;
    }

    public void setWeight(int w) {
        weight = w;
    }

    public final String getPerformanceSuggestion() {

        double prob = MCMCOperator.Utils.getAcceptanceProbability(this);
        double targetProb = getTargetAcceptanceProbability();
        dr.util.NumberFormatter formatter = new dr.util.NumberFormatter(5);
        double sf = OperatorUtils.optimizeScaleFactor(scaleFactor, prob, targetProb);
        if (prob < getMinimumGoodAcceptanceLevel()) {
            return "Try setting scaleFactor to about " + formatter.format(sf);
        } else if (prob > getMaximumGoodAcceptanceLevel()) {
            return "Try setting scaleFactor to about " + formatter.format(sf);
        } else return "";
    }


    public static dr.xml.XMLObjectParser PARSER = new dr.xml.AbstractXMLObjectParser() {

        public String getParserName() {
            return GMRF_GIBBS_OPERATOR;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {


            int mode = CoercableMCMCOperator.DEFAULT;

            if (xo.hasAttribute(AUTO_OPTIMIZE)) {
                if (xo.getBooleanAttribute(AUTO_OPTIMIZE)) {
                    mode = CoercableMCMCOperator.COERCION_ON;
                } else {
                    mode = CoercableMCMCOperator.COERCION_OFF;
                }
            }

            int weight = xo.getIntegerAttribute(WEIGHT);
            double scaleFactor = xo.getDoubleAttribute(SCALE_FACTOR);

            if (scaleFactor <= 0.0 || scaleFactor >= 1.0) {
                throw new XMLParseException("scaleFactor must be between 0.0 and 1.0");
            }

            XMLObject cxo = (XMLObject) xo.getChild(POPULATION_PARAMETER);
            Parameter populationSizeParameter = (Parameter) cxo.getChild(Parameter.class);

            cxo = (XMLObject) xo.getChild(PRECISION_PARAMETER);
            Parameter precisionParameter = (Parameter) cxo.getChild(Parameter.class);

            return new GMRFSkylineGibbsOperator(populationSizeParameter, precisionParameter, scaleFactor, weight, mode);

        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "This element returns a Gibbs operator for the joint distribution of the population sizes and precision parameter.";
        }

        public Class getReturnType() {
            return MCMCOperator.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                AttributeRule.newDoubleRule(SCALE_FACTOR),
                AttributeRule.newIntegerRule(WEIGHT),
                AttributeRule.newBooleanRule(AUTO_OPTIMIZE, true),
                new ElementRule(POPULATION_PARAMETER, new XMLSyntaxRule[]{
                        new ElementRule(Parameter.class)
                }),
                new ElementRule(PRECISION_PARAMETER, new XMLSyntaxRule[]{
                        new ElementRule(Parameter.class)
                })
        };

    };

}
