package dr.inference.model;

import dr.xml.*;

import java.util.ArrayList;

/**
 * A class that returns a log likelihood of a series of statistics.
 * If all the statistics are in increasing (decreasing) order then it returns 0.0 otherwise -INF.
 *
 * @author Marc Suchard
 */

public class MonotonicLikelihood extends BooleanLikelihood {

    public static final String MONOTONIC_LIKELIHOOD = "monotonicLikelihood";
    public static final String DATA = "data";
    public static final String STRICTLY = "strictlyMonotic";
    public static final String ORDER = "order";

    public MonotonicLikelihood(boolean strict, boolean increasing) {
        super();
        this.strict = strict;
        this.increasing = increasing;
    }

    public boolean getBooleanState() {

        double currentValue;
        if (increasing)
            currentValue = Double.NEGATIVE_INFINITY;
        else
            currentValue = Double.POSITIVE_INFINITY;

        for (Statistic statistic : dataList) {
            for (int j = 0; j < statistic.getDimension(); j++) {
                final double newValue = statistic.getStatisticValue(j);
                if (strict) {
                    if( (increasing  && newValue <= currentValue) ||
                        (!increasing && newValue >= currentValue) )
                        return true;
                } else { // not strict
                    if( (increasing  && newValue < currentValue) ||
                        (!increasing && newValue > currentValue) )
                        return true;
                }
                currentValue = newValue;
            }
        }
        return false;
    }

    /**
     * Reads a distribution likelihood from a DOM Document element.
     */
    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return MONOTONIC_LIKELIHOOD;
        }

        public Object parseXMLObject(XMLObject xo) {

            BooleanLikelihood likelihood = new BooleanLikelihood();

            for (int i = 0; i < xo.getChildCount(); i++) {
                if (xo.getChild(i) instanceof BooleanStatistic) {
                    likelihood.addData((BooleanStatistic) xo.getChild(i));
                }
            }

            return likelihood;
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "";
        }

        public Class getReturnType() {
            return MonotonicLikelihood.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                new ElementRule(Statistic.class, 1, Integer.MAX_VALUE)
        };

    };

    private boolean strict;
    private boolean increasing;
}

