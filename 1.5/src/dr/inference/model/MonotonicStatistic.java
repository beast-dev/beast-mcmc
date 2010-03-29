package dr.inference.model;

import dr.xml.*;

import java.util.ArrayList;

/**
 * A class that returns a log likelihood of a series of statistics.
 * If all the statistics are in increasing (decreasing) order then it returns 0.0 otherwise -INF.
 *
 * @author Marc Suchard
 */

public class MonotonicStatistic extends BooleanStatistic {

    public static final String MONOTONIC_STATISTIC = "monotonicStatistic";
    public static final String STRICTLY = "strictlyMonotic";
    public static final String ORDER = "order";

    public MonotonicStatistic(boolean strict, boolean increasing) {
        super(MONOTONIC_STATISTIC);
        this.strict = strict;
        this.increasing = increasing;
    }

    public int getDimension() {
        return 1;
    }

    public void addStatistic(Statistic stat) {
        dataList.add(stat);
    }

    public boolean getBoolean(int dim) {

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
                        return false;
                } else { // not strict
                    if( (increasing  && newValue < currentValue) ||
                        (!increasing && newValue > currentValue) )
                        return false;
                }
                currentValue = newValue;
            }
        }
        return true;
    }

    /**
     * Reads a distribution likelihood from a DOM Document element.
     */
    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return MONOTONIC_STATISTIC;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            boolean increasing = true;
            String order = xo.getAttribute(ORDER, "increasing");
            if (order.compareToIgnoreCase("decreasing") == 0)
                increasing = false;

            boolean strictly = xo.getAttribute(STRICTLY, false);

            MonotonicStatistic monotonicStatistic = new MonotonicStatistic(strictly, increasing);

            for(int i=0; i<xo.getChildCount(); i++) {
                if (xo.getChild(i) instanceof Statistic) {
                    monotonicStatistic.addStatistic((Statistic)xo.getChild(i));
                }
            }

            return monotonicStatistic;
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "";
        }

        public Class getReturnType() {
            return MonotonicStatistic.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                AttributeRule.newStringRule(ORDER, true),
                AttributeRule.newBooleanRule(STRICTLY,true),
                new ElementRule(Statistic.class, 1, Integer.MAX_VALUE)
        };

    };

    private boolean strict;
    private boolean increasing;
    protected ArrayList<Statistic> dataList = new ArrayList<Statistic>();

}

