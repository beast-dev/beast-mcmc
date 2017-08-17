/*
 * MonotonicStatistic.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.inference.model;

import dr.xml.*;

import java.util.ArrayList;

/**
 * A class that returns a log likelihood of a series of statistics.
 * If all the statistics are in increasing (decreasing) order then it returns 0.0 otherwise -INF.
 *
 * @author Marc Suchard
 */

public class MonotonicStatistic extends Statistic.Abstract implements BooleanStatistic {

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

    /**
     * @return boolean result of test.
     */
    public double getStatisticValue(int dim) {
        return getBoolean(dim) ? 1.0 : 0.0;
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

