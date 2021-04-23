/*
 * RateStatisticParser.java
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

package dr.evomodelxml.epidemiology;

import dr.evolution.tree.Tree;
import dr.evolution.util.Taxa;
import dr.evolution.util.Units;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.epidemiology.EpidemiologyStatistic;
import dr.evomodel.tree.RateStatistic;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Parameter;
import dr.inference.model.Statistic;
import dr.xml.*;

import java.util.Locale;

/**
 */
public class EpidemiologyStatisticParser extends AbstractXMLObjectParser {

    public static final String EPIDEMIOLOGY_STATISTIC = "epidemiologyStatistic";
    public static final String GROWTH_RATE = "growthRate";
    public static final String TIME_UNITS = "timeUnits";
    public static final String DOUBLING_TIME = "doublingTime";
    public static final String R0 = "R0";
    public static final String SERIAL_INTERVAL = "serialInterval";
    public static final String MEAN = "mean";
    public static final String STDEV = "stdev";

    public String getParserName() {
        return EPIDEMIOLOGY_STATISTIC;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        final String name = xo.getAttribute(Statistic.NAME, xo.getId());

        XMLObject cxo = xo.getChild(GROWTH_RATE);
        final Parameter growthRate = (Parameter)cxo.getChild(Parameter.class);
        Units.Type timeUnits = Units.Type.valueOf(cxo.getAttribute(TIME_UNITS, "years").toUpperCase());

        if (getParserName().equals(DOUBLING_TIME)) {
            return new EpidemiologyStatistic(name, growthRate, timeUnits);
        }

        // R0
        if (xo.hasChildNamed(SERIAL_INTERVAL)) {
            cxo = xo.getChild(SERIAL_INTERVAL);
            double mean = cxo.getAttribute(MEAN, 0.0);
            double stdev = cxo.getAttribute(STDEV, 0.0);

            return new EpidemiologyStatistic(name, growthRate, timeUnits, mean, stdev);
        } else {
            throw new XMLParseException("R0 statistic requires a serialInterval element");
        }
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    @Override
    public String[] getParserNames() {
        return new String[] {
                DOUBLING_TIME, R0
        };
    }

    public String getParserDescription() {
        return "A statistic that transforms the doubling time or R0 from a growth rate";
    }

    public Class getReturnType() {
        return EpidemiologyStatistic.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(GROWTH_RATE, new XMLSyntaxRule[] {
                    new StringAttributeRule(TIME_UNITS, "time units for the growth rate (default: years)", true),
                    new ElementRule(Parameter.class)
            },  "The growth rate parameter"),
            new ElementRule(SERIAL_INTERVAL, new XMLSyntaxRule[] {
                    AttributeRule.newDoubleRule(MEAN, false, "the mean of the gamma distribution modelling serial interval"),
                    AttributeRule.newDoubleRule(STDEV, false, "the stdev of the gamma distribution modelling serial interval"),
            },  "The serial interval distribution", true)
    };

}
