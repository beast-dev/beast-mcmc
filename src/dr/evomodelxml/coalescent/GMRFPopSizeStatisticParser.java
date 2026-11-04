/*
 * GMRFPopSizeStatisticParser.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
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
 *
 */

package dr.evomodelxml.coalescent;

import dr.evomodel.coalescent.GMRFPopSizeStatistic;
import dr.evomodel.coalescent.GMRFSkyrideLikelihood;
import dr.xml.*;

/**
 *
 */
public class GMRFPopSizeStatisticParser extends AbstractXMLObjectParser {

    public final static String TIMES = "time";
    public final static String FROM = "from";
    public final static String TO = "to";
    public final static String NUMBER_OF_INTERVALS = "number";
    public final static String GMRF_POP_SIZE_STATISTIC = "gmrfPopSizeStatistic";

    public String getParserName() {
        return GMRF_POP_SIZE_STATISTIC;
    }
    
    public String getParserDescription() {
        return "The pop sizes at the given times";
    }

    public Class getReturnType() {
        return GMRFPopSizeStatistic.class;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        double[] times;

        if (xo.hasAttribute(FROM)) {
            double from = xo.getDoubleAttribute(FROM);
            double to = xo.getDoubleAttribute(TO);
            int number = xo.getIntegerAttribute(NUMBER_OF_INTERVALS);

            double length = (to - from) / number;

            times = new double[number + 1];

            for (int i = 0; i < times.length; i++) {
                times[i] = from + i * length;
            }
        } else {
            times = xo.getDoubleArrayAttribute(TIMES);
        }

        GMRFSkyrideLikelihood gsl = (GMRFSkyrideLikelihood) xo.getChild(GMRFSkyrideLikelihood.class);

        return new GMRFPopSizeStatistic(times, gsl);
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                AttributeRule.newDoubleRule(FROM, true),
                AttributeRule.newDoubleRule(TO, true),
                AttributeRule.newIntegerRule(NUMBER_OF_INTERVALS, true),
                AttributeRule.newDoubleArrayRule(TIMES, true),

                new ElementRule(GMRFSkyrideLikelihood.class),
        };
    }
    
}
