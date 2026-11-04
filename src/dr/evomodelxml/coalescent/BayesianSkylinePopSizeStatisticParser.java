/*
 * BayesianSkylinePopSizeStatisticParser.java
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

import dr.evomodel.coalescent.BayesianSkylineLikelihood;
import dr.evomodel.coalescent.BayesianSkylinePopSizeStatistic;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

/**
 */
public class BayesianSkylinePopSizeStatisticParser extends AbstractXMLObjectParser {
    public static final String TIME = "time";
    public static final String BAYESIAN_SKYLINE_POP_SIZE_STATISTIC = "generalizedSkylinePopSizeStatistic";


    public String getParserDescription() {
        return "The pop sizes at the given times";
    }

    public Class getReturnType() {
        return BayesianSkylinePopSizeStatistic.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return null;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        double time = xo.getDoubleAttribute(TIME);

        BayesianSkylineLikelihood bsl =
                (BayesianSkylineLikelihood) xo.getChild(BayesianSkylineLikelihood.class);

        return new BayesianSkylinePopSizeStatistic(time, bsl);
    }

    public String getParserName() {
        return BAYESIAN_SKYLINE_POP_SIZE_STATISTIC;
    }

}
