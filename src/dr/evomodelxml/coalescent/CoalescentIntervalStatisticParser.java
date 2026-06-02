/*
 * CoalescentIntervalStatisticParser.java
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

import dr.evomodel.coalescent.CoalescentIntervalProvider;
import dr.evomodel.coalescent.CoalescentIntervalStatistic;
import dr.xml.*;


public class CoalescentIntervalStatisticParser extends AbstractXMLObjectParser {

    public static final String COALESCENT_INTERVAL_STATISTIC = "coalescentIntervalStatistic";

    public String getParserDescription() {
        return "";
    }

    public Class getReturnType() {
        return CoalescentIntervalStatistic.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                new ElementRule(CoalescentIntervalProvider.class),
        };
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        CoalescentIntervalProvider coalescent = (CoalescentIntervalProvider) xo.getChild(CoalescentIntervalProvider.class);
        return new CoalescentIntervalStatistic(coalescent);
    }

    public String getParserName() {
        return COALESCENT_INTERVAL_STATISTIC;
    }

}
