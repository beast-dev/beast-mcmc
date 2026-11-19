/*
 * LineageCountStatisticParser.java
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

import dr.evolution.coalescent.IntervalList;
import dr.evomodel.coalescent.AbstractCoalescentLikelihood;
import dr.evomodel.coalescent.LineageCountStatistic;
import dr.evomodel.coalescent.OldAbstractCoalescentLikelihood;
import dr.xml.*;

/**
 */
public class LineageCountStatisticParser extends AbstractXMLObjectParser {

    public static final String LINEAGE_COUNT_STATISTIC = "lineageCountStatistic";


    public String getParserDescription() {
        return "";
    }

    public Class getReturnType() {
        return LineageCountStatistic.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                new XORRule(new XMLSyntaxRule[]{
                        new ElementRule(IntervalList.class),
                        new ElementRule(AbstractCoalescentLikelihood.class),
                        new ElementRule(OldAbstractCoalescentLikelihood.class)
                })
        };
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        OldAbstractCoalescentLikelihood old =
                (OldAbstractCoalescentLikelihood) xo.getChild(OldAbstractCoalescentLikelihood.class);
        if (old != null) {
            return new LineageCountStatistic(old);
        }
        AbstractCoalescentLikelihood likelihood =
                (AbstractCoalescentLikelihood) xo.getChild(AbstractCoalescentLikelihood.class);
        if (likelihood != null) {
            return new LineageCountStatistic(likelihood.getIntervalList());
        }
        return new LineageCountStatistic((IntervalList) xo.getChild(IntervalList.class));
    }

    public String getParserName() {
        return LINEAGE_COUNT_STATISTIC;
    }

}
