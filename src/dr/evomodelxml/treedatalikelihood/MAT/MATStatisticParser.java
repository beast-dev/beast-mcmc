/*
 * MATStatisticParser.java
 *
 * Copyright © 2026, the BEAST Development Team.
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
 */

package dr.evomodelxml.treedatalikelihood.MAT;

import dr.evomodel.treedatalikelihood.DataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.MAT.MATDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.inference.model.Statistic;
import dr.xml.*;
import dr.evomodel.treedatalikelihood.MAT.MATStatistic;

/**
 * Parses a statistic reporting the number of mutations and ambiguities in a MAT.
 */
public class MATStatisticParser extends AbstractXMLObjectParser {

    public static final String MAT_STATISTIC = "matStatistic";

    @Override
    public String getParserName() {
        return MAT_STATISTIC;
    }

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        String name = xo.getAttribute(Statistic.NAME, xo.getId());
        return new MATStatistic(name, MutationHistoryResamplerParser.getDelegate(xo));
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[] {
                AttributeRule.newStringRule(Statistic.NAME, true),
                new ElementRule(TreeDataLikelihood.class)
        };
    }

    @Override
    public String getParserDescription() {
        return "Reports the number of mutations and ambiguities (per site) in a mutation annotated tree.";
    }

    @Override
    @SuppressWarnings("rawtypes")
    public Class getReturnType() {
        return MATStatistic.class;
    }
}
