/*
 * ReflectiveHMCStatisticsParser.java
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

package dr.inferencexml.operators.hmc;

import dr.inference.operators.hmc.ReflectiveHMCStatistics;
import dr.inference.operators.hmc.ReflectiveHamiltonianMonteCarloOperator;
import dr.xml.*;

/**
 * Parses &lt;reflectiveHMCStatistics&gt;, the log columns of a reflective HMC operator's counters.
 */
public class ReflectiveHMCStatisticsParser extends AbstractXMLObjectParser {

    public static final String STATISTICS = "reflectiveHMCStatistics";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        ReflectiveHamiltonianMonteCarloOperator operator =
                (ReflectiveHamiltonianMonteCarloOperator) xo.getChild(ReflectiveHamiltonianMonteCarloOperator.class);
        return new ReflectiveHMCStatistics(operator, xo.hasId() ? xo.getId() : "rhmc");
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{new ElementRule(ReflectiveHamiltonianMonteCarloOperator.class)};
    }

    @Override
    public String getParserDescription() {
        return "Log columns with a reflective HMC operator's event and rejection counters";
    }

    @Override
    public Class getReturnType() {
        return ReflectiveHMCStatistics.class;
    }

    @Override
    public String getParserName() {
        return STATISTICS;
    }
}
