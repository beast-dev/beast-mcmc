/*
 * GradientErrorLoggerParser.java
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

package dr.inferencexml.hmc;

import dr.inference.hmc.GradientErrorLogger;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

import static dr.inference.hmc.GradientErrorLogger.Statistic;

/**
 * @author Marc A. Suchard
 * @author Andy Magee
 * @author Andrew Holbrook
 */
public class GradientErrorLoggerParser extends AbstractXMLObjectParser {

    private final static String PARSER_NAME = "gradientErrorLogger";
    private final static String STATISTIC = "statistic";
    private final static String NAME = "name";

    @Override
    public String getParserName() {
        return PARSER_NAME;
    }

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        GradientWrtParameterProvider gradient = (GradientWrtParameterProvider)
                xo.getChild(GradientWrtParameterProvider.class);

        List<Statistic> statistics = parseStatistics(xo);

        return new GradientErrorLogger(gradient, statistics);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(GradientWrtParameterProvider.class),
            new ElementRule(STATISTIC, new XMLSyntaxRule[]{
                    AttributeRule.newStringRule(NAME),
            }, 1, Integer.MAX_VALUE),
    };

    @Override
    public String getParserDescription() {
        return null;
    }

    @Override
    public Class getReturnType() {
        return GradientErrorLogger.class;
    }

    private List<Statistic> parseStatistics(XMLObject xo) throws XMLParseException {
        List<Statistic> statistics = new ArrayList<>();
        for (XMLObject cxo : xo.getAllChildren(STATISTIC)) {
            String name = cxo.getStringAttribute(NAME);
            Statistic statistic = Statistic.parse(name);
            if (statistic == null) {
                throw new XMLParseException("Unknown statistic '" + name + "'");
            }
            statistics.add(statistic);
        }
        return statistics;
    }
}
