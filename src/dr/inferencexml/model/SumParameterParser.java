/*
 * ProductParameterParser.java
 *
 * Copyright (c) 2002-2012 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.inferencexml.model;

import dr.inference.model.Statistic;
import dr.inference.model.SumParameter;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

/**
 */
public class SumParameterParser extends AbstractXMLObjectParser {

    public static final String SUM_ALL = "sumAll";
    public static final String SUM_PARAMETER = "sumParameter";

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        List<Statistic> statisticList = new ArrayList<>();
        int dim = -1;
        for (int i = 0; i < xo.getChildCount(); ++i) {
            Statistic s = (Statistic) xo.getChild(i);
            if (dim == -1) {
                dim = s.getDimension();
            } else {
                if (s.getDimension() != dim) {
                    throw new XMLParseException("All statistics/parameters in sum '" + xo.getId() + "' must be the same length");
                }
            }
            statisticList.add(s);
        }

        boolean sumAll = statisticList.size() == 1;
        if (xo.hasAttribute(SUM_ALL)) {
            sumAll = xo.getBooleanAttribute(SUM_ALL);
        }

        if (sumAll && statisticList.size() > 1) {
            throw new XMLParseException("To sum all the elements, only one parameter should be given");
        }
        if (!sumAll && statisticList.size() < 2) {
            throw new XMLParseException("For an element-wise sum, more than one parameter should be given");
        }

        return new SumParameter(statisticList);
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(Statistic.class,1,Integer.MAX_VALUE),
    };

    public String getParserDescription() {
        return "A element-wise sum of statistics or parameters.";
    }

    public Class getReturnType() {
        return Parameter.class;
    }

    public String getParserName() {
        return SUM_PARAMETER;
    }
}
