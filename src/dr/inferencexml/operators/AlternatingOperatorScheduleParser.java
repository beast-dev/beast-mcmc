/*
 * AlternatingOperatorScheduleParser.java
 *
 * Copyright (c) 2002-2016 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.inferencexml.operators;

import dr.inference.operators.AlternatingOperatorSchedule;
import dr.inference.operators.OperatorSchedule;
import dr.xml.*;

/**
 *
 */
public class AlternatingOperatorScheduleParser extends AbstractXMLObjectParser {

    public static final String ALTERNATING_OPERATORS = "alternatingOperators";
    public static final String SCHEDULE = "schedule";
    public static final String COUNT = "count";

    public String getParserName() {
        return ALTERNATING_OPERATORS;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        AlternatingOperatorSchedule schedule = new AlternatingOperatorSchedule();

        for (int i = 0; i < xo.getChildCount(); i++) {
            XMLObject cxo = (XMLObject)xo.getChild(i);
            long count = cxo.getLongIntegerAttribute(COUNT);
            schedule.addOperatorSchedule((OperatorSchedule)cxo.getChild(OperatorSchedule.class), count);
        }

        return schedule;
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(SCHEDULE, new XMLSyntaxRule[] {
                    AttributeRule.newLongIntegerRule(COUNT),
                    new ElementRule(OperatorSchedule.class)
            }, 1, Integer.MAX_VALUE)
    };

    public String getParserDescription() {
        return "An operator scheduler that alternates between nested schedulers";
    }

    public Class getReturnType() {
        return AlternatingOperatorSchedule.class;
    }
    
}
