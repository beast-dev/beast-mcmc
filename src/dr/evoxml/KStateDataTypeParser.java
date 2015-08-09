/*
 * KStateDataTypeParser.java
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

package dr.evoxml;

import dr.evolution.datatype.DataType;
import dr.evolution.datatype.GeneralDataType;
import dr.xml.*;

import java.util.ArrayList;
import java.util.Collection;

/**
 * Package: KStateDataTypeParser
 * Description:
 * <p/>
 * <p/>
 * Created by
 *
 * @author Alexander V. Alekseyenko (alexander.alekseyenko@gmail.com)
 *         Date: Oct 13, 2009
 *         Time: 3:42:14 PM
 */
public class KStateDataTypeParser extends AbstractXMLObjectParser {
    public static final String K_STATE_DATATYPE = "kStateType";
    public static final String STATE_COUNT = "stateCount";
    public static final String START_WITH = "startWith";

    //public static XMLObjectParser PARSER=new KStateDataTypeParser();

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        int k = xo.getIntegerAttribute(STATE_COUNT);
        int sw = 0;
        if (xo.hasAttribute(START_WITH))
            sw = xo.getIntegerAttribute(START_WITH);
        Collection<String> states = new ArrayList<String>();

        System.err.println(states.toArray().toString());

        for (int i = sw; i < k + sw; ++i) {
            states.add(Integer.toString(i));
        }
        System.err.println("Constructing " + k + "-state datatype");
        return new GeneralDataType(states);
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newIntegerRule(STATE_COUNT, false),
            AttributeRule.newIntegerRule(START_WITH, true)};

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    public String getParserDescription() {
        return "Parser for k-state model.";
    }

    public Class getReturnType() {
        return DataType.class;
    }

    public String getParserName() {
        return K_STATE_DATATYPE;
    }

}
