/*
 * LKJTransformParser.java
 *
 * Copyright (c) 2002-2018 Alexei Drummond, Andrew Rambaut and Marc Suchard
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


package dr.util;

import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Xiang Ji
 * @author Marc Suchard
 */
public class MultipleTransformParser extends AbstractXMLObjectParser {

    private static final String NAME = "multipleTransform";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        List<Transform> array = new ArrayList<Transform>();
        for (Object xoc : xo.getChildren()) {
            if (xoc instanceof Transform) {
                array.add((Transform) xoc);
            } else if (xoc instanceof Transform.ParsedTransform) {
                array.add(((Transform.ParsedTransform) xoc).transform);
            }
        }
        return new Transform.MultipleTransform(array);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(Transform.class, 1, Integer.MAX_VALUE),
    };

    @Override
    public String getParserDescription() {
        return null;
    }

    @Override
    public Class getReturnType() {
        return Transform.MultipleTransform.class;
    }

    @Override
    public String getParserName() {
        return NAME;
    }
}
