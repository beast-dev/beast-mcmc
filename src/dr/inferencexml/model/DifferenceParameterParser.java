/*
 * DifferenceParameterParser.java
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

package dr.inferencexml.model;

import dr.inference.model.DifferenceParameter;
import dr.inference.model.Parameter;
import dr.inference.model.ProductParameter;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

/**
 */
public class DifferenceParameterParser extends AbstractXMLObjectParser {

    public static final String DIFFERENCE_PARAMETER = "differenceParameter";

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        List<Parameter> paramList = new ArrayList<Parameter>();
        int dim = -1;
        for (int i = 0; i < xo.getChildCount(); ++i) {
            Parameter parameter = (Parameter) xo.getChild(i);
            if (dim == -1) {
                dim = parameter.getDimension();
            } else {
                if (parameter.getDimension() != dim) {
                    throw new XMLParseException("All parameters in difference '" + xo.getId() + "' must be the same length");
                }
            }
            paramList.add(parameter);
        }

        return new DifferenceParameter(paramList);
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(Parameter.class,1,Integer.MAX_VALUE),
    };

    public String getParserDescription() {
        return "A element-wise difference of parameters.";
    }

    public Class getReturnType() {
        return Parameter.class;
    }

    public String getParserName() {
        return DIFFERENCE_PARAMETER;
    }
}
