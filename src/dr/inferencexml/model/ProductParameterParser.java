/*
 * ProductParameterParser.java
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

package dr.inferencexml.model;

import dr.inference.model.Parameter;
import dr.inference.model.ProductParameter;
import dr.inference.model.ScaledParameter;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

/**
 *
 */
public class ProductParameterParser extends AbstractXMLObjectParser {

    public static final String PRODUCT_PARAMETER = "productParameter";
    public static final String SCALE = "scale";

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        List<Parameter> paramList = new ArrayList<Parameter>();
        int dim = -1;
        for (int i = 0; i < xo.getChildCount(); ++i) {
            if (xo.getChild(i) instanceof Parameter) {
                Parameter parameter = (Parameter) xo.getChild(i);
                if (dim == -1) {
                    dim = parameter.getDimension();
                } else {
                    if (parameter.getDimension() != dim) {
                        throw new XMLParseException("All parameters in product '" + xo.getId() +
                                "' must be the same length");
                    }
                }
                paramList.add(parameter);
            }
        }

        ProductParameter prodParam = new ProductParameter(paramList);

        if (xo.hasChildNamed(SCALE)) {
            Parameter scaleParam = (Parameter) xo.getChild(SCALE).getChild(Parameter.class);

            if (scaleParam.getDimension() != 1) {
                throw new XMLParseException("The scale parameter must be one-dimensional.");
            }

            Parameter vecParam;

            if (paramList.size() == 1) { //Don't pass a product parameter if you don't need to
                vecParam = paramList.get(0);
            } else {
                vecParam = prodParam;
            }


            return new ScaledParameter(scaleParam, vecParam);
        } else {
            return prodParam;
        }

    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(Parameter.class, 1, Integer.MAX_VALUE),
            new ElementRule(SCALE, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }, true)
    };

    public String getParserDescription() {
        return "A element-wise product of parameters.";
    }

    public Class getReturnType() {
        return Parameter.class;
    }

    public String getParserName() {
        return PRODUCT_PARAMETER;
    }
}
