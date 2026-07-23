/*
 * ElementWiseMatrixMultiplicationParser.java
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


import dr.inference.model.ElementWiseMatrixMultiplicationParameter;
import dr.inference.model.MatrixParameter;
import dr.inference.model.MatrixParameterInterface;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

/**
 * Created by max on 11/30/15.
 */
public class ElementWiseMatrixMultiplicationParser extends AbstractXMLObjectParser {
    public final static String ELEMENT_WISE_MATRIX_MULTIPLICATION_PARAMETER="elementWiseMatrixMultiplicationParameter";
    public final static String NAME="name";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        final String name = xo.hasId() ? xo.getId() : null;
        MatrixParameterInterface[] matList=new MatrixParameterInterface[xo.getChildCount()];
        for (int i = 0; i <xo.getChildCount(); i++) {
            matList[i]=(MatrixParameterInterface) xo.getChild(i);
        }

        return new ElementWiseMatrixMultiplicationParameter(name, matList);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[0];
    }

    @Override
    public String getParserDescription() {
        return "Returns element wise matrix multiplication of a series of matrices";
    }

    @Override
    public Class getReturnType() {
        return ElementWiseMatrixMultiplicationParameter.class;
    }

    @Override
    public String getParserName() {
        return ELEMENT_WISE_MATRIX_MULTIPLICATION_PARAMETER;
    }
}
