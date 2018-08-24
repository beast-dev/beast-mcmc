/*
 * CompoundMatrixParameter.java
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

package dr.inference.model;

import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Marc A. Suchard
 */
public class CompoundMatrixParameter extends MatrixParameter {

    public CompoundMatrixParameter(String name, List<MatrixParameterInterface> matrices) {
        super(name, compoundMatrices(matrices));
    }

    private static Parameter[] compoundMatrices(List<MatrixParameterInterface> matrices) {
        int length = 0;
        for (MatrixParameterInterface matrix : matrices) {
            length += matrix.getColumnDimension();
        }

        Parameter[] parameters = new Parameter[length];
        int index = 0;
        for (MatrixParameterInterface matrix : matrices) {
            for (int i = 0; i < matrix.getColumnDimension(); ++i) {
                parameters[index] = matrix.getParameter(i);
                ++index;
            }
        }
        return parameters;
    }

    public final static String COMPOUND_MATRIX_PARAMETER = "compoundMatrixParameter";

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return COMPOUND_MATRIX_PARAMETER;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            List<MatrixParameterInterface> matrices = new ArrayList<MatrixParameterInterface>();

            for (int i = 0; i < xo.getChildCount(); ++i) {
                matrices.add((MatrixParameterInterface) xo.getChild(i));
            }

            final String name = xo.hasId() ? xo.getId() : null;

            return new CompoundMatrixParameter(name, matrices);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "A compound matrix parameter constructed from its component parameters.";
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                new ElementRule(MatrixParameter.class, 1, Integer.MAX_VALUE),
        };

        public Class getReturnType() {
            return CompoundMatrixParameter.class;
        }
    };
}
