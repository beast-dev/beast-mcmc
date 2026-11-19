/*
 * DifferentiableSubstitutionModelDerivativeReport.java
 *
 * Copyright (c) 2002-2024 the BEAST Development Team
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
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 *
 */

package dr.xml.unittest;

import dr.evomodel.substmodel.DifferentiableSubstitutionModel;
import dr.evomodel.substmodel.DifferentialMassProvider;
import dr.inference.model.Parameter;
import dr.inference.model.Statistic;
import dr.math.matrixAlgebra.WrappedMatrix;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.AttributeRule;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

/**
 * Numeric probe for a differentiable substitution model's infinitesimal
 * differential matrix with respect to one parameter dimension.
 */
public class DifferentiableSubstitutionModelDerivativeReport extends Statistic.Abstract {

    private final DifferentiableSubstitutionModel model;
    private final DifferentialMassProvider.DifferentialWrapper.WrtParameter wrt;
    private final int stateCount;

    public DifferentiableSubstitutionModelDerivativeReport(DifferentiableSubstitutionModel model,
                                                           Parameter parameter,
                                                           int dimension) {
        super("differentiableSubstitutionModelDerivativeReport");
        this.model = model;
        this.wrt = model.factory(parameter, dimension);
        this.stateCount = model.getDataType().getStateCount();
    }

    @Override
    public int getDimension() {
        return stateCount * stateCount;
    }

    @Override
    public double getStatisticValue(int dim) {
        if (dim < 0 || dim >= getDimension()) {
            throw new IllegalArgumentException("Derivative matrix dimension index out of range");
        }
        return getMatrix().get(dim);
    }

    @Override
    public String getDimensionName(int dim) {
        return getStatisticName() + Integer.toString(dim + 1);
    }

    @Override
    public String getReport() {
        StringBuilder sb = new StringBuilder();
        sb.append("differentialMatrix:");
        WrappedMatrix matrix = getMatrix();
        for (int i = 0; i < getDimension(); ++i) {
            sb.append(' ').append(matrix.get(i));
        }
        sb.append('\n');
        return sb.toString();
    }

    private WrappedMatrix getMatrix() {
        return model.getInfinitesimalDifferentialMatrix(wrt);
    }

    private static final String PARSER_NAME = "differentiableSubstitutionModelDerivativeReport";
    private static final String DIMENSION = "dim";

    public static final AbstractXMLObjectParser PARSER = new AbstractXMLObjectParser() {

        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            DifferentiableSubstitutionModel model =
                    (DifferentiableSubstitutionModel) xo.getChild(DifferentiableSubstitutionModel.class);
            Parameter parameter = (Parameter) xo.getChild(Parameter.class);
            int dimension = xo.getIntegerAttribute(DIMENSION);

            if (dimension < 0 || dimension >= parameter.getDimension()) {
                throw new XMLParseException(PARSER_NAME + " dim is outside the parameter dimension");
            }

            return new DifferentiableSubstitutionModelDerivativeReport(model, parameter, dimension);
        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return new XMLSyntaxRule[]{
                    AttributeRule.newIntegerRule(DIMENSION),
                    new ElementRule(DifferentiableSubstitutionModel.class),
                    new ElementRule(Parameter.class)
            };
        }

        @Override
        public String getParserDescription() {
            return "Reports the infinitesimal differential matrix for a differentiable substitution model.";
        }

        @Override
        public Class getReturnType() {
            return DifferentiableSubstitutionModelDerivativeReport.class;
        }

        @Override
        public String getParserName() {
            return PARSER_NAME;
        }
    };
}
