/*
 * LogRateSubstitutionModelRestoreTest.java
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

package dr.xml.unittest;

import dr.evomodel.substmodel.DifferentialMassProvider;
import dr.evomodel.substmodel.LogRateSubstitutionModel;
import dr.inference.model.Parameter;
import dr.math.matrixAlgebra.WrappedMatrix;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.AttributeRule;
import dr.xml.ElementRule;
import dr.xml.Reportable;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

/**
 * XML unittest helper for checking that LogRateSubstitutionModel invalidates its
 * derivative-matrix cache after a rejected/restored proposal.
 */
public class LogRateSubstitutionModelRestoreTest implements Reportable {

    private final LogRateSubstitutionModel model;
    private final Parameter parameter;
    private final int dimension;
    private final double proposalOffset;

    public LogRateSubstitutionModelRestoreTest(LogRateSubstitutionModel model,
                                               Parameter parameter,
                                               int dimension,
                                               double proposalOffset) {
        this.model = model;
        this.parameter = parameter;
        this.dimension = dimension;
        this.proposalOffset = proposalOffset;
    }

    @Override
    public String getReport() {

        if (dimension < 0 || dimension >= parameter.getDimension()) {
            throw new IllegalArgumentException("Dimension index out of range");
        }

        final double[] originalValues = parameter.getParameterValues();
        final DifferentialMassProvider.DifferentialWrapper.WrtParameter wrt =
                model.factory(parameter, dimension);

        final double[] reference = copy(model.getInfinitesimalDifferentialMatrix(wrt));

        boolean stored = false;
        double[] proposal = null;
        double[] restored;

        try {
            model.storeModelState();
            stored = true;

            final double[] proposalValues = originalValues.clone();
            proposalValues[dimension] += proposalOffset;
            setParameter(proposalValues);

            proposal = copy(model.getInfinitesimalDifferentialMatrix(wrt));

            model.restoreModelState();
            stored = false;

            restored = copy(model.getInfinitesimalDifferentialMatrix(wrt));

        } finally {
            if (stored && !model.isValidState()) {
                model.restoreModelState();
            }
        }

        StringBuilder sb = new StringBuilder();
        sb.append("LogRateSubstitutionModel restore derivative cache\n");
        sb.append("maxAbsDifference: ").append(maxAbsDifference(reference, restored)).append('\n');
        sb.append("proposalMaxAbsDifference: ").append(maxAbsDifference(reference, proposal)).append('\n');
        appendVector(sb, "reference", reference);
        appendVector(sb, "restored", restored);

        return sb.toString();
    }

    private void setParameter(double[] values) {
        for (int i = 0; i < values.length; ++i) {
            parameter.setParameterValueQuietly(i, values[i]);
        }
        parameter.fireParameterChangedEvent();
    }

    private double[] copy(WrappedMatrix matrix) {
        final double[] values = new double[matrix.getDim()];
        for (int i = 0; i < values.length; ++i) {
            values[i] = matrix.get(i);
        }
        return values;
    }

    private double maxAbsDifference(double[] lhs, double[] rhs) {
        double max = 0.0;
        for (int i = 0; i < lhs.length; ++i) {
            max = Math.max(max, Math.abs(lhs[i] - rhs[i]));
        }
        return max;
    }

    private void appendVector(StringBuilder sb, String label, double[] values) {
        sb.append(label).append(":");
        for (double value : values) {
            sb.append(' ').append(value);
        }
        sb.append('\n');
    }

    private static final String PARSER_NAME = "logRateSubstitutionModelRestoreTest";
    private static final String DIMENSION = "dim";
    private static final String PROPOSAL_OFFSET = "proposalOffset";

    public static final AbstractXMLObjectParser PARSER = new AbstractXMLObjectParser() {

        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            LogRateSubstitutionModel model = (LogRateSubstitutionModel) xo.getChild(LogRateSubstitutionModel.class);
            Parameter parameter = (Parameter) xo.getChild(Parameter.class);
            int dimension = xo.getAttribute(DIMENSION, 0);
            double proposalOffset = xo.getAttribute(PROPOSAL_OFFSET, 1.0);

            return new LogRateSubstitutionModelRestoreTest(model, parameter, dimension, proposalOffset);
        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return new XMLSyntaxRule[]{
                    AttributeRule.newIntegerRule(DIMENSION, true),
                    AttributeRule.newDoubleRule(PROPOSAL_OFFSET, true),
                    new ElementRule(LogRateSubstitutionModel.class),
                    new ElementRule(Parameter.class)
            };
        }

        @Override
        public String getParserDescription() {
            return "Checks LogRateSubstitutionModel derivative-cache behavior after restore.";
        }

        @Override
        public Class getReturnType() {
            return LogRateSubstitutionModelRestoreTest.class;
        }

        @Override
        public String getParserName() {
            return PARSER_NAME;
        }
    };
}
