/*
 * NodeHeightTransformTest.java
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

package dr.evomodel.treedatalikelihood.discrete;

import dr.inference.model.Parameter;
import dr.math.MultivariateFunction;
import dr.math.NumericalDerivative;
import dr.xml.*;

/**
 * @author Marc A. Suchard
 * @author Xiang Ji
 */
public class NodeHeightTransformTest implements Reportable{

    private final NodeHeightTransform nodeHeightTransform;
    private final NodeHeightGradientForDiscreteTrait nodeHeightGradient;
    private final Parameter ratios;

    public NodeHeightTransformTest(NodeHeightTransform nodeHeightTransform,
                                   NodeHeightGradientForDiscreteTrait nodeHeightGradient,
                                   Parameter ratios) {

        this.nodeHeightTransform = nodeHeightTransform;
        this.nodeHeightGradient = nodeHeightGradient;
        this.ratios = ratios;

    }

    @Override
    public String getReport() {
        String message = nodeHeightGradient.getReport();
        double[] gradient = nodeHeightGradient.getGradientLogDensity();
        double[] updatedGradient = nodeHeightTransform.updateGradientLogDensity(gradient, nodeHeightGradient.getParameter().getParameterValues(), 0, gradient.length);
        double[] numericGradient = NumericalDerivative.diagonalHessian(numeric1, ratios.getParameterValues());

        StringBuilder sb = new StringBuilder();
        sb.append("\nGradient Peeling: ").append(new dr.math.matrixAlgebra.Vector(updatedGradient));
        sb.append("\nGradient Numeric: ").append(new dr.math.matrixAlgebra.Vector(numericGradient));
        return message + sb.toString();
    }

    protected MultivariateFunction numeric1 = new MultivariateFunction() {
        @Override
        public double evaluate(double[] argument) {

            for (int i = 0; i < argument.length; ++i) {
                ratios.setParameterValueQuietly(i, argument[i]);
            }
            ratios.fireParameterChangedEvent();

//            treeDataLikelihood.makeDirty();
            return nodeHeightGradient.getLikelihood().getLogLikelihood();
        }

        @Override
        public int getNumArguments() {
            return ratios.getDimension();
        }

        @Override
        public double getLowerBound(int n) {
            return 0;
        }

        @Override
        public double getUpperBound(int n) {
            return 1.0;
        }
    };

    private static final String NODE_HEIGHT_TRANSFORM_TEST = "nodeHeightTransformTest";

    public static AbstractXMLObjectParser PARSER = new AbstractXMLObjectParser() {
        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            NodeHeightTransform transform = (NodeHeightTransform) xo.getChild(NodeHeightTransform.class);
            NodeHeightGradientForDiscreteTrait gradient = (NodeHeightGradientForDiscreteTrait) xo.getChild(NodeHeightGradientForDiscreteTrait.class);
            Parameter ratios = (Parameter) xo.getChild(Parameter.class);
            return new NodeHeightTransformTest(transform, gradient, ratios);
        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return new XMLSyntaxRule[]{
                    new ElementRule(NodeHeightTransform.class),
                    new ElementRule(NodeHeightGradientForDiscreteTrait.class),
                    new ElementRule(Parameter.class),
            };
        }

        @Override
        public String getParserDescription() {
            return null;
        }

        @Override
        public Class getReturnType() {
            return NodeHeightTransformTest.class;
        }

        @Override
        public String getParserName() {
            return NODE_HEIGHT_TRANSFORM_TEST;
        }
    };
}
