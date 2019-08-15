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
import dr.util.Transform;
import dr.xml.*;

/**
 * @author Marc A. Suchard
 * @author Xiang Ji
 */
public class NodeHeightTransformTest implements Reportable{

    private final NodeHeightTransform nodeHeightTransform;
    private final NodeHeightGradientForDiscreteTrait nodeHeightGradient;
    private final Parameter ratios;
    private final Transform.ComposeMultivariable logitRatioTransform;

    public NodeHeightTransformTest(NodeHeightTransform nodeHeightTransform,
                                   NodeHeightGradientForDiscreteTrait nodeHeightGradient,
                                   Parameter ratios) {

        this.nodeHeightTransform = nodeHeightTransform;
        this.nodeHeightGradient = nodeHeightGradient;
        this.ratios = ratios;
        Transform.Array logitTransforms = new Transform.Array(Transform.LOGIT, ratios.getDimension(), null);
        this.logitRatioTransform = new Transform.ComposeMultivariable(logitTransforms, nodeHeightTransform);
    }

    @Override
    public String getReport() {
        String message = nodeHeightGradient.getReport();
        double[] gradient = nodeHeightGradient.getGradientLogDensity();
        double[] updatedUnweightedGradient = nodeHeightTransform.updateGradientUnWeightedLogDensity(gradient, nodeHeightGradient.getParameter().getParameterValues(), 0, gradient.length);
        double[] numericUnweightedGradient = NumericalDerivative.gradient(numericUnweighted, ratios.getParameterValues());
        double[] updatedWeightedGradient = nodeHeightTransform.updateGradientLogDensity(gradient, nodeHeightGradient.getParameter().getParameterValues(), 0, gradient.length);
        double[] numericWeightedGradient = NumericalDerivative.gradient(numericWeighted, ratios.getParameterValues());

        StringBuilder sb = new StringBuilder();
        sb.append("\nGradient wrt Unweighted LogLikelihood:");
        sb.append("\nPeeling: ").append(new dr.math.matrixAlgebra.Vector(updatedUnweightedGradient));
        sb.append("\nNumeric: ").append(new dr.math.matrixAlgebra.Vector(numericUnweightedGradient));
        sb.append("\nGradient wrt Weighted LogLikelihood:");
        sb.append("\nPeeling: ").append(new dr.math.matrixAlgebra.Vector(updatedWeightedGradient));
        sb.append("\nNumeric: ").append(new dr.math.matrixAlgebra.Vector(numericWeightedGradient));


        double[] updatedMultipleWeightedGradient = logitRatioTransform.updateGradientLogDensity(gradient, nodeHeightTransform.getNodeHeights().getParameterValues(), 0, nodeHeightTransform.getNodeHeights().getDimension());
        double[] numericMultipleWeightedGradient = NumericalDerivative.gradient(numericMultipleWeighted,
                logitRatioTransform.transform(nodeHeightTransform.getNodeHeights().getParameterValues(), 0, nodeHeightTransform.getNodeHeights().getDimension()));

        sb.append("\nGradient wrt Multiple Weighted LogLikelihood:");
        sb.append("\nPeeling: ").append(new dr.math.matrixAlgebra.Vector(updatedMultipleWeightedGradient));
        sb.append("\nNumeric: ").append(new dr.math.matrixAlgebra.Vector(numericMultipleWeightedGradient));

        return message + sb.toString();
    }

    protected MultivariateFunction numericUnweighted = new MultivariateFunction() {
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

    protected MultivariateFunction numericWeighted = new MultivariateFunction() {
        @Override
        public double evaluate(double[] argument) {

            for (int i = 0; i < argument.length; ++i) {
                ratios.setParameterValueQuietly(i, argument[i]);
            }
            ratios.fireParameterChangedEvent();

//            treeDataLikelihood.makeDirty();
            return nodeHeightGradient.getLikelihood().getLogLikelihood() - nodeHeightTransform.getLogJacobian(argument);
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

    protected MultivariateFunction numericMultipleWeighted = new MultivariateFunction() {
        @Override
        public double evaluate(double[] argument) {

            double[] inverseValues = logitRatioTransform.inverse(argument, 0, argument.length);
            Parameter nodeHeights = nodeHeightTransform.getNodeHeights();

            for (int i = 0; i < inverseValues.length; ++i) {
                nodeHeights.setParameterValueQuietly(i, inverseValues[i]);
            }
            nodeHeightTransform.getNodeHeights().fireParameterChangedEvent();

            nodeHeightGradient.getLikelihood().makeDirty();
            final double result = nodeHeightGradient.getLikelihood().getLogLikelihood() - logitRatioTransform.getLogJacobian(inverseValues, 0, argument.length);
            return result;
        }

        @Override
        public int getNumArguments() {
            return nodeHeightTransform.getNodeHeights().getDimension();
        }

        @Override
        public double getLowerBound(int n) {
            return Double.NEGATIVE_INFINITY;
        }

        @Override
        public double getUpperBound(int n) {
            return Double.POSITIVE_INFINITY;
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
