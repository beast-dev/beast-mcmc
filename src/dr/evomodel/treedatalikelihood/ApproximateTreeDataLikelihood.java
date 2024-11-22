/*
 * TreeDataLikelihood.java
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

package dr.evomodel.treedatalikelihood;

import dr.evolution.tree.Tree;
import dr.evomodel.treedatalikelihood.discrete.MaximizerWrtParameter;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.hmc.HessianWrtParameterProvider;
import dr.inference.hmc.JointGradient;
import dr.inference.model.*;
import dr.inference.operators.hmc.NumericalHessianFromGradient;
import dr.math.MultivariateFunction;
import dr.math.matrixAlgebra.WrappedVector;
import dr.util.Transform;
import dr.xml.*;
import dr.math.NumericalDerivative;

import java.util.List;

import static dr.math.matrixAlgebra.ReadableVector.Utils.setParameter;

/**
 * @author Alexander Fisher
 */

public class ApproximateTreeDataLikelihood extends AbstractModelLikelihood {
    private double marginalLikelihood;
    private MaximizerWrtParameter maximizer;
    private Parameter parameter;
    private Likelihood likelihood;
    private boolean likelihoodKnown = false;
    private final HessianWrtParameterProvider hessianWrtParameterProvider;


    // begin parser stuff
    public static final String APPROXIMATE_LIKELIHOOD = "approximateTreeDataLikelihood";

    // end parser stuff
    public ApproximateTreeDataLikelihood(MaximizerWrtParameter maximizer) {
        super(APPROXIMATE_LIKELIHOOD);

        this.maximizer = maximizer;
        this.likelihood = maximizer.getLikelihood();
        final GradientWrtParameterProvider gradient = maximizer.getGradient();
        this.parameter = gradient.getParameter();
        this.marginalLikelihoodConst = (parameter.getDimension() - 1) * Math.log(2 * Math.PI);
        // todo: get Numerical Hessian.
        if (maximizer.getTransform() != null) {
            this.hessianWrtParameterProvider = constructHessian();
        } else if (isGradientProvidingHessian(gradient)) {
            this.hessianWrtParameterProvider = (HessianWrtParameterProvider) gradient;
        } else {
            this.hessianWrtParameterProvider = new NumericalHessianFromGradient(gradient);
        }
        updateParameterMAP();
        updateMarginalLikelihood();
        addVariable(parameter);
    }

    private boolean isGradientProvidingHessian(GradientWrtParameterProvider gradient) {
        boolean isHessianProvider = false;
        if (gradient instanceof HessianWrtParameterProvider) {
            if (gradient instanceof JointGradient) {
                JointGradient jointGradient = (JointGradient) gradient;
                boolean isNotHessianProvider = false;
                for (GradientWrtParameterProvider gradientWrtParameterProvider : jointGradient.getDerivativeList()) {
                    if (!(gradientWrtParameterProvider instanceof HessianWrtParameterProvider)) {
                        isNotHessianProvider = true;
                    }
                }
                isHessianProvider = !isNotHessianProvider;
            } else {
                isHessianProvider = true;
            }
        }
        return isHessianProvider;
    }

    private HessianWrtParameterProvider constructHessian() {
        GradientWrtParameterProvider gradientWrtParameterProvider = new GradientWrtParameterProvider() {

            private TransformedMultivariateParameter transformedParameter = new TransformedMultivariateParameter(parameter, (Transform.MultivariableTransform) maximizer.getTransform());

            @Override
            public Likelihood getLikelihood() {
                throw new RuntimeException("should not be called");
            }

            @Override
            public Parameter getParameter() {
                return transformedParameter;
            }

            @Override
            public int getDimension() {
                return transformedParameter.getDimension();
            }

            @Override
            public double[] getGradientLogDensity() {
                double[] untransformedGradient = maximizer.getGradient().getGradientLogDensity();
                return maximizer.getTransform().updateGradientLogDensity(untransformedGradient, parameter.getParameterValues(), 0, parameter.getDimension());
            }
        };

        return new NumericalHessianFromGradient(gradientWrtParameterProvider);
    }

    private void updateMarginalLikelihood() {
        double[] diagonalHessian = hessianWrtParameterProvider.getDiagonalHessianLogDensity();
        double logDiagonalDeterminant = 0;
        for (int i = 0; i < parameter.getDimension(); i++) {
            logDiagonalDeterminant += Math.log(Math.abs(diagonalHessian[i]));
        }
        // 2pi^{-k/2} * det(Sigma)^{-1/2} * likelihood(map) * prior(map)
        this.marginalLikelihood = marginalLikelihoodConst + 0.5 * logDiagonalDeterminant + likelihood.getLogLikelihood()
         + (maximizer.getTransform() == null ? 0 : maximizer.getTransform().logJacobian(parameter.getParameterValues(), 0, parameter.getDimension()));
        likelihoodKnown = true;
    }

    private final double marginalLikelihoodConst;

    private void updateParameterMAP() {
        maximizer.maximize();
    }


    // **************************************************************
    // XMLObjectParser
    // **************************************************************

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return APPROXIMATE_LIKELIHOOD;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            MaximizerWrtParameter maximizer =
                    (MaximizerWrtParameter) xo.getChild(MaximizerWrtParameter.class);

            return new ApproximateTreeDataLikelihood(maximizer);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        @Override
        public String getParserDescription() {
            return "Approximates the marginal likelihood of the data given the tree using Laplace approximation";
        }

        @SuppressWarnings("rawtypes")
        @Override
        public Class getReturnType() {
            return ApproximateTreeDataLikelihood.class;
        }

        private final XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                new ElementRule(MaximizerWrtParameter.class)
        };
    };

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {

    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {

    }

    @Override
    protected void storeState() {

    }

    @Override
    protected void restoreState() {

    }

    @Override
    protected void acceptState() {

    }

    @Override
    public Model getModel() {
        return null;
    }

    @Override
    public double getLogLikelihood() {
        if (!likelihoodKnown) {
            updateParameterMAP();
            updateMarginalLikelihood();
        }
        return marginalLikelihood;
    }

    @Override
    public void makeDirty() {
        likelihoodKnown = false;
    }
}