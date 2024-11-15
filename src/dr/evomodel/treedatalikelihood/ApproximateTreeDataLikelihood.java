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

import dr.evomodel.treedatalikelihood.discrete.MaximizerWrtParameter;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.util.Transform;
import dr.xml.*;
import dr.math.NumericalDerivative;

/**
 * @author Alexander Fisher
 */

public class ApproximateTreeDataLikelihood {
    private double marginalLikelihood;
    private double[] parameterMAP;
    private MaximizerWrtParameter maximizer;
    private double[] numericalHessian;
    private Parameter parameter;


    // begin parser stuff
    public static final String APPROXIMATE_LIKELIHOOD = "approximateTreeDataLikelihood";

    // end parser stuff
    public ApproximateTreeDataLikelihood(MaximizerWrtParameter maximizer) {

        this.maximizer = maximizer;
        this.parameter = maximizer.getGradient().getParameter();
        this.numericalHessian = new double[parameter.getDimension()];
        // todo: get Numerical Hessian.
//    NumericalDerivative.getNumericalHessian();
        updateParameterMAP();
        updateMarginalLikelihood();

    }

    private void updateMarginalLikelihood() {
        double diagonalDeterminant = 1;
        for (int i = 0; i < parameter.getDimension(); i++) {
            diagonalDeterminant *= numericalHessian[i];
        }
        // 2pi^{-k/2} * det(Sigma)^{-1/2} * likelihood(map) * prior(map)
        // todo: eval posterior(map)
        // todo: log likelihood
        this.marginalLikelihood = 2 / (Math.pow(Math.PI, -1 * parameter.getDimension() / 2) * Math.sqrt(diagonalDeterminant));
    }

    private void updateParameterMAP() {
        this.parameterMAP = maximizer.getMinimumPoint(true);
    }

    public double getMarginalLikelihood() {
        return marginalLikelihood;
    }

    public double[] getParameterMAP() {
        return parameterMAP;
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
}